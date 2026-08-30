"""
NAME: pyp_db.py
VERSION: originally 2.0.2; current package version is PyPedal.__version__
AUTHOR: John B. Cole, PhD
LICENSE: LGPL
Modified for the PyPedal 4.0 Python 3 release by Geert Degezelle, 2025-2026.
See CHANGELOG.md. SPDX-License-Identifier: LGPL-2.1-or-later.

Database backend rewritten from ADOdb to sqlite3 during the Python 3 /
PyPedal 4 modernization. Original author remains John B. Cole.
"""

import logging
import sqlite3
from . import pyp_io, pyp_nrm, pyp_utils

def connect_to_database(pedobj):
    """
    Opens a connection to a user-specified SQLite 3 database.
    """
    try:
        conn = sqlite3.connect(pedobj.kw['database_file'])
        logging.info(f"Connected to SQLite database {pedobj.kw['database_file']}.")
        return conn
    except sqlite3.Error as e:
        logging.error(f"Unable to connect to SQLite database {pedobj.kw['database_file']}: {e}")
        return None


def create_pedigree_table(pedobj, conn=None, drop=False):
    """
    Creates a new pedigree table in the database. Drops existing table if drop=True.
    """
    created_conn = False
    if conn is None:
        conn = connect_to_database(pedobj)
        created_conn = True

    if not conn:
        logging.error("No valid database connection.")
        return False

    if does_table_exist(pedobj, conn):
        if drop:
            table_drop_rows(pedobj, conn)
            logging.warning(f"Dropped existing table rows in {pedobj.kw['database_table']}.")
        else:
            logging.info(f"Table {pedobj.kw['database_table']} already exists. No action taken.")
            if created_conn:
                conn.close()
            return False

    try:
        sql = f'''
        CREATE TABLE {pedobj.kw['database_table']} (
            animalID INTEGER PRIMARY KEY,
            animalName TEXT,
            sireID INTEGER,
            sireName TEXT,
            damID INTEGER,
            damName TEXT,
            generation REAL,
            infGeneration REAL,
            birthyear INTEGER,
            sex TEXT,
            coi REAL,
            founder TEXT,
            ancestor TEXT,
            originalID TEXT,
            renumberedID INTEGER,
            pedigreeComp REAL,
            breed TEXT,
            age REAL,
            alive TEXT,
            num_sons INTEGER,
            num_daus INTEGER,
            num_unk INTEGER,
            herd INTEGER,
            originalHerd TEXT,
            gencoeff REAL,
            alleles TEXT,
            userField TEXT
        );
        '''
        conn.execute(sql)
        conn.commit()
        logging.info(f"Table {pedobj.kw['database_table']} created successfully.")
        return True
    except sqlite3.Error as e:
        logging.error(f"Failed to create table {pedobj.kw['database_table']}: {e}")
        return False
    finally:
        if created_conn:
            conn.close()


def does_table_exist(pedobj, conn=None):
    """
    Checks if the specified table exists in the database.
    """
    created_conn = False
    if conn is None:
        conn = connect_to_database(pedobj)
        created_conn = True

    if not conn:
        logging.error("No valid database connection.")
        return False

    try:
        sql = f"SELECT name FROM sqlite_master WHERE type='table' AND name='{pedobj.kw['database_table']}';"
        cursor = conn.execute(sql)
        result = cursor.fetchone()
        return result is not None
    except sqlite3.Error as e:
        logging.error(f"Error checking table existence: {e}")
        return False
    finally:
        if created_conn:
            conn.close()


def table_count_rows(pedobj, conn=None):
    """
    Returns the number of rows in a table.
    """
    created_conn = False
    if conn is None:
        conn = connect_to_database(pedobj)
        created_conn = True

    if not conn:
        logging.error("No valid database connection.")
        return 0

    try:
        if does_table_exist(pedobj, conn):
            sql = f"SELECT COUNT(*) FROM {pedobj.kw['database_table']};"
            cursor = conn.execute(sql)
            count = cursor.fetchone()[0]
            logging.info(f"Table {pedobj.kw['database_table']} has {count} rows.")
            return count
        else:
            logging.warning(f"Table {pedobj.kw['database_table']} does not exist.")
            return 0
    except sqlite3.Error as e:
        logging.error(f"Error counting rows: {e}")
        return 0
    finally:
        if created_conn:
            conn.close()


def table_drop_rows(pedobj, conn=None):
    """
    Deletes all rows from an existing table.
    """
    created_conn = False
    if conn is None:
        conn = connect_to_database(pedobj)
        created_conn = True

    if not conn:
        logging.error("No valid database connection.")
        return False

    try:
        if does_table_exist(pedobj, conn):
            sql = f"DELETE FROM {pedobj.kw['database_table']};"
            conn.execute(sql)
            conn.commit()
            logging.info(f"All rows deleted from table {pedobj.kw['database_table']}.")
            return True
        else:
            logging.warning(f"Table {pedobj.kw['database_table']} does not exist.")
            return False
    except sqlite3.Error as e:
        logging.error(f"Error deleting rows: {e}")
        return False
    finally:
        if created_conn:
            conn.close()


def populate_pedigree_table(pedobj, conn=None):
    """
    Populates the pedigree table with data from a PyPedal pedigree object.
    """
    created_conn = False
    if conn is None:
        conn = connect_to_database(pedobj)
        created_conn = True

    if not conn:
        logging.error("No valid database connection.")
        return False

    try:
        if not does_table_exist(pedobj, conn):
            if not create_pedigree_table(pedobj, conn):
                logging.error("Failed to create table for population.")
                return False

        for record in pedobj.pedigree:
            alleles = '__'.join(record.alleles)
            sql = f'''
            INSERT INTO {pedobj.kw['database_table']} (
                animalID, animalName, sireID, sireName, damID, damName, generation, 
                infGeneration, birthyear, sex, coi, founder, ancestor, originalID, 
                renumberedID, pedigreeComp, breed, age, alive, num_sons, num_daus, 
                num_unk, herd, originalHerd, gencoeff, alleles, userField
            ) VALUES (
                ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?
            );
            '''
            conn.execute(sql, (
                record.animalID, record.name, record.sireID, record.sireName,
                record.damID, record.damName, record.gen, record.igen,
                record.by, record.sex, record.fa, record.founder, record.ancestor,
                record.originalID, record.renumberedID, record.pedcomp,
                record.breed, record.age, record.alive, len(record.sons),
                len(record.daus), len(record.unks), record.herd,
                record.originalHerd, record.gencoeff, alleles, str(record.userField)
            ))
        conn.commit()
        logging.info("Table populated successfully.")
        return True
    except sqlite3.Error as e:
        logging.error(f"Error populating table: {e}")
        return False
    finally:
        if created_conn:
            conn.close()


def delete_table(pedobj, conn=None):
    """Drop the pedigree table if it exists."""
    created_conn = False
    if conn is None:
        conn = connect_to_database(pedobj)
        created_conn = True

    if not conn:
        logging.error("No valid database connection.")
        return False

    try:
        if does_table_exist(pedobj, conn):
            conn.execute(f"DROP TABLE {pedobj.kw['database_table']};")
            conn.commit()
            logging.info(f"Dropped table {pedobj.kw['database_table']}.")
            return True
        logging.warning(f"Table {pedobj.kw['database_table']} does not exist.")
        return False
    except sqlite3.Error as e:
        logging.error(f"Error dropping table: {e}")
        return False
    finally:
        if created_conn:
            conn.close()


# Historical camelCase names used by NewPedigree and examples.
connectToDatabase = connect_to_database
doesTableExist = does_table_exist
deleteTable = delete_table
tableDropRows = table_drop_rows
createPedigreeTable = create_pedigree_table
populatePedigreeTable = populate_pedigree_table
