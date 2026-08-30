#!/usr/bin/env python3
import csv
import sys

def correct_dob(dob, line_number, log_entries):
    """ Corrects the DOB format based on the specified rules and logs changes. """
    original_dob = dob.strip()  # Remove any leading/trailing spaces

    # Rule: Replace "00000000" (or any all-zero input) with "01010001"
    if dob == "00000000" or set(dob) == {"0"}:
        corrected_dob = "01010001"
        log_entries.append(f"Line {line_number}: Corrected DOB '{original_dob}' to '{corrected_dob}' (all zeros case)")
        return corrected_dob

    # Rule: If DOB is 7 characters long (missing leading zero in day)
    if len(dob) == 7:
        corrected_dob = dob[:2] + "0" + dob[2:]  # Insert zero before the day part
        log_entries.append(f"Line {line_number}: Corrected DOB '{original_dob}' to '{corrected_dob}' (missing leading zero in day)")
        return corrected_dob

    # Rule: If DOB is 6 characters or less, assume YYYY is correct and prepend '0101'
    if len(dob) <= 6:
        corrected_dob = "0101" + dob[-4:]  # Keeping only the year part
        log_entries.append(f"Line {line_number}: Corrected DOB '{original_dob}' to '{corrected_dob}' (short DOB case)")
        return corrected_dob

    # If DOB is already 8 characters, return as-is
    return dob

def correct_parent_id(value, line_number, column_name, log_entries):
    """ Converts sire/dam ID from '-1' to '0' and logs the correction. """
    if value.strip() == "-1":
        log_entries.append(f"Line {line_number}: Corrected {column_name} '-1' to '0'")
        return "0"
    return value

def process_csv(input_file, output_file, log_file):
    """ Reads the input CSV, corrects errors, writes the corrected output CSV, and logs changes. """
    log_entries = []

    with open(input_file, newline='', encoding='utf-8') as infile, \
         open(output_file, 'w', newline='', encoding='utf-8') as outfile, \
         open(log_file, 'w', encoding='utf-8') as logf:
        
        reader = csv.reader(infile)
        writer = csv.writer(outfile)

        headers = next(reader)
        writer.writerow(headers)  # Write the header row

        for line_number, row in enumerate(reader, start=2):  # Line 1 is header, so data starts at line 2
            if len(row) < 4:
                log_entries.append(f"Line {line_number}: Skipped due to insufficient columns")
                continue  # Skip lines that don't have enough columns

            row[1] = correct_parent_id(row[1], line_number, "SireId", log_entries)
            row[2] = correct_parent_id(row[2], line_number, "DamId", log_entries)
            row[-1] = correct_dob(row[-1], line_number, log_entries)

            writer.writerow(row)

        # Write log file
        for entry in log_entries:
            logf.write(entry + "\n")

    print(f"Processing complete. Output written to {output_file}. Log written to {log_file}.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input.csv output.csv log.txt")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_csv = sys.argv[2]
    log_txt = sys.argv[3]

    process_csv(input_csv, output_csv, log_txt)
