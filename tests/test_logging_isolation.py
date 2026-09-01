"""PyPedal logging must not hijack the process root logger."""
import logging
import os
import subprocess
import sys
import textwrap

from _pedhelpers import REPO


def _run(script, extra_env=None):
    env = os.environ.copy()
    env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
    if extra_env:
        env.update(extra_env)
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=REPO,
        env=env,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise AssertionError(
            "subprocess failed\nSTDOUT:\n%s\nSTDERR:\n%s"
            % (result.stdout, result.stderr)
        )
    return result


def test_import_does_not_configure_the_root_logger():
    script = textwrap.dedent(
        """
        import logging
        root = logging.getLogger()
        before_handlers = list(root.handlers)
        before_level = root.level
        import PyPedal
        root = logging.getLogger()
        assert list(root.handlers) == before_handlers, list(root.handlers)
        assert root.level == before_level
        names = [type(h).__name__ for h in logging.getLogger("PyPedal").handlers]
        assert "NullHandler" in names
        print("ok")
        """
    )
    result = _run(script)
    assert "ok" in result.stdout
    assert "INFO:root:" not in result.stderr


def test_pedigree_logfile_is_written_with_a_pypedal_record(tmp_path):
    pedfile = tmp_path / "t.ped"
    pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    script = textwrap.dedent(
        f"""
        import os
        os.chdir({str(tmp_path)!r})
        from PyPedal.pyp_newclasses import load_pedigree
        ped = load_pedigree(options={{
            "pedfile": {str(pedfile)!r},
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
            "renumber": True,
        }})
        logfile = ped.kw["logfile"]
        assert os.path.exists(logfile), logfile
        text = open(logfile, encoding="utf-8").read()
        assert "Logfile" in text and "instantiated" in text, text
        print(logfile)
        """
    )
    result = _run(script)
    logfile = result.stdout.strip().splitlines()[-1]
    assert os.path.exists(logfile)


def test_quiet_does_not_emit_root_info_on_stderr(tmp_path):
    pedfile = tmp_path / "t.ped"
    pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    script = textwrap.dedent(
        f"""
        from PyPedal.pyp_newclasses import load_pedigree
        load_pedigree(options={{
            "pedfile": {str(pedfile)!r},
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
            "renumber": True,
        }})
        print("ok")
        """
    )
    result = _run(script)
    assert "ok" in result.stdout
    assert "INFO:root:" not in result.stderr
    assert "INFO:PyPedal" not in result.stderr
    assert "SNP data is not available" not in result.stderr


def test_host_root_handlers_are_not_removed(tmp_path):
    pedfile = tmp_path / "t.ped"
    pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    script = textwrap.dedent(
        f"""
        import logging
        host = logging.StreamHandler()
        host._host_owned = True
        logging.getLogger().addHandler(host)
        import PyPedal
        from PyPedal.pyp_newclasses import load_pedigree
        load_pedigree(options={{
            "pedfile": {str(pedfile)!r},
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
            "renumber": True,
        }})
        root_handlers = logging.getLogger().handlers
        assert host in root_handlers, root_handlers
        assert getattr(host, "_host_owned", False)
        print("ok")
        """
    )
    result = _run(script)
    assert "ok" in result.stdout


def test_second_pedigree_replaces_the_owned_logfile_handler(tmp_path):
    first = tmp_path / "a.ped"
    second = tmp_path / "b.ped"
    first.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    second.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    script = textwrap.dedent(
        f"""
        import logging
        import os
        from PyPedal.pyp_newclasses import (
            load_pedigree, _PYPEDAL_OWNED_HANDLER, PYPEDAL_LOGGER_NAME,
        )
        ped_a = load_pedigree(options={{
            "pedfile": {str(first)!r},
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
            "renumber": True,
        }})
        log_a = ped_a.kw["logfile"]
        ped_b = load_pedigree(options={{
            "pedfile": {str(second)!r},
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
            "renumber": True,
        }})
        log_b = ped_b.kw["logfile"]
        pkg = logging.getLogger(PYPEDAL_LOGGER_NAME)
        owned = [h for h in pkg.handlers if getattr(h, _PYPEDAL_OWNED_HANDLER, False)]
        assert len(owned) == 1, owned
        assert os.path.abspath(owned[0].baseFilename) == os.path.abspath(log_b)
        text_a = open(log_a, encoding="utf-8").read()
        text_b = open(log_b, encoding="utf-8").read()
        assert log_a not in text_b
        assert "instantiated" in text_a
        assert "instantiated" in text_b
        print("ok")
        """
    )
    result = _run(script)
    assert "ok" in result.stdout


def test_close_owned_handlers_leaves_host_and_root_alone(tmp_path):
    """Test utility must close only PyPedal-owned logfile handlers."""
    from _pedhelpers import close_owned_pypedal_log_handlers
    from PyPedal.pyp_newclasses import (
        PYPEDAL_LOGGER_NAME,
        _PYPEDAL_OWNED_HANDLER,
        load_pedigree,
    )

    pedfile = tmp_path / "t.ped"
    pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    root = logging.getLogger()
    package = logging.getLogger(PYPEDAL_LOGGER_NAME)
    host_root = logging.StreamHandler()
    host_pkg = logging.StreamHandler()
    root.addHandler(host_root)
    package.addHandler(host_pkg)
    try:
        load_pedigree(
            options={
                "pedfile": str(pedfile),
                "pedformat": "asd",
                "messages": "quiet",
                "pedigree_summary": 0,
                "renumber": True,
            }
        )
        owned_before = [
            handler
            for handler in package.handlers
            if getattr(handler, _PYPEDAL_OWNED_HANDLER, False)
        ]
        assert owned_before
        logfile = owned_before[0].baseFilename
        close_owned_pypedal_log_handlers()
        owned_after = [
            handler
            for handler in package.handlers
            if getattr(handler, _PYPEDAL_OWNED_HANDLER, False)
        ]
        assert owned_after == []
        assert host_pkg in package.handlers
        assert host_root in root.handlers
        os.remove(logfile)
    finally:
        root.removeHandler(host_root)
        package.removeHandler(host_pkg)
        host_root.close()
        host_pkg.close()
        close_owned_pypedal_log_handlers()
