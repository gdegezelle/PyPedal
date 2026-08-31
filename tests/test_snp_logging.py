"""Missing SNP data is a normal condition, not an error."""
import os
import subprocess
import sys
import textwrap

from _pedhelpers import REPO


def _run(script):
    env = os.environ.copy()
    env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
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


def test_missing_snp_data_is_logged_at_debug(tmp_path):
    pedfile = tmp_path / "t.ped"
    pedfile.write_text("1 0 0\n2 0 0\n3 1 2\n", encoding="utf-8")
    message = "SNP data is not available. No IDs need to be renumbered."
    script = textwrap.dedent(
        f"""
        from PyPedal.pyp_newclasses import load_pedigree
        ped = load_pedigree(options={{
            "pedfile": {str(pedfile)!r},
            "pedformat": "asd",
            "messages": "quiet",
            "pedigree_summary": 0,
            "renumber": True,
        }})
        text = open(ped.kw["logfile"], encoding="utf-8").read()
        print(text)
        """
    )
    result = _run(script)
    log_text = result.stdout
    assert message in log_text
    debug_line = [
        line for line in log_text.splitlines() if message in line
    ]
    assert debug_line, log_text
    assert "DEBUG" in debug_line[0]
    assert "ERROR" not in debug_line[0]
    assert message not in result.stderr
