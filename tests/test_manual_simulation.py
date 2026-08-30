"""Simulation and example-script product checks."""
from test_examples_integration import KNOWN_FAILING
from test_manual_pages import USER

OPTIONS = USER / "configuration.md"
EXAMPLES = USER / "recipes.md"
LIMITATIONS = USER / "limitations.md"


def test_ph6_options_do_not_claim_simulate_early_return():
    text = OPTIONS.read_text(encoding="utf-8")
    assert "returns** before" not in text
    assert "returns before" not in text


def test_ph6_examples_do_not_list_inbreeding3_4_as_known_failing():
    text = EXAMPLES.read_text(encoding="utf-8")
    assert "IndexError" not in text
    assert "metadata` stays `{}`" not in text


def test_ph6_limitations_name_dyad_census_bound():
    text = LIMITATIONS.read_text(encoding="utf-8")
    assert "dyad_census" in text
    assert "NetworkX 3" in text


def test_ph6_known_failing_is_empty_after_dyad_cleanup():
    assert KNOWN_FAILING == {}
