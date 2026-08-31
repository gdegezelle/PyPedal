"""User-facing inbreeding summaries round for display; stored F is unchanged."""
from PyPedal.pyp_io import format_display_coefficient, summary_inbreeding
from PyPedal.pyp_nrm import compute_inbreeding_stats


GRIFFON_RESIDUE = -4.9960036108132044e-15
GRIFFON_MAX = 0.546875


def test_tiny_negative_displays_as_unsigned_zero():
    assert format_display_coefficient(GRIFFON_RESIDUE) == "0.000000"
    assert format_display_coefficient(GRIFFON_RESIDUE).startswith("0.")
    assert "-0.000000" not in format_display_coefficient(GRIFFON_RESIDUE)
    assert format_display_coefficient(-0.0) == "0.000000"


def test_exact_max_displays_without_residue():
    assert format_display_coefficient(0.546875000000005) == "0.546875"
    assert format_display_coefficient(GRIFFON_MAX) == "0.546875"


def test_summary_rounds_display_and_keeps_raw_metadata():
    metadata = {
        "all": {
            "f_count": 98001,
            "f_min": GRIFFON_RESIDUE,
            "f_max": GRIFFON_MAX,
            "f_avg": 0.09313044278029989,
        },
        "nonzero": {
            "f_count": 84442,
            "f_min": 1e-16,
            "f_max": GRIFFON_MAX,
        },
    }
    text = summary_inbreeding(metadata)
    assert metadata["all"]["f_min"] == GRIFFON_RESIDUE
    assert "\tf_count\t98001" in text
    assert "\tf_count\t84442" in text
    assert "\tf_min\t0.000000" in text
    assert "-0.000000" not in text
    assert "\tf_max\t0.546875" in text
    assert "Animals with computed F > 0:" in text
    assert "Animals with non-zero CoI" not in text


def test_compute_inbreeding_stats_keeps_unrounded_residue():
    fx = {1: GRIFFON_RESIDUE, 2: 0.0, 3: GRIFFON_MAX}
    stats = compute_inbreeding_stats(fx)
    assert stats["all"]["f_min"] == GRIFFON_RESIDUE
    assert stats["all"]["f_max"] == GRIFFON_MAX
    assert stats["nonzero"]["f_count"] == 1
    text = summary_inbreeding(stats)
    assert "0.000000" in text
    assert "-0.000000" not in text
