"""GUI preview caps, busy-state intent, and neutral zero-inbreeding wording."""
from types import SimpleNamespace

from PyPedal.pyp_app import (
    GUI_PREVIEW_ROWS,
    GuiProgressBridge,
    _format_inbreeding,
    _list_animals,
    format_preview_caption,
    gui_control_states,
    gui_progress_mode,
)


def _animal(i):
    return SimpleNamespace(
        animalID=i,
        sireID=0,
        damID=0,
        by=2000,
        sex="u",
        name=f"a{i}",
    )


def test_list_animals_caps_preview_and_reports_total():
    pedigree = SimpleNamespace(
        pedigree=[_animal(i) for i in range(1, GUI_PREVIEW_ROWS + 2)]
    )
    text = _list_animals(pedigree)
    assert format_preview_caption(GUI_PREVIEW_ROWS, GUI_PREVIEW_ROWS + 1) in text
    assert text.count("\n") >= GUI_PREVIEW_ROWS
    assert f"{GUI_PREVIEW_ROWS + 1:>8}" not in text.split("Name", 1)[-1]


def test_list_animals_does_not_caption_small_pedigrees():
    pedigree = SimpleNamespace(pedigree=[_animal(1), _animal(2)])
    text = _list_animals(pedigree)
    assert "Showing" not in text
    assert "a1" in text and "a2" in text


def test_inbreeding_preview_caps_and_mentions_result_file():
    fx = {i: 0.125 for i in range(1, 601)}
    result = {
        "metadata": {
            "all": {"f_count": 600, "f_min": 0.125, "f_max": 0.125, "f_avg": 0.125},
            "nonzero": {"f_count": 600, "f_min": 0.125, "f_max": 0.125, "f_avg": 0.125},
        },
        "fx": fx,
    }
    text = _format_inbreeding(result, result_file="/tmp/dogs_inbreeding.dat")
    assert len(fx) == 600
    assert format_preview_caption(500, 600) in text
    assert "Full coefficients are in dogs_inbreeding.dat" in text
    assert "  500: 0.125000" in text
    assert "  501: 0.125000" not in text


def test_preview_caption_for_griffon_scale():
    assert format_preview_caption(500, 98001) == "Showing 500 of 98,001"


def test_zero_inbreeding_wording_is_neutral():
    result = {
        "metadata": {
            "all": {"f_count": 3, "f_min": 0.0, "f_max": 0.0, "f_avg": 0.0},
            "nonzero": {"f_count": 0, "f_min": 0.0, "f_max": 0.0, "f_avg": 0.0},
        },
        "fx": {1: 0.0, 2: 0.0, 3: 0.0},
    }
    text = _format_inbreeding(result)
    assert "No inbreeding in this pedigree (every coefficient is 0)." in text
    assert "new_lacy.ped" not in text
    assert "mrode.ped" not in text
    assert "hartlandclark.ped" not in text


def test_busy_disables_open_and_analyses_but_not_about():
    busy = gui_control_states(True)
    idle = gui_control_states(False)
    assert busy == {"open": False, "analyses": False, "about": True}
    assert idle == {"open": True, "analyses": True, "about": True}


def test_gui_progress_bridge_does_not_configure_widgets():
    calls = []

    class Widget:
        def configure(self, **kwargs):
            calls.append(kwargs)

    widget = Widget()
    bridge = GuiProgressBridge()
    bridge(1, 10)
    bridge(10, 10)
    assert bridge.latest == (10, 10)
    assert calls == []
    widget.configure(mode="determinate")
    assert calls == [{"mode": "determinate"}]


def test_gui_progress_mode_switches_on_known_total():
    assert gui_progress_mode(5, None) == ("indeterminate", None)
    assert gui_progress_mode(5, 10) == ("determinate", 0.5)

