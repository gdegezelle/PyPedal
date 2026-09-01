"""Direct tests of the private parse helpers. Not a public API."""
from PyPedal._pyp_parse import (
    PEDFORMAT_ABSENT,
    PedigreeRecordSource,
    build_pedformat_locations,
    canonicalize_pedformat,
    implicit_parent_locations,
    iter_implicit_parent_tokens,
)

CODES = [
    "a", "s", "d", "g", "x", "b", "f", "r", "n",
    "y", "l", "e", "p", "A", "S", "D", "L", "Z",
    "h", "H", "u", "T", "P", "G", "Y",
]


def test_canonicalize_keeps_case_and_maps_Z_and_unknown():
    canonical, events = canonicalize_pedformat("asdxZqA", CODES)
    assert canonical == ["a", "s", "d", "x", ".", ".", "A"]
    assert events == [("Z", "Z"), ("invalid", "q")]


def test_build_locations_asd():
    locations, critical, collision, _debug = build_pedformat_locations(
        list("asd"), alleles_sepchar="/", sepchar=" "
    )
    assert critical == 0
    assert collision is False
    assert locations["animal"] == 0
    assert locations["sire"] == 1
    assert locations["dam"] == 2
    assert locations["sex"] == PEDFORMAT_ABSENT
    assert locations["herd"] == PEDFORMAT_ABSENT


def test_build_locations_string_identity():
    locations, critical, _collision, _debug = build_pedformat_locations(
        list("ASD"), alleles_sepchar="/", sepchar=","
    )
    assert critical == 0
    assert locations["animal"] == 0
    assert locations["sire"] == 1
    assert locations["dam"] == 2


def test_missing_required_codes_increment_critical_once():
    _locations, critical, _collision, _debug = build_pedformat_locations(
        list("sd"), alleles_sepchar="/", sepchar=" "
    )
    assert critical == 1
    _locations, critical, _collision, _debug = build_pedformat_locations(
        list("xy"), alleles_sepchar="/", sepchar=" "
    )
    assert critical == 3


def test_h_without_H_leaves_herd_absent():
    locations, critical, _collision, _debug = build_pedformat_locations(
        list("asdh"), alleles_sepchar="/", sepchar=" "
    )
    assert critical == 0
    assert locations["herd"] == PEDFORMAT_ABSENT


def test_H_populates_herd():
    locations, _critical, _collision, _debug = build_pedformat_locations(
        list("asdH"), alleles_sepchar="/", sepchar=" "
    )
    assert locations["herd"] == 3


def test_alleles_sepchar_collision_disables_alleles():
    locations, _critical, collision, _debug = build_pedformat_locations(
        list("asdL"), alleles_sepchar=" ", sepchar=" "
    )
    assert collision is True
    assert locations["alleles"] == PEDFORMAT_ABSENT


def test_record_source_file_and_text_and_db(tmp_path):
    path = tmp_path / "tiny.ped"
    path.write_text("1 0 0\n2 0 0\n", encoding="utf-8")

    class _Log:
        def warning(self, *args, **kwargs):
            pass

        def info(self, *args, **kwargs):
            pass

    file_src = PedigreeRecordSource(str(path))
    file_lines = []
    while True:
        line = file_src.readline(len(file_lines), _Log())
        if line is False or not line or line.strip() == "":
            break
        file_lines.append(line.strip())
    file_src.close()

    text_src = PedigreeRecordSource(str(path), textstream="1 0 0\n2 0 0\n")
    text_lines = []
    while True:
        line = text_src.readline(len(text_lines), _Log())
        if line is False:
            break
        if not line or line.strip() == "":
            break
        text_lines.append(line.strip())
    text_src.close()

    db_src = PedigreeRecordSource(str(path), dbstream=[("1", "0", "0"), ("2", "0", "0")])
    db_lines = []
    while True:
        line = db_src.readline(len(db_lines), _Log())
        if line is False:
            break
        db_lines.append(line)
    db_src.close()

def test_record_source_known_total_is_cheap(tmp_path):
    path = tmp_path / "tiny.ped"
    path.write_text("1 0 0\n2 0 0\n", encoding="utf-8")
    file_src = PedigreeRecordSource(str(path))
    assert file_src.known_total is None
    file_src.close()
    text_src = PedigreeRecordSource(str(path), textstream="1 0 0\n2 0 0\n")
    assert text_src.known_total == 2
    text_src.close()
    db_src = PedigreeRecordSource(
        str(path), dbstream=[("1", "0", "0"), ("2", "0", "0")]
    )
    assert db_src.known_total == 2
    db_src.close()



def test_text_source_without_trailing_newline_drops_last_record(tmp_path):
    class _Log:
        def warning(self, *args, **kwargs):
            pass

        def info(self, *args, **kwargs):
            pass

    src = PedigreeRecordSource(str(tmp_path / "x.ped"), textstream="1 0 0\n2 0 0")
    lines = []
    while True:
        line = src.readline(len(lines), _Log())
        if line is False:
            break
        if not line or line.strip() == "":
            break
        lines.append(line.strip())
    src.close()
    assert lines == ["1 0 0"]


def test_implicit_parent_tokens_sire_dam_both_and_shared():
    tokens = list(
        iter_implicit_parent_tokens(
            sires={1: 1},
            dams={},
            idmap={2: 2},
            pedformat="asd",
            missing_parent=0,
            missing_name="Unknown_Name",
        )
    )
    assert tokens == [("sire", 1)]

    tokens = list(
        iter_implicit_parent_tokens(
            sires={},
            dams={3: 3},
            idmap={2: 2},
            pedformat="asd",
            missing_parent=0,
            missing_name="Unknown_Name",
        )
    )
    assert tokens == [("dam", 3)]

    tokens = list(
        iter_implicit_parent_tokens(
            sires={1: 1},
            dams={2: 2},
            idmap={3: 3},
            pedformat="asd",
            missing_parent=0,
            missing_name="Unknown_Name",
        )
    )
    assert tokens == [("sire", 1), ("dam", 2)]

    tokens = list(
        iter_implicit_parent_tokens(
            sires={1: 1},
            dams={1: 1},
            idmap={3: 3, 4: 4},
            pedformat="asd",
            missing_parent=0,
            missing_name="Unknown_Name",
        )
    )
    assert tokens == [("sire", 1)]


def test_implicit_parent_tokens_string_ids():
    tokens = list(
        iter_implicit_parent_tokens(
            sires={"sireX": "sireX"},
            dams={},
            idmap={"child": "child"},
            pedformat="ASD",
            missing_parent=0,
            missing_name="Unknown_Name",
        )
    )
    assert tokens == [("sire", "sireX")]


def test_implicit_parent_locations_zero_sire_dam():
    locations = {"animal": 0, "sire": 1, "dam": 2, "sex": 3}
    null = implicit_parent_locations(locations)
    assert null["animal"] == 0
    assert null["sire"] == 1
    assert null["dam"] == 2
    assert null["sex"] == PEDFORMAT_ABSENT
    assert locations["sex"] == 3
