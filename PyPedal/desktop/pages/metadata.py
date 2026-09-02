"""Read-only pedigree metadata page."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QFormLayout, QLabel, QVBoxLayout, QWidget

from PyPedal.application import PedigreeSession

_EMPTY = "No pedigree is open."
_ABSENT = "—"


def _text(value: object) -> str:
    if value is None or value == "":
        return _ABSENT
    return str(value)


class MetadataPage(QWidget):
    """Labeled metadata from the loaded pedigree snapshot."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._banner = QLabel(_EMPTY)
        self._banner.setObjectName("metadata_empty")
        self._labels: dict[str, QLabel] = {}
        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        for key, title in (
            ("name", "Name"),
            ("filename", "File"),
            ("pedcode", "Format"),
            ("num_records", "Records"),
            ("num_unique_sires", "Unique sires"),
            ("num_unique_dams", "Unique dams"),
            ("num_unique_gens", "Unique generations"),
            ("num_unique_years", "Unique years"),
            ("num_unique_founders", "Unique founders"),
            ("num_unique_herds", "Unique herds"),
            ("num_implicit_parents", "Implicit parents"),
            ("renumbered", "Renumbered"),
            ("snp_count", "SNP count"),
            ("num_unique_fields", "Unique user fields"),
        ):
            label = QLabel(_ABSENT)
            label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            label.setObjectName(f"metadata_{key}")
            self._labels[key] = label
            form.addRow(title, label)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.addWidget(self._banner)
        layout.addLayout(form)
        self.show_empty()

    def show_empty(self) -> None:
        self._banner.setText(_EMPTY)
        self._banner.show()
        for label in self._labels.values():
            label.setText(_ABSENT)

    def show_session(self, session: PedigreeSession) -> None:
        pedigree = session.pedigree
        if pedigree is None:
            self.show_empty()
            return
        self._banner.hide()
        metadata = pedigree.metadata
        values = {
            "name": getattr(metadata, "name", None),
            "filename": session.source_path,
            "pedcode": getattr(metadata, "pedcode", None),
            "num_records": getattr(metadata, "num_records", None),
            "num_unique_sires": getattr(metadata, "num_unique_sires", None),
            "num_unique_dams": getattr(metadata, "num_unique_dams", None),
            "num_unique_gens": getattr(metadata, "num_unique_gens", None),
            "num_unique_years": getattr(metadata, "num_unique_years", None),
            "num_unique_founders": getattr(metadata, "num_unique_founders", None),
            "num_unique_herds": getattr(metadata, "num_unique_herds", None),
            "num_implicit_parents": getattr(metadata, "num_implicit_parents", None),
            "renumbered": getattr(metadata, "renumbered", None),
            "snp_count": getattr(metadata, "snp_count", None),
            "num_unique_fields": getattr(metadata, "num_unique_fields", None),
        }
        pedcode = str(values["pedcode"] or "")
        if "u" not in pedcode:
            values["num_unique_fields"] = None
        for key, label in self._labels.items():
            label.setText(_text(values.get(key)))
