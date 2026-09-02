"""Read-only pedigree metadata page."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFormLayout,
    QLabel,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from PyPedal.application import PedigreeSession

_EMPTY = "No pedigree is open."
_HINT = "Use File -> Open… to open a pedigree."
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
        self._hint = QLabel(_HINT)
        self._hint.setObjectName("metadata_open_hint")
        self._hint.setWordWrap(True)
        self._labels: dict[str, QLabel] = {}

        self._form_host = QWidget()
        self._form_host.setObjectName("metadata_form")
        self._form_host.setSizePolicy(
            QSizePolicy.Policy.Preferred,
            QSizePolicy.Policy.Maximum,
        )
        self._form = QFormLayout(self._form_host)
        self._form.setContentsMargins(0, 0, 0, 0)
        self._form.setHorizontalSpacing(12)
        self._form.setVerticalSpacing(6)
        self._form.setFormAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self._form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self._form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
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
            label.setWordWrap(key == "filename")
            self._labels[key] = label
            self._form.addRow(title, label)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)
        layout.addWidget(self._banner)
        layout.addWidget(self._hint)
        layout.addWidget(self._form_host)
        layout.addStretch(1)
        self.show_empty()

    def show_empty(self) -> None:
        self._banner.setText(_EMPTY)
        self._banner.show()
        self._hint.show()
        self._form_host.hide()
        for label in self._labels.values():
            label.setText(_ABSENT)

    def show_session(self, session: PedigreeSession) -> None:
        pedigree = session.pedigree
        if pedigree is None:
            self.show_empty()
            return
        self._banner.hide()
        self._hint.hide()
        self._form_host.show()
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
