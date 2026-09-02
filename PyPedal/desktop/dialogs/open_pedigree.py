"""Pedigree open dialog using application-layer options."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QWidget,
)

from PyPedal.application import PedigreeOpenOptions, normalize_sepchar
from PyPedal.desktop.settings import (
    DEFAULT_PEDFORMAT,
    DEFAULT_SEPARATOR,
    DesktopSettings,
)


class OpenPedigreeDialog(QDialog):
    """Compact path/format/separator/renumber form."""

    def __init__(
        self,
        settings: DesktopSettings,
        parent: QWidget | None = None,
        initial_path: Path | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Open pedigree")
        self._settings = settings

        self.path_edit = QLineEdit()
        self.path_edit.setObjectName("open_path")
        browse = QPushButton("Browse…")
        browse.setObjectName("open_browse")
        browse.clicked.connect(self._browse)

        path_row = QHBoxLayout()
        path_row.addWidget(self.path_edit)
        path_row.addWidget(browse)

        self.format_edit = QLineEdit(settings.last_pedformat() or DEFAULT_PEDFORMAT)
        self.format_edit.setObjectName("open_pedformat")
        self.separator_edit = QLineEdit()
        self.separator_edit.setObjectName("open_separator")
        self.separator_edit.setPlaceholderText("space")
        last_sep = settings.last_separator()
        if last_sep != DEFAULT_SEPARATOR:
            self.separator_edit.setText(last_sep)
        self.renumber_box = QCheckBox("Renumber")
        self.renumber_box.setObjectName("open_renumber")
        self.renumber_box.setChecked(settings.last_renumber())

        if initial_path is not None:
            self.path_edit.setText(str(initial_path))

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        form = QFormLayout(self)
        form.addRow("File", path_row)
        form.addRow("Format", self.format_edit)
        form.addRow("Separator", self.separator_edit)
        form.addRow(self.renumber_box)
        form.addRow(buttons)

    def _browse(self) -> None:
        directory = self._settings.last_directory()
        start = str(directory) if directory is not None else ""
        path, _filter = QFileDialog.getOpenFileName(
            self,
            "Open pedigree",
            start,
            "Pedigree files (*.ped);;All files (*)",
        )
        if path:
            self.path_edit.setText(path)

    def selected_path(self) -> Path | None:
        text = self.path_edit.text().strip()
        return Path(text) if text else None

    def selected_options(self) -> PedigreeOpenOptions:
        return PedigreeOpenOptions(
            pedformat=self.format_edit.text().strip() or DEFAULT_PEDFORMAT,
            separator=normalize_sepchar(self.separator_edit.text()),
            renumber=self.renumber_box.isChecked(),
            messages="quiet",
            pedigree_summary=0,
        ).normalized()
