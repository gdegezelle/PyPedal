"""Isolated QSettings access for the Qt desktop.

Presentation widgets do not read the platform registry or plist directly.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from PySide6.QtCore import QByteArray, QSettings

ORGANIZATION_NAME = "PyPedal"
APPLICATION_NAME = "PyPedal"
MAX_RECENT_FILES = 10

_KEY_GEOMETRY = "window/geometry"
_KEY_STATE = "window/state"
_KEY_DIRECTORY = "open/lastDirectory"
_KEY_PEDFORMAT = "open/lastPedformat"
_KEY_SEPARATOR = "open/lastSeparator"
_KEY_RENUMBER = "open/lastRenumber"
_KEY_RECENT = "open/recentFiles"

DEFAULT_PEDFORMAT = "asd"
DEFAULT_SEPARATOR = " "
DEFAULT_RENUMBER = True


class DesktopSettings:
    """Typed facade over a ``QSettings`` store."""

    def __init__(self, settings: QSettings | None = None) -> None:
        self._qs = settings or QSettings(ORGANIZATION_NAME, APPLICATION_NAME)

    @property
    def qsettings(self) -> QSettings:
        return self._qs

    def window_geometry(self) -> QByteArray | None:
        value = self._qs.value(_KEY_GEOMETRY)
        return value if isinstance(value, QByteArray) else None

    def set_window_geometry(self, geometry: QByteArray) -> None:
        self._qs.setValue(_KEY_GEOMETRY, geometry)

    def window_state(self) -> QByteArray | None:
        value = self._qs.value(_KEY_STATE)
        return value if isinstance(value, QByteArray) else None

    def set_window_state(self, state: QByteArray) -> None:
        self._qs.setValue(_KEY_STATE, state)

    def last_directory(self) -> Path | None:
        value = self._qs.value(_KEY_DIRECTORY, "")
        text = str(value or "")
        return Path(text) if text else None

    def set_last_directory(self, directory: Path) -> None:
        self._qs.setValue(_KEY_DIRECTORY, str(directory))

    def last_pedformat(self) -> str:
        value = self._qs.value(_KEY_PEDFORMAT, DEFAULT_PEDFORMAT)
        text = str(value or "").strip()
        return text or DEFAULT_PEDFORMAT

    def set_last_pedformat(self, pedformat: str) -> None:
        self._qs.setValue(_KEY_PEDFORMAT, pedformat)

    def last_separator(self) -> str:
        value = self._qs.value(_KEY_SEPARATOR, DEFAULT_SEPARATOR)
        if value is None:
            return DEFAULT_SEPARATOR
        return str(value)

    def set_last_separator(self, separator: str) -> None:
        self._qs.setValue(_KEY_SEPARATOR, separator)

    def last_renumber(self) -> bool:
        value = self._qs.value(_KEY_RENUMBER, DEFAULT_RENUMBER)
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            return value.lower() in {"1", "true", "yes"}
        return bool(value) if value is not None else DEFAULT_RENUMBER

    def set_last_renumber(self, renumber: bool) -> None:
        self._qs.setValue(_KEY_RENUMBER, bool(renumber))

    def recent_files(self) -> list[Path]:
        raw = self._qs.value(_KEY_RECENT, [])
        if raw is None:
            return []
        if isinstance(raw, str):
            items = [raw] if raw else []
        elif isinstance(raw, Sequence) and not isinstance(raw, (bytes, bytearray)):
            items = [str(item) for item in raw]
        else:
            items = [str(raw)]
        paths: list[Path] = []
        seen: set[str] = set()
        for item in items:
            if not item or item in seen:
                continue
            seen.add(item)
            paths.append(Path(item))
            if len(paths) >= MAX_RECENT_FILES:
                break
        return paths

    def remember_successful_open(
        self,
        source: Path,
        pedformat: str,
        separator: str,
        renumber: bool,
    ) -> None:
        """Record a successful load: options plus most-recent-first path list."""
        self.set_last_directory(source.parent)
        self.set_last_pedformat(pedformat)
        self.set_last_separator(separator)
        self.set_last_renumber(renumber)
        resolved = str(source)
        current = [str(path) for path in self.recent_files() if str(path) != resolved]
        self._qs.setValue(_KEY_RECENT, [resolved, *current][:MAX_RECENT_FILES])
        self._qs.sync()

    def remove_recent_file(self, source: Path) -> None:
        target = str(source)
        remaining = [str(path) for path in self.recent_files() if str(path) != target]
        self._qs.setValue(_KEY_RECENT, remaining)
        self._qs.sync()
