"""Allow `python -m PyPedal` to launch the desktop GUI.

SPDX-License-Identifier: LGPL-2.1-or-later
"""

from PyPedal.desktop.main import main

if __name__ == "__main__":
    raise SystemExit(main())
