#!/usr/bin/env sh
# Launch the GUI on macOS / Linux. Keep this file next to gui.py and run
# ./run_gui.sh  (make it executable once with: chmod +x run_gui.sh).
cd "$(dirname "$0")" || exit 1

PY="${PYTHON:-}"
[ -z "$PY" ] && [ -n "$VIRTUAL_ENV" ] && [ -x "$VIRTUAL_ENV/bin/python" ] && PY="$VIRTUAL_ENV/bin/python"
[ -z "$PY" ] && [ -x ".venv/bin/python" ] && PY=".venv/bin/python"
[ -z "$PY" ] && command -v python3 >/dev/null 2>&1 && PY="python3"
[ -z "$PY" ] && command -v python  >/dev/null 2>&1 && PY="python"

if [ -z "$PY" ]; then
    echo "No Python interpreter found. Install Python 3.9 or newer." >&2
    exit 1
fi
echo "Using Python: $PY"

if ! "$PY" -c "import numpy, scipy, matplotlib" >/dev/null 2>&1; then
    echo "Some required packages are missing (numpy / scipy / matplotlib)."
    printf "Install them now with pip? [Y/n] "
    read -r reply
    case "$reply" in
        [Nn]*) echo "Cannot start without them. Exiting."; exit 1 ;;
    esac
    "$PY" -m pip install -r requirements.txt || {
        echo "Installation failed. Try: $PY -m pip install --user -r requirements.txt" >&2
        exit 1
    }
fi

exec "$PY" gui.py
