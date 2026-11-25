#!/usr/bin/env sh
# Launch the GUI on macOS / Linux. Keep this file next to gui.py and run
# ./run_gui.sh  (make it executable once with: chmod +x run_gui.sh).
#
# Looks through the Python interpreters on this machine and prefers one that
# already has numpy, scipy and matplotlib, rather than the first on PATH.
cd "$(dirname "$0")" || exit 1

candidates=""
[ -n "$PYTHON" ] && candidates="$PYTHON"
[ -n "$VIRTUAL_ENV" ] && [ -x "$VIRTUAL_ENV/bin/python" ] && candidates="$candidates $VIRTUAL_ENV/bin/python"
[ -x ".venv/bin/python" ] && candidates="$candidates .venv/bin/python"
for name in python3 python; do
    path=$(command -v "$name" 2>/dev/null) && candidates="$candidates $path"
done
for path in "$HOME/anaconda3/bin/python" "$HOME/miniconda3/bin/python" \
            "$HOME/miniforge3/bin/python" /opt/homebrew/bin/python3 /usr/local/bin/python3; do
    [ -x "$path" ] && candidates="$candidates $path"
done

PY=""
FALLBACK=""
for path in $candidates; do
    "$path" -c "import sys" >/dev/null 2>&1 || continue
    [ -z "$FALLBACK" ] && FALLBACK="$path"
    if "$path" -c "import numpy, scipy, matplotlib" >/dev/null 2>&1; then
        PY="$path"
        break
    fi
done

if [ -z "$PY" ]; then
    if [ -z "$FALLBACK" ]; then
        echo "No Python interpreter found. Install Python 3.9 or newer." >&2
        exit 1
    fi
    PY="$FALLBACK"
    echo "Using Python: $PY"
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
else
    echo "Using Python: $PY"
fi

exec "$PY" gui.py
