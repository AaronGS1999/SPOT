#!/bin/bash
# SPOT launcher (macOS .app and Linux). On first launch, or after an update, it creates the engine
# environment with the bundled micromamba (Python, Qt, BLAST, minimap2, SyRI, Gepard and Java) and
# installs SPOT into it; afterwards it only starts the app.
set -u

RES="$(cd "$(dirname "${BASH_SOURCE[0]}")/../Resources" && pwd)"
VERSION="$(cat "$RES/version")"
if [ "$(uname -s)" = "Darwin" ]; then
    SUPPORT="$HOME/Library/Application Support/SPOT"
else
    SUPPORT="${XDG_DATA_HOME:-$HOME/.local/share}/spot"
fi
ENV="$SUPPORT/env"
LOG="$SUPPORT/install.log"
MARK="$ENV/.spot-$VERSION"
mkdir -p "$SUPPORT"

notify() {
    if [ "$(uname -s)" = "Darwin" ]; then
        osascript -e "display notification \"$1\" with title \"SPOT\"" >/dev/null 2>&1 || true
    elif command -v notify-send >/dev/null 2>&1; then
        notify-send SPOT "$1" || true
    fi
}

fail() {
    if [ "$(uname -s)" = "Darwin" ]; then
        osascript -e "display alert \"SPOT could not be set up\" message \"$1\n\nDetails: $LOG\" as critical" \
            >/dev/null 2>&1
    else
        echo "SPOT could not be set up: $1 (details: $LOG)" >&2
    fi
    exit 1
}

if [ ! -f "$MARK" ]; then
    case "$(uname -s)-$(uname -m)" in
        Darwin-arm64) MM="$RES/micromamba/osx-arm64" ;;
        Darwin-x86_64) MM="$RES/micromamba/osx-64" ;;
        Linux-x86_64) MM="$RES/micromamba/linux-64" ;;
        Linux-aarch64) MM="$RES/micromamba/linux-aarch64" ;;
        *) fail "unsupported system $(uname -s) $(uname -m)" ;;
    esac
    [ -x "$MM" ] || fail "micromamba for $(uname -m) is missing from the app"
    notify "Setting up SPOT $VERSION. The first launch downloads its tools and takes a few minutes."
    {
        echo "== $(date) setting up SPOT $VERSION"
        export MAMBA_ROOT_PREFIX="$SUPPORT/mamba"
        if [ -d "$ENV/conda-meta" ]; then
            "$MM" install -y -p "$ENV" -f "$RES/environment-app.yml"
        else
            "$MM" create -y -p "$ENV" -f "$RES/environment-app.yml"
        fi
    } >>"$LOG" 2>&1 || fail "the tools could not be installed (is the computer online?)"
    "$ENV/bin/python" -m pip install --no-deps --force-reinstall "$RES"/spot-*.whl >>"$LOG" 2>&1 \
        || fail "SPOT could not be installed"
    rm -f "$ENV"/.spot-*
    touch "$MARK"
    notify "SPOT is ready."
fi

export PATH="$ENV/bin:$PATH"
[ -d "$ENV/lib/jvm" ] && export JAVA_HOME="$ENV/lib/jvm"
exec "$ENV/bin/python" -m spot.gui "$@"
