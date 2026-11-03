#!/bin/bash
# Install SPOT for this user: copies the app to ~/.local/opt/spot and adds it to the applications menu.
set -eu
SRC="$(cd "$(dirname "$0")" && pwd)"
DEST="${XDG_DATA_HOME:-$HOME/.local/share}/../opt/spot"
mkdir -p "$DEST" "$HOME/.local/share/applications" "$HOME/.local/share/icons"
cp -R "$SRC/bin" "$SRC/Resources" "$DEST/"
DEST="$(cd "$DEST" && pwd)"
cp "$SRC/Resources/spot.png" "$HOME/.local/share/icons/spot.png"
cat > "$HOME/.local/share/applications/spot.desktop" <<DESKTOP
[Desktop Entry]
Type=Application
Name=SPOT
Comment=Synteny-based Positional Orthology Test
Exec="$DEST/bin/spot-app" %f
Icon=$HOME/.local/share/icons/spot.png
Terminal=false
Categories=Science;Biology;
MimeType=application/x-yaml;
DESKTOP
echo "SPOT installed. Open it from the applications menu."
