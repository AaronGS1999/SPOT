"""Build the SPOT desktop app from this repository.

    python packaging/build.py macos  -> dist/SPOT.app   (Apple Silicon and Intel; zip it or put it in a dmg)
    python packaging/build.py linux  -> dist/SPOT-linux/ (x86_64 and aarch64, with install.sh)

The app holds the SPOT wheel, the launcher, the environment file and micromamba for each supported CPU;
the engine environment itself is created on the user's computer at first launch. Needs network access to
conda.anaconda.org (micromamba) and Pillow (icon).
"""

from __future__ import annotations

import io
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "packaging"
DIST = ROOT / "dist"
MICROMAMBA = "1.5.10-0"
PLATFORMS = {"macos": ["osx-arm64", "osx-64"], "linux": ["linux-64", "linux-aarch64"]}


def version() -> str:
    ns: dict[str, str] = {}
    exec((ROOT / "src" / "spot" / "__init__.py").read_text(), ns)
    return ns["__version__"]


def micromamba(platform: str, dest: Path) -> None:
    url = f"https://conda.anaconda.org/conda-forge/{platform}/micromamba-{MICROMAMBA}.tar.bz2"
    data = urllib.request.urlopen(url, timeout=120).read()
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:bz2") as tar:
        member = tar.getmember("bin/micromamba")
        fh = tar.extractfile(member)
        assert fh is not None
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(fh.read())
    dest.chmod(0o755)


def wheel(dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [sys.executable, "-m", "pip", "wheel", "--no-deps", "-q", "-w", str(dest), str(ROOT)], check=True
    )
    return next(dest.glob("spot-*.whl"))


def resources(res: Path, target: str) -> None:
    res.mkdir(parents=True, exist_ok=True)
    (res / "version").write_text(version() + "\n")
    shutil.copy(PKG / "environment-app.yml", res / "environment-app.yml")
    wheel(res)
    for platform in PLATFORMS[target]:
        micromamba(platform, res / "micromamba" / platform)


def icon_png(dest: Path) -> None:
    shutil.copy(ROOT / "src" / "spot" / "gui" / "assets" / "spot.png", dest)


def macos() -> Path:
    from PIL import Image

    app = DIST / "SPOT.app"
    shutil.rmtree(app, ignore_errors=True)
    contents = app / "Contents"
    (contents / "MacOS").mkdir(parents=True)
    plist = (PKG / "macos" / "Info.plist").read_text().replace("@VERSION@", version())
    (contents / "Info.plist").write_text(plist)
    (contents / "PkgInfo").write_text("APPL????")
    launcher = contents / "MacOS" / "SPOT"
    shutil.copy(PKG / "launcher.sh", launcher)
    launcher.chmod(0o755)
    res = contents / "Resources"
    resources(res, "macos")
    Image.open(ROOT / "src" / "spot" / "gui" / "assets" / "spot.png").save(
        res / "SPOT.icns",
        sizes=[(16, 16), (32, 32), (64, 64), (128, 128), (256, 256), (512, 512), (1024, 1024)],
    )
    return app


def linux() -> Path:
    out = DIST / "SPOT-linux"
    shutil.rmtree(out, ignore_errors=True)
    (out / "bin").mkdir(parents=True)
    launcher = out / "bin" / "spot-app"
    shutil.copy(PKG / "launcher.sh", launcher)
    launcher.chmod(0o755)
    resources(out / "Resources", "linux")
    icon_png(out / "Resources" / "spot.png")
    shutil.copy(PKG / "linux" / "install.sh", out / "install.sh")
    (out / "install.sh").chmod(0o755)
    return out


if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "macos"
    print({"macos": macos, "linux": linux}[target]())
