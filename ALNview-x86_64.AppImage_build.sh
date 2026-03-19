#!/usr/bin/env bash
set -euo pipefail

# This comment must stay as it is build instructions for the reader
# podman run --rm -it -v "$PWD":/work -w /work almalinux:8 bash ./build-el8.sh

# Linux AppImage for ALNVIEW (CentOS 8 compatible build environment)

QT_VER="${QT_VER:-6.9.2}"
QT_SERIES="${QT_VER%.*}"              # e.g. 6.9
CACHE_DIR="${CACHE_DIR:-/work/_cache}"
NPROC="${NPROC:-$(nproc)}"

mkdir -p "$CACHE_DIR"

download_cached() {
  # download_cached <dest_path> <url>
  local dest="$1"
  local url="$2"
  if [ -s "$dest" ]; then
    echo "Using cached: $dest"
  else
    echo "Downloading: $url"
    wget -c -O "$dest" "$url"
  fi
}

# -------------------- Repos & deps --------------------
# (Tip: comment out dnf update to reduce repeated downloads)
# dnf -y update

dnf -y install dnf-plugins-core ca-certificates wget

# Enable optional repos (names differ by distro)
dnf config-manager --set-enabled powertools 2>/dev/null || true
dnf config-manager --set-enabled PowerTools 2>/dev/null || true
dnf config-manager --set-enabled crb 2>/dev/null || true
dnf -y makecache

# EPEL (needed on EL8 for xcb-util-cursor)
if ! rpm -q epel-release >/dev/null 2>&1; then
  dnf -y install https://dl.fedoraproject.org/pub/epel/epel-release-latest-8.noarch.rpm
fi
dnf -y makecache

# Toolchain + ninja
dnf -y install gcc-toolset-10 gcc-toolset-10-gcc gcc-toolset-10-gcc-c++ ninja-build

# gcc-toolset enable scripts can touch unset vars; set -u breaks that
set +u
source /opt/rh/gcc-toolset-10/enable
set -u

gcc --version | head -n1

# Build deps (+ fonts + fontconfig tools)
dnf -y install \
  git make perl python3 tar xz \
  cmake pkgconf-pkg-config \
  zlib-devel \
  libX11-devel libXext-devel libXrender-devel libXrandr-devel libXi-devel libXcursor-devel libXinerama-devel \
  libxcb-devel \
  xcb-util-devel xcb-util-wm-devel xcb-util-image-devel xcb-util-keysyms-devel xcb-util-renderutil-devel \
  xcb-util-cursor xcb-util-cursor-devel \
  libxkbcommon-devel libxkbcommon-x11-devel \
  mesa-libGL-devel mesa-libEGL-devel \
  freetype freetype-devel \
  fontconfig fontconfig-devel \
  dejavu-fonts-common dejavu-sans-fonts dejavu-sans-mono-fonts \
  file squashfs-tools patchelf

# -------------------- Fetch ALNVIEW --------------------
if [ -d ALNVIEW/.git ]; then
  git -C ALNVIEW fetch --depth 1 origin
  git -C ALNVIEW reset --hard origin/main 2>/dev/null || git -C ALNVIEW reset --hard origin/master
else
  git clone --depth 1 https://github.com/thegenemyers/ALNVIEW
fi

# viewer.pro hard-codes macx-clang; guard it so qmake works on Linux
perl -i -pe 's/^QMAKE_SPEC\s*=\s*macx-clang/macx: QMAKE_SPEC = macx-clang/' ALNVIEW/viewer.pro

# Upstream hard-codes the macOS font family "Monaco" in several places.
# On Linux that can resolve poorly or not at all, and in a self-contained AppImage
# it defeats font portability. Patch the source to use bundled DejaVu fonts instead:
#   - UI/default font      -> DejaVu Sans
#   - alignment/code views -> DejaVu Sans Mono
python3 - <<'PY_PATCH'
from pathlib import Path
p = Path('ALNVIEW/main_window.cpp')
s = p.read_text()
s = s.replace('setFont(QFont(tr("Monaco")));', 'setFont(QFont(QStringLiteral("DejaVu Sans")));')
s = s.replace('QFont(tr("Monaco"),11)', 'QFont(QStringLiteral("DejaVu Sans Mono"),11)')
s = s.replace('QFont(tr("Monaco"))', 'QFont(QStringLiteral("DejaVu Sans Mono"))')
p.write_text(s)
PY_PATCH

# The initial startup dialog is a custom Qt dialog (OpenDialog), not only the
# later native file picker. So we must make Qt itself load bundled fonts before
# that dialog is created. We do that in main.cpp with addApplicationFont().
python3 - <<'PY_MAINPATCH'
from pathlib import Path
p = Path('ALNVIEW/main.cpp')
s = p.read_text()
if '#include <QFontDatabase>' not in s:
    s = s.replace('#include <QApplication>', '#include <QApplication>\n#include <QFontDatabase>\n#include <QCoreApplication>\n#include <QFont>')
old = '  QApplication app(argc, argv);'
new = """  QCoreApplication::setAttribute(Qt::AA_DontUseNativeDialogs);
  QApplication app(argc, argv);

  QString fontdir = QCoreApplication::applicationDirPath() + QStringLiteral("/../share/fonts/truetype/dejavu/");
  int sansId = QFontDatabase::addApplicationFont(fontdir + QStringLiteral("DejaVuSans.ttf"));
  int monoId = QFontDatabase::addApplicationFont(fontdir + QStringLiteral("DejaVuSansMono.ttf"));

  QString sansFamily = QStringLiteral("DejaVu Sans");
  if (sansId >= 0) {
    QStringList fams = QFontDatabase::applicationFontFamilies(sansId);
    if (!fams.isEmpty())
      sansFamily = fams.first();
  }

  if (monoId < 0)
    fprintf(stderr, "Warning: could not load bundled DejaVu Sans Mono\\n");

  app.setFont(QFont(sansFamily));"""
if old not in s:
    raise SystemExit('Could not find QApplication construction in ALNVIEW/main.cpp')
s = s.replace(old, new, 1)
p.write_text(s)
PY_MAINPATCH

# -------------------- Build Qt (qtbase only; cached) --------------------
rm -rf qt-src qt-install qt-build
mkdir -p qt-src qt-install

QT_BASE_TARBALL="$CACHE_DIR/qtbase-everywhere-src-${QT_VER}.tar.xz"
QT_BASE_URL="https://download.qt.io/official_releases/qt/${QT_SERIES}/${QT_VER}/submodules/qtbase-everywhere-src-${QT_VER}.tar.xz"
download_cached "$QT_BASE_TARBALL" "$QT_BASE_URL"

tar -C qt-src -xf "$QT_BASE_TARBALL"
QT_SRC_DIR="$PWD/qt-src/qtbase-everywhere-src-${QT_VER}"
QT_PREFIX="$PWD/qt-install/${QT_VER}"

# Clean build dir every run to avoid generator mismatch (Makefiles vs Ninja)
rm -rf qt-build
mkdir -p qt-build
cd qt-build

"$QT_SRC_DIR/configure" \
  -prefix "$QT_PREFIX" \
  -opensource -confirm-license \
  -release -nomake examples -nomake tests \
  -xcb

cmake --build . --parallel "$NPROC"
cmake --install .

export PATH="$QT_PREFIX/bin:$PATH"
qmake -v

# Ensure xcb platform plugin exists (linuxdeploy-plugin-qt needs it)
test -f "$QT_PREFIX/plugins/platforms/libqxcb.so"

# -------------------- Build ALNVIEW --------------------
cd "$PWD/../ALNVIEW"
qmake viewer.pro
make -j"$NPROC"
cd "$PWD/.."

# -------------------- linuxdeploy (cached) --------------------
export APPIMAGE_EXTRACT_AND_RUN=1
export QMAKE="$QT_PREFIX/bin/qmake"

LINUXDEPLOY="$CACHE_DIR/linuxdeploy-x86_64.AppImage"
LINUXDEPLOY_QT="$CACHE_DIR/linuxdeploy-plugin-qt-x86_64.AppImage"
download_cached "$LINUXDEPLOY"    "https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage"
download_cached "$LINUXDEPLOY_QT" "https://github.com/linuxdeploy/linuxdeploy-plugin-qt/releases/download/continuous/linuxdeploy-plugin-qt-x86_64.AppImage"
chmod +x "$LINUXDEPLOY" "$LINUXDEPLOY_QT"

# Desktop file (validator-clean)
cp -f "$PWD/ALNVIEW/images/open.png" "$PWD/alnview.png"
cat > alnview.desktop <<'DESK'
[Desktop Entry]
Type=Application
Name=ALNview
Exec=ALNview
Icon=alnview
Categories=Science;
Terminal=false
DESK

# Create AppDir (NO appimage yet; we’ll add fonts + AppRun first)
rm -rf AppDir
"$LINUXDEPLOY" \
  --appdir AppDir \
  --executable "$PWD/ALNVIEW/ALNview" \
  --desktop-file "$PWD/alnview.desktop" \
  --icon-file "$PWD/alnview.png" \
  --plugin qt

# -------------------- Bundle fonts + private fontconfig --------------------
# Note: Qt standard dialogs (e.g. QFileDialog) follow the user locale. If the host
# locale is Japanese/Chinese/Korean, their labels/buttons may require CJK glyphs.
# This AppImage keeps its font payload small by forcing Qt messages to English in
# AppRun, so bundled Latin fonts are sufficient and the bundle stays independent
# from host fonts.
# Make the AppImage self-contained for fonts:
#   - ship the actual TTF files needed by ALNview
#   - use a private fonts.conf that ignores host font directories
#   - generate a private font cache inside AppDir
#   - point Qt at the bundled font directory as a fallback
mkdir -p AppDir/usr/share/fonts/truetype/dejavu
cp -a \
  /usr/share/fonts/dejavu/DejaVuSans.ttf \
  /usr/share/fonts/dejavu/DejaVuSans-Bold.ttf \
  /usr/share/fonts/dejavu/DejaVuSans-Oblique.ttf \
  /usr/share/fonts/dejavu/DejaVuSans-BoldOblique.ttf \
  /usr/share/fonts/dejavu/DejaVuSansMono.ttf \
  /usr/share/fonts/dejavu/DejaVuSansMono-Bold.ttf \
  /usr/share/fonts/dejavu/DejaVuSansMono-Oblique.ttf \
  /usr/share/fonts/dejavu/DejaVuSansMono-BoldOblique.ttf \
  AppDir/usr/share/fonts/truetype/dejavu/

mkdir -p AppDir/etc/fonts/conf.d AppDir/var/cache/fontconfig

# Fully private fontconfig setup:
# - <reset-dirs/> prevents scanning host font directories.
# - prefix="relative" keeps everything inside AppDir without FONTCONFIG_SYSROOT.
# - aliases force generic monospace/sans-serif requests to our bundled fonts.
cat > AppDir/etc/fonts/fonts.conf <<'EOF_FONTS'
<?xml version="1.0"?>
<fontconfig>
  <description>ALNview bundled font configuration</description>

  <reset-dirs />
  <dir prefix="relative">../../usr/share/fonts</dir>
  <cachedir prefix="relative">../../var/cache/fontconfig</cachedir>

  <alias>
    <family>sans-serif</family>
    <prefer>
      <family>DejaVu Sans</family>
    </prefer>
  </alias>

  <alias>
    <family>sans</family>
    <prefer>
      <family>DejaVu Sans</family>
    </prefer>
  </alias>

  <alias>
    <family>monospace</family>
    <prefer>
      <family>DejaVu Sans Mono</family>
    </prefer>
  </alias>

  <alias>
    <family>mono</family>
    <prefer>
      <family>DejaVu Sans Mono</family>
    </prefer>
  </alias>

  <match target="pattern">
    <test qual="any" name="family">
      <string>Monaco</string>
    </test>
    <edit name="family" mode="assign" binding="strong">
      <string>DejaVu Sans Mono</string>
    </edit>
  </match>
</fontconfig>
EOF_FONTS

# Build a cache that belongs to the AppImage, not to the build host.
# HOME/XDG_CACHE_HOME are redirected so no host cache is consulted or polluted.
mkdir -p AppDir/tmp-home AppDir/var/cache/xdg
HOME="$PWD/AppDir/tmp-home" \
XDG_CACHE_HOME="$PWD/AppDir/var/cache/xdg" \
FONTCONFIG_PATH="$PWD/AppDir/etc/fonts" \
FONTCONFIG_FILE="$PWD/AppDir/etc/fonts/fonts.conf" \
fc-cache -fsv "$PWD/AppDir/usr/share/fonts"
rm -rf AppDir/tmp-home

# Licenses (optional)
mkdir -p AppDir/usr/share/licenses/ALNview
if [ -d /usr/share/licenses/dejavu-fonts-common ]; then
  cp -a /usr/share/licenses/dejavu-fonts-common/. AppDir/usr/share/licenses/ALNview/ || true
fi

# -------------------- Safe AppRun (never clobber the ELF binary) --------------------
# linuxdeploy may have created AppRun as a symlink; replace it atomically.
cat > AppDir/AppRun.__new <<'EOF_APPRUN'
#!/bin/sh
set -eu
APPDIR="${APPDIR:-$(cd "$(dirname "$0")" && pwd)}"

# Use only the bundled font configuration and bundled fonts.
export FONTCONFIG_PATH="$APPDIR/etc/fonts"
export FONTCONFIG_FILE="$APPDIR/etc/fonts/fonts.conf"
export QT_QPA_FONTDIR="$APPDIR/usr/share/fonts/truetype/dejavu"
export XDG_CACHE_HOME="$APPDIR/var/cache/xdg"

# Keep the user's locale untouched. Fonts are loaded explicitly by Qt at startup,
# so even the first custom Qt dialog does not depend on host font discovery.

exec "$APPDIR/usr/bin/ALNview" "$@"
EOF_APPRUN
chmod 0755 AppDir/AppRun.__new
mv -f AppDir/AppRun.__new AppDir/AppRun

# -------------------- appimagetool (cached) --------------------
APPIMAGETOOL="$CACHE_DIR/appimagetool-x86_64.AppImage"
download_cached "$APPIMAGETOOL" "https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage"
chmod +x "$APPIMAGETOOL"

"$APPIMAGETOOL" AppDir

# -------------------- Sanity checks --------------------
APPIMAGE="$(ls -1 ALNview-*.AppImage 2>/dev/null | head -n1 || true)"
if [ -n "$APPIMAGE" ]; then
  echo "Built AppImage: $APPIMAGE"
  ./$APPIMAGE --appimage-extract >/dev/null
  echo "ALNview binary type:"
  file squashfs-root/usr/bin/ALNview || true
  echo "Bundled fonts:"
  ls -1 squashfs-root/usr/share/fonts/truetype/dejavu/*.ttf || true
  echo "Bundled fontconfig:"
  sed -n '1,200p' squashfs-root/etc/fonts/fonts.conf || true
fi
