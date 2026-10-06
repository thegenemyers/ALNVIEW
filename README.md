
# ALNview: A Qt-based alignment viewer for FASTGA

## _Author:  Gene Myers_
## _First:   April 10, 2025_

In brief, ALNview allows you to view the alignments in a .1aln as a "dot" plot between the two source
genomes, where each alignnment is viewed as a line segment.  You can view several different .1aln files that
are between the same two genomes as a set of layers that can be turned on an off.  For each layers
the thickness and color of the lines is under your control.  There is also a special layer that show
a true k-mer dot plot where you can control the size of k and the color of the dots.  This special
layer only becomes visible when the field of view in both dimensions is less than 1Mbp.

You can zoom by selecting regions or pressing up/down buttons. You can also pick alignment segments
which displays the coordinates, length, and iid of the alignment and gives you the option of requesting
to see the actual alignment in a secondary window.  And more ...

## Installation

Choose your platform: [Linux](#linux), [Windows](#windows-wsl-2), or
[macOS](#macos).

### Linux

Building from source requires Qt **6.9 or newer**, a C/C++ compiler, CMake,
and zlib development headers. The build and application startup have been
verified on Ubuntu 26.04 with Qt 6.10.

On Ubuntu 26.04, or another Ubuntu/Debian release providing Qt 6.9 or newer:

```sh
sudo apt-get update
sudo apt-get install git build-essential cmake qt6-base-dev zlib1g-dev
git clone https://github.com/Malfoy/ALNVIEW.git
cd ALNVIEW
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
./build/ALNview
```

Older distributions may provide an earlier Qt version. In that case, install
Qt 6.9 or newer separately and pass its installation prefix to CMake with
`-DCMAKE_PREFIX_PATH=/path/to/Qt/6.x/gcc_64`.

### Windows (WSL 2)

There is currently no native Windows installer. The source depends on Unix
APIs, so the Windows route is to build and run the Linux application in
WSL 2. This route has not yet been tested with ALNview.

WSL GUI applications require Windows 11 or Windows 10 build 19044 or newer.
See Microsoft's [WSL installation guide](https://learn.microsoft.com/en-us/windows/wsl/install)
and [GUI application guide](https://learn.microsoft.com/en-us/windows/wsl/tutorials/gui-apps).

Open PowerShell as Administrator and run:

```powershell
wsl --install
wsl --update
```

Restart if prompted, then open Ubuntu from the Start menu and complete its
initial setup. Follow the [Linux instructions](#linux) inside the Ubuntu
terminal, including the Qt 6.9 minimum version requirement. The default WSL
Ubuntu release may have an older Qt; use a newer distribution or install a
newer Qt as described above. Run `./build/ALNview` in that terminal to open
the application on the Windows desktop.

Windows files are accessible inside WSL under `/mnt/c`, for example
`/mnt/c/Users/your-name/Documents/alignments`.

### macOS

For the prebuilt application, download [ALNview.dmg](ALNview.dmg), open the
disk image, and copy ALNview to your Applications folder. Launch it from
Applications. The bundled application was supplied by the original author;
it has not been retested as part of the Linux build changes.

To build from source instead, install Apple's command-line tools with
`xcode-select --install`. With [Homebrew](https://brew.sh/) installed, run:

```sh
brew install cmake qt zlib
git clone https://github.com/Malfoy/ALNVIEW.git
cd ALNVIEW
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$(brew --prefix qt)" \
  -DZLIB_ROOT="$(brew --prefix zlib)"
cmake --build build -j 4
./build/ALNview
```

Use Qt 6.9 or newer; the current [Homebrew Qt package](https://formulae.brew.sh/formula/qt)
provides Qt 6. The macOS source build instructions have not been tested here.

## Opening alignments

The application opens a dialog for selecting a FASTGA `.1aln` alignment file.
Keep the associated source genome files available to enable sequence alignment
and k-mer dot plots.

## Alternative build and installation

On Linux or macOS, you can also use Qt's qmake in a separate directory:

```sh
mkdir -p build/qmake
cd build/qmake
qmake6 ../../viewer.pro
make -j 4
```

On macOS, Qt may provide `qmake` instead of `qmake6`; use the executable from
your Qt 6 installation. The qmake macOS build may produce an application bundle.

To install the CMake build to a chosen prefix, run
`cmake --install build --prefix /your/install/prefix`.
