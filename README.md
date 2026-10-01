
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

The prebuilt `.dmg` is for Apple computers. Linux users can build the source
with Qt 6.9.0 or higher, a C/C++ compiler, CMake, and zlib development headers.

## Build and run on Linux

On Ubuntu or Debian with Qt 6.9 or newer in their repositories:

```sh
sudo apt-get install build-essential cmake qt6-base-dev zlib1g-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
./build/ALNview
```

The application opens a dialog for selecting a FASTGA `.1aln` alignment file.
Keep the associated source genome files available to enable sequence alignment
and k-mer dot plots.

Alternatively, build with Qt's qmake in a separate directory:

```sh
mkdir -p build/qmake
cd build/qmake
qmake6 ../../viewer.pro
make -j 4
./ALNview
```

To install the CMake build to a chosen prefix, run
`cmake --install build --prefix /your/install/prefix`.
