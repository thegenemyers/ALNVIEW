
# ALNview: A Qt-based alignment viewer for FASTGA

## _Author:  Gene Myers_
## _First:   April 10, 2025_

In brief, ALNview allows you to view the alignments in a .1aln as a "dot" plot between the two source
genomes, where each align is an alignment segment.  You can view several different .1aln files that
are between the same two genomes as a set of layers that can be turned on an off.  For each layers
the thickness and color of the lines is under your control.  There is also a special layer that show
a true k-mer dot plot where you can control the size of k and the color of the dots.  This special
layer only becomes visible when the field of view in both dimentions is less than 1Mbp.

You can zoom be selecting regions or pressing up/down buttons. You can also pick alignment segments
which displays the coordinates, length, and iid of the alignment and gives you the option of requesting
to see the actual alignment in a secondary window.  And more ...

ALNview is currently only available as a prebuilt, binary .dmg for Apple computers.  We also give
you all the source files so the ambitious (or desperate :-) ) user can build it for other operating
systems using Qt 6.9.0 or higher.  Indeed if you make a binary image for a Windows or Unix machine
I would be very happy if you made a pull request with the relevant file and I will place it here.
