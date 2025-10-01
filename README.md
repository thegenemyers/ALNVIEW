
# ALNview: A Qt-based alignment viewer for FASTGA

## _Author:  Gene Myers_
## _First:   April 10, 2025_

In brief, DaViewer allows you to view any subset of the information in a given .las
file of local alignments computed by the daligner or the forth coming damapper, and any
track information associated with the read database(s) that were compared to produce
the local alignments (LAs).  The complete document can be found at the
[Dazzler blog](https://dazzlerblog.wordpress.com/command-guides/daviewer).

To build DaViewer you will need to download the latest (free) version of the Qt library
which you can get here.  Once Qt is installed you will need to set the shell variable
QTDIR to the path of the library which depends on where you chose to install it.  I set
this variable in the shell file .login at my root so that every time I open a window,
Qt is known.  Another good spot to set this shell variable is .cshrc (for bash) of your
shell's equivalent.

To make the viewer you must first utter "qmake viewer.pro" in the directory to which
you downloaded the module.  qmake is a Qt program that constructs a Makefile (much like
configure) from the specification given in viewer.pro.  You then simply say "make" and
hopefully the compiles all go off without complaint.  You should then find the DaViewer
application in the folder, which you can launch by double clicking it.
