#ifndef STICKS
#define STICKS

#include "gene_core.h"
#include "GDB.h"
#include "hash.h"
#include "ONElib.h"


  //  C-structs for frame, view, and focus

typedef struct
  { double x, y;
    double w, h;
  } Frame;

typedef struct
  { int64 x, y;
    int64 w, h;
  } View;

typedef struct
  { int64 x, y;
  } Focus;


  //  Routine for coordinate display

#define FORMAT_n   0
#define FORMAT_c   1
#define FORMAT_s   2
#define FORMAT_s_c 3
#define FORMAT_i   4
#define FORMAT_i_c 5

int64 divide_bar(int64 t);

int64 digits(int64 t, char **suf, int *prec);

char *Map_Coord(GDB *gdb, int64 coord, int64 coord2, int format, int64 width);


  //  Data structures and routines for Quad Trees

#define BLK_SIZE 100000   //  Quad tree blocks ~4.6MB

typedef struct
  { int   length;
    int   depth;
    int   idx[8];
  } QuadLeaf;

typedef struct _qnode
  { int            length;
    int            depth;
    struct _qnode *quads[4];
  } QuadNode;


  //  Data structures and routines for managing a "plot" of layers

#define MAX_LAYERS 5

typedef struct
  { int64 abeg, aend;
    int64 bbeg, bend;
    int16 iid;
    int16 mark;
    int   idx;
  } DotSegment;

typedef struct
  { int         nref;
    char       *name;
    OneFile    *input;
    int64       novls;
    int         tspace;
    DotSegment *segs;
    QuadNode   *qtree;
    QuadNode   *blocks;
  } DotLayer;

typedef struct
  { int        nref;
    GDB        gdb;
    char      *hash;
    char      *name;
  } DotGDB;

typedef struct
  { int64        alen;
    int64        blen;
    DotGDB      *db1;
    DotGDB      *db2;
    int          nlays;
    DotLayer    *layers[MAX_LAYERS];
    void        *dotmemory;
  } DotPlot;

DotPlot *createPlot(char *alnPath, int lCut, int iCut, int sCut, DotPlot *plot);

DotPlot *copyPlot(DotPlot *plot);

QuadLeaf *Plot_Layer(DotPlot *plot, int ilay, Frame *query);

void Free_List(QuadLeaf *list);

void Free_DotPlot(DotPlot *plot);


  //  Data structures and routines to support alignment generation & display

char *create_alignment(DotPlot *plot, DotLayer *layer, DotSegment *seg, char **title);

#endif
