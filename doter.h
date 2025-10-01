#ifndef DOTER
#define DOTER

#include "GDB.h"
#include "sticks.h"

#define MAX_DOTPLOT 1000000

void *dotplot_memory();

typedef struct 
  { uint64 code;
    int    pos; 
  } Tuple;

typedef struct
  { int    ahit, brun;
    int   *aplot;
    Tuple *blist;
  } Dots;

Dots *dotplot(DotPlot *plot, int kmer, View *view);

// Dots *dotplot(DotPlot *plot, int kmer, View *view, int rectW, int rectH, uint8 **raster);

#endif
