#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "GDB.h"
#include "doter.h"
#include "sticks.h"

#undef  DEBUG_CHECK
#undef  DEBUG_STATS
#undef  DEBUG_SHOW
#undef  DEBUG_FILL

static int build_vector(int64 len, char *seq, int kmer, Tuple *list)
{ uint64 Kmask, Cumber[5];
  int64  p, q, k, km1;
  uint64 u, c, x;

  km1  = kmer-1;
  len -= km1;
 
  if (kmer == 32)
    Kmask = 0xffffffffffffffffllu;
  else
    Kmask = (0x1llu << 2*kmer) - 1;

  Cumber[0] = (0x3llu << 2*km1);
  Cumber[1] = (0x2llu << 2*km1);
  Cumber[2] = (0x1llu << 2*km1);
  Cumber[3] = 0x0llu;
  Cumber[4] = 0x0llu;

  c = u = 0;
  for (p = 0; p < km1; p++)
    { x = seq[p];
      c = (c << 2) | x;
      u = (u >> 2) | Cumber[x];
    }
  q = 0;
  k = 0;
  seq += km1;
  for (p = 0; p < len; p++)
    { x = seq[p];
      if (x >= 4)
        { k = p+kmer;
          c = u = 0;
        }
      else
        { c = ((c << 2) | x) & Kmask;
          u = (u >> 2) | Cumber[x];
          if (p >= k)
            { if (u < c)
                list[q].code = u;
              else
                list[q].code = c;
              list[q++].pos = p;
            }
        }
    }
  return ((int) q);
}

static int merge(int brun, Tuple *blist, int arun, Tuple *alist)
{ int    *aplot;
  int     i, j;
  int     al, bl;
  int     al2, bl2;
  int     y;
  uint64  kb, lc, anull;

  bl = brun;
  al = arun;

  lc = alist[al-1].code;
  for (bl2 = bl-1; bl2 >= 0; bl2--)
    if (blist[bl2].code < lc)
      break;
  bl2 += 1;
  for (al2 = al-2; al2 >= 0; al2--)
    if (alist[al2].code < lc)
      break;
  al2 += 1;

  aplot = (int *) alist;
  anull = 2*al-1;

  y = 0;
  i = j = 0;
  while (i < bl2)
    { kb = blist[i].code;
      while (alist[j].code < kb)
        j += 1;

      if (alist[j].code == kb)
        { blist[i++].code = y;
          while (blist[i].code == kb)
            blist[i++].code = y;
          aplot[y++] = alist[j++].pos;
          while (alist[j].code == kb)
            aplot[y++] = alist[j++].pos;
          aplot[y++] = -1;
        }
      else
        { blist[i++].code = anull;
          while (blist[i].code == kb)
            blist[i++].code = anull;
        }
    }

  if (bl2 < bl && blist[i].code == lc)
    { blist[i++].code = y;
      while (i < bl && blist[i].code == lc)
        blist[i++].code = y;
      while (al2 < al)
        aplot[y++] = alist[al2++].pos;
      aplot[y++] = -1;
    }
  while (i < bl)
    blist[i++].code = anull;
  aplot[anull] = -1;

  return (y);
}

static int TSORT(const void *l, const void *r)
{ Tuple *x = (Tuple *) l;
  Tuple *y = (Tuple *) r;
  if (x->code < y->code)
    return (-1);
  if (x->code > y->code)
    return (1);
  return (0);
}

static void map(GDB *gdb, int64 coord, int *cps, int *pos)
{ GDB_SCAFFOLD *scf;
  GDB_CONTIG   *ctg;
  int           i, s, c;

  scf = gdb->scaffolds;
  ctg = gdb->contigs;

  for (i = 0; i < gdb->nscaff; i++)       //  binary search instead ?
    if (ctg[scf[i].fctg].sbeg > coord)
      break;
  s = i-1;

  for (i = scf[s].fctg; i < scf[s].ectg; i++)
    if (ctg[i].sbeg > coord)
      break;
  c = i-1;

  *cps = c;
  *pos = coord - ctg[c].sbeg;
}

static int runOfN(int len, char *seq)
{ int i;

  for (i = 0; i < len; i++)
    seq[i] = 4;
  return (len);
}

static char *build_string(GDB *gdb, int64 vbeg, int64 vend, char *seq)
{ GDB_CONTIG *ctg;
  int cpb, cpe;
  int beg, end;
  int s, c, x;

  ctg = gdb->contigs;

  map(gdb,vbeg,&cpb,&beg);
  map(gdb,vend,&cpe,&end);

  // printf(" %d,%d - %d,%d\n",cpb,beg,cpe,end);

  if (cpb == cpe)
    seq = Get_Contig_Piece(gdb,cpb,beg,end,NUMERIC,seq);
  else
    { seq = Get_Contig_Piece(gdb,cpb,beg,ctg[cpb].clen,NUMERIC,seq);
      s = ctg[cpb].scaf;
      x = ctg[cpb].clen - beg;
      for (c = cpb+1; c < cpe; c++)
        { if (ctg[c].scaf == s)
            x += runOfN(ctg[c].sbeg-(ctg[c-1].sbeg+ctg[c-1].clen),seq+x);
          Get_Contig(gdb,c,NUMERIC,seq+x);
          s = ctg[c].scaf;
          x += ctg[c].clen;
        }
      if (ctg[cpe].scaf == s)
        x += runOfN(ctg[cpe].sbeg-(ctg[cpe-1].sbeg+ctg[cpe-1].clen),seq+x);
      Get_Contig_Piece(gdb,c,0,end,NUMERIC,seq+x);
    }
  return (seq);
}

void *dotplot_memory()
{ return (malloc((sizeof(Tuple)+1)*2*MAX_DOTPLOT + 8 + sizeof(Dots))); }

Dots *dotplot(DotPlot *plot, int kmer, View *view)
{ Dots  *dot   = (Dots *) plot->dotmemory;
  Tuple *alist = (Tuple *) (dot+1);
  Tuple *blist = alist + MAX_DOTPLOT;
  char  *aseq  = (char *) (blist+MAX_DOTPLOT);
  char  *bseq  = aseq + (MAX_DOTPLOT+4);
  int   *aplot = (int *) alist;

  int    arun, brun, ahit;

  int64 vX = view->x;
  int64 vY = view->y;
  int64 vW = view->w;
  int64 vH = view->h;

  // double xa = (rectH-44.)/vH;
  // double xb = 22.;
  // double ya = (rectW-44.)/vW;
  // double yb = 22.;

  // printf(" %lld-%lld vs %lld-%lld %d\n",vX,vX+vW,vY,vY+vH,kmer);

  aseq = build_string(&(plot->db1->gdb),vX,vX+vW,aseq);
  bseq = build_string(&(plot->db2->gdb),vY,vY+vH,bseq);

  arun = build_vector(vW,aseq,kmer,alist);
  brun = build_vector(vH,bseq,kmer,blist);

  qsort(alist,arun,sizeof(Tuple),TSORT);
  qsort(blist,brun,sizeof(Tuple),TSORT);

#ifdef DEBUG_CHECK
  { int i;

    for (i = 1; i < arun; i++)
      if (alist[i].code < alist[i-1].code)
        printf("Not sorted\n");

    for (i = 1; i < brun; i++)
      if (blist[i].code < blist[i-1].code)
        printf("Not sorted\n");
  }
#endif

  ahit = merge(brun,blist,arun,alist);

#ifdef DEBUG_SHOW
  { int i, j;

    printf("Scan Lines:\n");
    for (i = 0; i < brun; i++)
      { printf(" %5d:",blist[i].pos);
        j = blist[i].code;
        while (aplot[j] >= 0)
          { printf(" %5d",aplot[j]);
            j += 1;
          }
        printf("\n");
      }
  }
#endif

#ifdef DEBUG_CHECK
  { int i, j;
    int amax = vW - (kmer-1);

    for (i = 0; i < brun; i++)
      { j = blist[i].code;
        if (j >= ahit && j != 2*arun-1)
          printf("bcode val out of bounds\n");
        while (aplot[j] >= 0)
          { if (aplot[j] > amax)
              printf("aplot val out of bounds\n");
            j += 1;
          }
        if (aplot[j] != -1)
          printf("aplot terminator is not -1\n");
      }
  }
#endif
  
#ifdef DEBUG_STATS
  { int i, j;
    int nz, nel;

    printf("Scan Stats:\n");
    nz = nel = 0;
    for (i = 0; i < brun; i++)
      { j = blist[i].code;
        if (aplot[j] >= 0)
          { nz += 1;
            while (aplot[j] >= 0)
              { nel += 1;
                j += 1;
              }
          }
      }

    printf("  Non-Zero lines: %d (out of %d)\n",nz,blen);
    printf("  Av/line = %.1f (out of %d)\n",(1.*nel)/alen,alen);
    printf("  Density = %.3f%%\n",((100.*nel)/alen)/blen);
  }
#endif

  dot->ahit  = ahit;
  dot->brun  = brun;
  dot->aplot = aplot;
  dot->blist = blist;

/*
  // printf("Paint = (%lld,%lld) %lld x %lld into %d x %d\n",vX,vY,vW,vH,rectW,rectH);

  { int i, x;

    for (i = 0; i < ahit; i++)
      { x = aplot[i];
        if (x >= 0)
          aplot[i] = ((int) floor(xa*x+xb));
      }
  }

  { int i, k, u;
    uint8 *ras;

    for (i = 0; i < brun; i++)
      { ras = raster[((int) floor(ya*blist[i].pos+yb))];
        k = blist[i].code;
        while (1)
          { u = aplot[k++];
            if (u < 0) 
              break;
            ras[u] = 255;
          }
      }
  }
*/

  return (dot);
}
