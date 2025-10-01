#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>

#include "hash.h"
#include "sticks.h"
#include "doter.h"
#include "ONElib.h"
#include "alncode.h"
#include "gene_core.h"
#include "GDB.h"

#undef DEBUG_ADD
#undef DEBUG_FIND
#undef DEBUG_LAYER



/*******************************************************************************************
*
*   COORDINATE DISPLAY
*
*******************************************************************************************/

static int64 units[] = { 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000,
                       100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000, 20000000,
                       50000000, 100000000, 200000000, 500000000, 1000000000, 2000000000,
                       5000000000, 10000000000 };

static int64 den3[] = { 0, 1, 1, 1, 1000, 1000, 1000, 1000000, 1000000, 1000000,  1000000000,
                        1000000000, 1000000000, 1000000000000 };

static char* suf3[] = { "", "", "", "", "k", "k", "k", "M", "M", "M", "G", "G", "G", "T" };

static int   pre3[] = {  0,  0,  0,  0,   2,   1,   0,   2,   1,   0,   2,   1,   0,  2  };

int64 digits(int64 max, char **suf, int *prec)
{ int   u;
  int64 w;

  w = 1;
  for (u = 1; u < 31; u++)
    { w *= 10;
      if (w > max)
        break;
    }
  *suf  = suf3[u];
  *prec = pre3[u];
  return (den3[u]);
}

int64 divide_bar(int64 t)
{ int u;

  for (u = 1; u < 31; u++)
    { if (units[u] > t)
        break;
    }
  u -= 1;
  return (units[u]);
}


static void map(GDB *gdb, int format, int64 coord, int *sp, int *cp)
{ GDB_SCAFFOLD *scf;
  GDB_CONTIG   *ctg;
  int           i, s, c;

  scf = gdb->scaffolds;
  ctg = gdb->contigs;

  for (i = 0; i < gdb->nscaff; i++)                         //  binary search instead ?
    if (ctg[scf[i].fctg].sbeg > coord)
      break;
  s = i-1;

  if (format % 2 != 0)
    { for (i = scf[s].fctg; i < scf[s].ectg; i++)
        if (ctg[i].sbeg > coord)
          break;
      c = i-1;
    }
  else
    c = scf[s].fctg;

  *sp = s;
  *cp = c;
}

static char *write_coord(GDB *gdb, int format, int s, int c, int64 coord,
                         double den, int prec, char *suf, char *answer)
{ GDB_SCAFFOLD *scf;
  GDB_CONTIG   *ctg;
  double        den2;
  int           prec2;
  char         *suf2;

  scf = gdb->scaffolds;
  ctg = gdb->contigs;

  if (format % 2 == 0)
    { if (format > 0)
        if (s < 0)
          den2 = digits(scf[-(s+1)].slen,&suf2,&prec2);
        else
          den2 = digits(scf[s].slen,&suf2,&prec2);
      else
        den2 = den;
    }
  else
    if (c < 0)
      den2 = digits(ctg[-(c+1)].clen,&suf2,&prec2);
    else
      den2 = digits(ctg[c].clen,&suf2,&prec2);

  if (den2 < den)
    { den  = den2;
      prec = prec2;
      suf  = suf2;
    }

  switch (format)
  { case 0:
      sprintf(answer,"%.*f%s",prec,coord/den,suf);
      break;
    case 1:
      if (c < 0)
        sprintf(answer,"%.*f%s",prec,(coord-ctg[-(c+1)].sbeg)/den,suf);
      else
        sprintf(answer,".%d:%.*f%s",c+1,prec,(coord-ctg[c].sbeg)/den,suf);
      break;
    case 2:
      if (s < 0)
        sprintf(answer,"%.*f%s",prec,(coord-ctg[c].sbeg)/den,suf);
      else
        sprintf(answer,"@%d:%.*f%s",s+1,prec,(coord-ctg[c].sbeg)/den,suf);
      break;
    case 3:
      if (s < 0)
        if (c < 0)
          sprintf(answer,"%.*f%s",prec,(coord-ctg[-(c+1)].sbeg)/den,suf);
        else
          sprintf(answer,"%d:%.*f%s",(c-scf[-(s+1)].fctg)+1,prec,(coord-ctg[c].sbeg)/den,suf);
      else
        sprintf(answer,"@%d.%d:%.*f%s",s+1,(c-scf[s].fctg)+1,prec,(coord-ctg[c].sbeg)/den,suf);
      break;
    case 4:
      if (s < 0)
        sprintf(answer,"%.*f%s",prec,(coord-ctg[c].sbeg)/den,suf);
      else
        sprintf(answer,"@%s:%.*f%s",gdb->headers+scf[s].hoff,prec,(coord-ctg[c].sbeg)/den,suf);
      break;
    case 5:
      if (s < 0)
        if (c < 0)
          sprintf(answer,"%.*f%s",prec,(coord-ctg[-(c+1)].sbeg)/den,suf);
        else
          sprintf(answer,"%d:%.*f%s",(c-scf[-(s+1)].fctg)+1,prec,(coord-ctg[c].sbeg)/den,suf);
      else
        sprintf(answer,"@%s.%d:%.*f%s",gdb->headers+scf[s].hoff,(c-scf[s].fctg)+1,
                                       prec,(coord-ctg[c].sbeg)/den,suf);
      break;
    default:
      sprintf(answer,"NA ??");
      break;
  }

  return (answer);
}


char *Map_Coord(GDB *gdb, int64 coord1, int64 coord2, int format, int64 width)
{ static char answer1[5000];
  static char answer2[5000];

  double den;
  int    prec;
  char  *suf;
  int    s1, c1;
  int    s2, c2;

  if (coord1 >= 0)
    map(gdb,format,coord1,&s1,&c1);
  if (coord2 >= 0)
    map(gdb,format,coord2,&s2,&c2);

  den = digits(width,&suf,&prec);

  if (coord1 < 0)
    { write_coord(gdb,format,s2,c2,coord2,den,prec,suf,answer2);
      return (answer2);
    }

  write_coord(gdb,format,s1,c1,coord1,den,prec,suf,answer1);
  if (coord2 < 0)
    return (answer1);

  sprintf(answer1+strlen(answer1)," - ");

  if (format > 1)
    { if (s1 == s2)
        s2 = -(s2+1);
    }
  if (format%2 == 1)
    { if (c1 == c2)
        c2 = -(c2+1);
    }
  write_coord(gdb,format,s2,c2,coord2,den,prec,suf,answer1+strlen(answer1));
  return (answer1);
}
 

/*******************************************************************************************
*
*   QUAD TREE
*
*******************************************************************************************/

#define QNW 0
#define QNE 1
#define QSE 2
#define QSW 3

static DotSegment *SEGS;
static int         FREECNT;
static QuadNode   *BLOCKS;

typedef struct
  { double abeg, aend;
    double bbeg, bend;
  } Double_Box;

static QuadNode *New_Quad()
{ if (FREECNT >= BLK_SIZE)
    { QuadNode *block;

      block = malloc(sizeof(QuadNode)*(BLK_SIZE+1)); 
      block[BLK_SIZE].quads[0] = BLOCKS;
      block[BLK_SIZE].length   = 0;
      BLOCKS  = block;
      FREECNT = 0;
    }
  return (BLOCKS+FREECNT++);
}

static int BEG_QUAD(Double_Box *seg, double amid, double bmid)
{ if (seg->abeg < amid || (seg->abeg == amid && seg->aend <= amid))
    if (seg->bbeg < bmid || (seg->bbeg == bmid && seg->bend <= bmid))
      return (QNW);
    else
      return (QNE);
  else
    if (seg->bbeg < bmid || (seg->bbeg == bmid && seg->bend <= bmid))
      return (QSW);
    else
      return (QSE);
}

static int END_QUAD(Double_Box *seg, double amid, double bmid)
{ if (seg->aend < amid || (seg->aend == amid && seg->abeg <= amid))
    if (seg->bend < bmid || (seg->bend == bmid && seg->bbeg <= bmid))
      return (QNW);
    else
      return (QNE);
  else
    if (seg->bend < bmid || (seg->bend == bmid && seg->bbeg <= bmid))
      return (QSW);
    else
      return (QSE);
}

static int WCH_QUAD(double ac, double bc, double amid, double bmid)
{ if (ac < amid)
    if (bc < bmid)
      return (QNW);
    else
      return (QNE);
  else
    if (bc < bmid)
      return (QSW);
    else
      return (QSE);
}

static void QUAD_CUT(Double_Box *frame, double amid, double bmid, int quad)
{ if (quad < 2)
    frame->aend = amid;
  else
    frame->abeg = amid;
  if (quad % 3 == 0)
    frame->bend = bmid;
  else
    frame->bbeg = bmid;
}

static void Clip_Segment(Double_Box *seg, Double_Box *frame)
{ double t, x1, x2, y1, y2;
  int    flipx, flipy, inter;

  inter = 0;
  flipx = (seg->abeg > seg->aend);
  if (flipx)
    { x1 = seg->aend; x2 = seg->abeg;
      y1 = seg->bend; y2 = seg->bbeg;
    }
  else
    { x1 = seg->abeg; x2 = seg->aend;
      y1 = seg->bbeg; y2 = seg->bend;
    }
  flipy = (y1 > y2);
  if (flipy)
    { t = x1; x1 = x2; x2 = t;
      t = y1; y1 = y2; y2 = t;
    }
  if (y2 > frame->bend)
    { x2 = x1 + (x2 - x1) * (frame->bend - y1) / (y2 - y1);
      y2 = frame->bend;
      inter = 1;
    }
  if (y1 < frame->bbeg)
    { x1 = x1 + (x2 - x1) * (frame->bbeg - y1) / (y2 - y1);
      y1 = frame->bbeg;
      inter = 1;
    }
  if (flipy)
    { t = x1; x1 = x2; x2 = t;
      t = y1; y1 = y2; y2 = t;
    }
  if (x2 > frame->aend)
    { y2 = y1 + (y2 - y1) * (frame->aend - x1) / (x2 - x1);
      x2 = frame->aend;
      inter = 1;
    }
  if (x1 < frame->abeg)
    { y1 = y1 + (y2 - y1) * (frame->abeg - x1) / (x2 - x1);
      x1 = frame->abeg;
      inter = 1;
    }
  if (inter)
    { if (flipx)
        { seg->abeg = x2; seg->aend = x1;
          seg->bbeg = y2; seg->bend = y1;
        }
      else
        { seg->abeg = x1; seg->aend = x2;
          seg->bbeg = y1; seg->bend = y2;
        }
    }
}

static QuadNode *Add_To_Node(QuadNode *quad, Double_Box *frame, Double_Box *seg, int idx, int deep)
{ int    qb, qe;
  double amid, bmid;

#ifdef DEBUG_ADD
  printf("%*s(%.0f,%.0f)->(%.0f,%.0f) in box %.0f-%.0f vs %.0f-%.0f\n",2*deep,"",
         seg->abeg,seg->bbeg,seg->aend,seg->bend,
         frame->abeg,frame->aend,frame->bbeg,frame->bend); fflush(stdout);
  if (seg->abeg < frame->abeg || seg->aend > frame->aend ||
      seg->bbeg < frame->bbeg || seg->bend > frame->bend)
    printf("Not in Frame\n");
#endif

  if (quad == NULL)
    { quad = New_Quad();
      quad->length = 1;
      quad->depth  = deep;
      ((QuadLeaf *) quad)->idx[0] = idx;
#ifdef DEBUG_ADD
      printf("%*sSimple Add %d\n",2*deep,"",quad->length); fflush(stdout);
#endif
      return (quad);
    }

  if (quad->length >= 8)
    { int         i;
      QuadLeaf    leaf;
      Double_Box  new_frame;
      DotSegment *o;

#ifdef DEBUG_ADD
      printf("%*sOverfull\n",2*deep,""); fflush(stdout);
#endif
      leaf = *((QuadLeaf *) quad);
      quad->length = 0;
      for (i = 0; i < 4; i++)
        quad->quads[i] = NULL;
      new_frame = *frame;
      Add_To_Node(quad,&new_frame,seg,idx,deep);   // quad already split
      for (i = 0; i < leaf.length; i++)
        { new_frame = *frame;
          o = SEGS + leaf.idx[i];
          seg->abeg = o->abeg;
          seg->bbeg = o->bbeg;
          seg->aend = o->aend;
          seg->bend = o->bend;
          Clip_Segment(seg,frame);
          Add_To_Node(quad,&new_frame,seg,leaf.idx[i],deep);  // quad already split
        }
      return (quad);
    }

  if (quad->length > 0)
    { ((QuadLeaf *) quad)->idx[quad->length++] = idx;
#ifdef DEBUG_ADD
      printf("%*sSimple Add %d\n",2*deep,"",quad->length); fflush(stdout);
#endif
      return (quad);
    }

  amid = (frame->abeg+frame->aend)/2.;
  bmid = (frame->bbeg+frame->bend)/2.;
  qb = BEG_QUAD(seg,amid,bmid);
  qe = END_QUAD(seg,amid,bmid);
#ifdef DEBUG_ADD
  printf("%*s%.0f x %.0f -> %d %d\n",2*deep,"",amid,bmid,qb,qe); fflush(stdout);
#endif

  deep += 1;

  if (qb == qe)
    { QUAD_CUT(frame,amid,bmid,qb);
#ifdef DEBUG_ADD
      printf("%*sEasy 1\n",2*deep-2,""); fflush(stdout);
#endif
      quad->quads[qb] = Add_To_Node(quad->quads[qb],frame,seg,idx,deep);
      return (quad);
    }

  { Double_Box seg2, frame2;
    double     x, y;

    frame2 = *frame;
    seg2   = *seg;

    if (abs(qb-qe) % 2 == 1)
      { if (qb+qe == 3)
          { x  = (amid - seg->abeg) / (seg->aend-seg->abeg);
            seg->bend = seg2.bbeg = seg->bbeg + x*(seg->bend-seg->bbeg);
            seg->aend = seg2.abeg = amid;
          }
        else
          { x  = (bmid - seg->bbeg) / (seg->bend-seg->bbeg);
            seg->aend = seg2.abeg = seg->abeg + x*(seg->aend-seg->abeg);
            seg->bend = seg2.bbeg = bmid;
          }
#ifdef DEBUG_ADD
        printf("%*sModr 2\n",2*deep-2,""); fflush(stdout);
#endif
        QUAD_CUT(frame,amid,bmid,qb);
        QUAD_CUT(&frame2,amid,bmid,qe);
        quad->quads[qb] = Add_To_Node(quad->quads[qb],frame,seg,idx,deep);
        quad->quads[qe] = Add_To_Node(quad->quads[qe],&frame2,&seg2,idx,deep);
        return (quad);
      }

    x = (bmid - seg->bbeg) / (seg->bend-seg->bbeg);
    y = (amid - seg->abeg) / (seg->aend-seg->abeg);
    if (x == y)
      { seg->aend = seg2.abeg = amid;
        seg->bend = seg2.bbeg = seg->bbeg + x*(seg->bend-seg->bbeg);
#ifdef DEBUG_ADD
        printf("%*sX ???\n",2*deep-2,""); fflush(stdout);
#endif
        quad->quads[qb] = Add_To_Node(quad->quads[qb],frame,seg,idx,deep);
        quad->quads[qe] = Add_To_Node(quad->quads[qe],&frame2,&seg2,idx,deep);
        return (quad);
      }

    { Double_Box seg3, frame3;
      int        qm;

      frame3 = *frame;
      seg3   = *seg;
      if (x < y)
        { seg3.bend  = seg2.bbeg = seg->bbeg + y*(seg->bend-seg->bbeg);
          seg3.aend  = seg2.abeg = amid;
          seg->aend  = seg3.abeg = seg->abeg + x*(seg->aend-seg->abeg);
          seg->bend  = seg3.bbeg = bmid;
        }
      else
        { seg3.aend  = seg2.abeg = seg->abeg + x*(seg->aend-seg->abeg);
          seg3.bend  = seg2.bbeg = bmid;
          seg->bend  = seg3.bbeg = seg->bbeg + y*(seg->bend-seg->bbeg);
          seg->aend  = seg3.abeg = amid;
        }
      qm = WCH_QUAD((seg3.abeg+seg3.aend)/2.,(seg3.bbeg+seg3.bend)/2.,amid,bmid);
      QUAD_CUT(frame,amid,bmid,qb);
      QUAD_CUT(&frame2,amid,bmid,qe);
      QUAD_CUT(&frame3,amid,bmid,qm);
#ifdef DEBUG_ADD
      printf("%*sWorst 3 (%d) (%9g,%9g)\n",2*deep-2,"",qm,x,y); fflush(stdout);
#endif
      quad->quads[qb] = Add_To_Node(quad->quads[qb],frame,seg,idx,deep);
      quad->quads[qm] = Add_To_Node(quad->quads[qm],&frame3,&seg3,idx,deep);
      quad->quads[qe] = Add_To_Node(quad->quads[qe],&frame2,&seg2,idx,deep);
      return (quad);
    }
  }
}

static void Make_QuadTree(DotPlot *plot, int ilay)
{ QuadNode  *quad;
  Double_Box seg;
  Double_Box frame;
  int64      novl;
  int        i;

  SEGS    = plot->layers[ilay]->segs;
  BLOCKS  = NULL;
  FREECNT = BLK_SIZE;

  novl = plot->layers[ilay]->novls;
  quad = NULL;
  for (i = 0; i < novl; i++)
    { seg.abeg = SEGS[i].abeg;
      seg.aend = SEGS[i].aend;
      seg.bbeg = SEGS[i].bbeg;
      seg.bend = SEGS[i].bend;
      frame.abeg = 0.;
      frame.bbeg = 0.;
      frame.aend = plot->alen;
      frame.bend = plot->blen;
#ifdef DEBUG_ADD
      printf("Doing %d\n",i);
#endif
      quad = Add_To_Node(quad,&frame,&seg,i,0);
     }
  plot->layers[ilay]->qtree  = quad;
  plot->layers[ilay]->blocks = BLOCKS;
}

static char *QLabel[] = { "NW", "NE", "SE", "SW", " *" };

static void Show_QuadNode(QuadNode *quad, int dir, Double_Box *frame)
{ int i;

  printf("%*s %s:",2*quad->depth,"",QLabel[dir]);

  if (quad->length > 0)
    { for (i = 0; i < quad->length; i++)
        printf(" %d",((QuadLeaf *) quad)->idx[i]);
      printf("\n");
      fflush(stdout);
      return;
    }

  { Double_Box copy;
    double     amid, bmid;

    printf(" [%9g-%9g] x [%9g,%9g]\n",frame->abeg,frame->aend,frame->bbeg,frame->bend);
    fflush(stdout);
    amid = (frame->abeg+frame->aend)/2.;
    bmid = (frame->bbeg+frame->bend)/2.;
    for (i = 0; i < 4; i++)
      if (quad->quads[i] == NULL)
        printf("%*s %s: .\n",2*quad->depth+2,"",QLabel[i]);
      else
        { copy = *frame; 
          QUAD_CUT(&copy,amid,bmid,i);
          Show_QuadNode(quad->quads[i],i,&copy);
        }
  }
}

static void Show_QuadTree(DotPlot *plot, int ilay)
{ Double_Box frame;

  frame.abeg = 0.;
  frame.bbeg = 0.;
  frame.aend = plot->alen;
  frame.bend = plot->blen;

  if (plot->layers[ilay]->qtree == NULL)
    printf("Empty Tree\n");
  else
    Show_QuadNode(plot->layers[ilay]->qtree,4,&frame);
}

static int nquad;
static int nleaf;
static int nlists;
static int dhist[100];
static int phist[100];

static void Stat_QuadNode(QuadNode *quad)
{ int i;

  if (quad == NULL)
    return;

  nquad += 1;
  dhist[quad->depth] += 1;

  if (quad->length > 0)
    { nlists += quad->length;
      nleaf  += 1;
      for (i = 0; i < quad->length; i++)
        SEGS[((QuadLeaf *) quad)->idx[i]].mark += 1;
      return;
    }

  for (i = 0; i < 4; i++)
    Stat_QuadNode(quad->quads[i]);
}

static void Stat_QuadTree(DotPlot *plot, int ilay)
{ int64 novl;
  int   i;

  nquad = 0;
  nleaf = 0;
  nlists = 0;
  for (i = 0; i < 100; i++)
    dhist[i] = 0;

  SEGS = plot->layers[ilay]->segs;
  novl = plot->layers[ilay]->novls;

  Stat_QuadNode(plot->layers[ilay]->qtree);

  printf("\nQuad Stats:\n");
  printf("  %d nodes of which %d are leaves.\n",nquad,nleaf);
  printf("  An average of %.1f segs per leaf\n",(1.*nlists)/nleaf);
  printf("  An average of %.1f pieces per alignment\n",(1.*nlists)/novl);
  printf("  Occupies %lldMB of memory\n",(sizeof(QuadNode)*((int64) nquad))/(1<<20));

  printf("\nDepth Profile:\n");
  for (i = 99; i >= 0; i--)
    if (dhist[i] > 0)
      printf(" %2d: %8d\n",i,dhist[i]);
  fflush(stdout);

  for (i = 0; i < novl; i++)
    { if (SEGS[i].mark >= 100)
        phist[99] += 1;
      else
        phist[SEGS[i].mark] += 1;
      SEGS[i].mark = 0;
    }
  
  printf("\nSegement Fracture Profile:\n");
  for (i = 99; i >= 0; i--)
    if (phist[i] > 0)
      printf(" %2d: %8d\n",i,phist[i]);
  fflush(stdout);
}

static QuadLeaf *LIST;

void QuadNode_Find(QuadNode *quad, Double_Box *frame, Double_Box *query)
{ 
#ifdef DEBUG_FIND
  printf(" [%.0f-%.0f] x [%.0f-%.0f]\n",frame->abeg,frame->aend,frame->bbeg,frame->bend);
#endif

  if (quad == NULL)
    return;

  if (quad->length > 0)
    { int i, id;

#ifdef DEBUG_FIND
      printf("%*sLeaf:",2*quad->depth,"");
#endif
      for (i = 0; i < quad->length; i++)
        { id = ((QuadLeaf *) quad)->idx[i];
#ifdef DEBUG_FIND
          printf(" %d%s",id,SEGS[id].mark?"*":"");
#endif
          if (SEGS[id].mark == 0)
            { SEGS[id].mark = 1;
              LIST->idx[LIST->length++] = ((QuadLeaf *) quad)->idx[i];
              if (LIST->length >= 8)
                { LIST = (QuadLeaf *) New_Quad();
                  LIST->length = 0;
                }
            }
        }
#ifdef DEBUG_FIND
      printf("\n");
#endif
      return;
    }

  { double amid, bmid;
    double atmp, btmp;

    amid = (frame->abeg + frame->aend) / 2.;
    bmid = (frame->bbeg + frame->bend) / 2.;

    if (query->abeg < amid && query->bbeg < bmid)
      { atmp = frame->aend;
        btmp = frame->bend;
        frame->aend = amid;
        frame->bend = bmid;
#ifdef DEBUG_FIND
        printf("%*sNW:",2*quad->depth+2,"");
#endif
        QuadNode_Find(quad->quads[0],frame,query);
        frame->aend = atmp;
        frame->bend = btmp;
      }

    if (query->abeg < amid && query->bend > bmid)
      { atmp = frame->aend;
        btmp = frame->bbeg;
        frame->aend = amid;
        frame->bbeg = bmid;
#ifdef DEBUG_FIND
        printf("%*sNE:",2*quad->depth+2,"");
#endif
        QuadNode_Find(quad->quads[1],frame,query);
        frame->aend = atmp;
        frame->bbeg = btmp;
      }

    if (query->aend > amid && query->bend > bmid)
      { atmp = frame->abeg;
        btmp = frame->bbeg;
        frame->abeg = amid;
        frame->bbeg = bmid;
#ifdef DEBUG_FIND
        printf("%*sSE:",2*quad->depth+2,"");
#endif
        QuadNode_Find(quad->quads[2],frame,query);
        frame->abeg = atmp;
        frame->bbeg = btmp;
      }

    if (query->aend > amid && query->bbeg < bmid)
      { atmp = frame->abeg;
        btmp = frame->bend;
        frame->abeg = amid;
        frame->bend = bmid;
#ifdef DEBUG_FIND
        printf("%*sSW:",2*quad->depth+2,"");
#endif
        QuadNode_Find(quad->quads[3],frame,query);
        frame->abeg = atmp;
        frame->bend = btmp;
      }
  }
}

QuadLeaf *Plot_Layer(DotPlot *plot, int ilay, Frame *query)
{ Double_Box frame;
  Double_Box qbox;

  SEGS    = plot->layers[ilay]->segs;
  BLOCKS  = NULL;
  FREECNT = BLK_SIZE;

  LIST = (QuadLeaf *) New_Quad();
  LIST->length = 0;

  qbox.abeg = query->x;
  qbox.bbeg = query->y;
  qbox.aend = query->x + query->w;
  qbox.bend = query->y + query->h;

  frame.abeg = 0.;
  frame.bbeg = 0.;
  frame.aend = plot->alen;
  frame.bend = plot->blen;

#ifdef DEBUG_FIND
  printf("..:");
#endif
  QuadNode_Find(plot->layers[ilay]->qtree,&frame,&qbox);

  if (LIST->length > 0)
    { LIST = (QuadLeaf *) New_Quad();
      LIST->length = 0;
    }

  return ((QuadLeaf *) BLOCKS);
}


/*******************************************************************************************
*
*   CREATE MODEL
*
*******************************************************************************************/

static int NSORT(const void *l, const void *r)
{ int x = *((int *) l);
  int y = *((int *) r);

  return (y-x);
}

DotPlot *copyPlot(DotPlot *model)
{ DotPlot *plot;
  int      j;

  plot = malloc(sizeof(DotPlot));
  if (plot == NULL)
    { sprintf(EPLACE,"Cannot allocate plot record\n");
      return (NULL);
    }

  *plot = *model;
  plot->db1->nref += 1;
  plot->db2->nref += 1;
  for (j = 0; j < plot->nlays; j++)
    plot->layers[j]->nref += 1;
  return (plot);
}

static char *simplify_path(char *path)
{ char *r, *s, *n;

  r = path;
  for (s = path; *s != '\0'; s++)
    if (!isspace(*s))
      *r++ = *s;
  *r = '\0';

  r = path+1;
  s = path;
  while ((n = index(s+1,'/')) != NULL)
    { if (s[1] == '.' && s[2] == '.')
        { r -= 2;
          while (*r != '/')
            r -= 1;
          r += 1;
          s = n;
        }
      else if (s+1 == n || (s+2 == n && s[1] == '.'))
        s = n;
      else
        { while (s < n)
            *r++ = *++s;
        }
    }
  while (*s != '\0')
    *r++ = *++s; 

  return (path);
}

static int compare_GDB(GDB *old, GDB *new)
{ int64 sum;
  int   s, c, b, e;

  if (old->ncontig != new->ncontig)
    return (1);
  if (old->nscaff != new->nscaff)
    return (1);
  
  sum = 0;
  for (s = 0; s < old->nscaff; s++)
    { e = old->scaffolds[s].ectg;
      b = old->scaffolds[s].fctg;
      if (new->scaffolds[s].ectg != e)
        return (1);
      if (new->scaffolds[s].fctg != b)
        return (1);
      if (old->scaffolds[s].slen != new->scaffolds[s].slen)
        return (1);
      for (c = b; c < e; c++)
        { if (old->contigs[c].sbeg != new->contigs[c].sbeg + sum)
            return (1);
          if (old->contigs[c].clen != new->contigs[c].clen)
            return (1);
        }
      sum += old->scaffolds[s].slen;
    }

  return (0);
}

DotPlot *createPlot(char *alnPath, int lCut, int iCut, int sCut, DotPlot *model)
{ DotPlot    *plot;
  OneFile    *input;
  char       *src1_name, *src2_name, *cpath;
  DotGDB     *db1, *db2;
  GDB        *gdb1, *gdb2;
    Hash_Table   *hash1, *hash2;
    GDB_CONTIG   *contigs1, *contigs2;
    GDB_SCAFFOLD *scaffs1, *scaffs2;
    int           nscaff1, nscaff2;
  int   tspace;
  int64 novl;

  if (model == NULL)
    { plot = malloc(sizeof(DotPlot));
      if (plot == NULL)
        { sprintf(EPLACE,"Cannot allocate plot record\n");
          return (NULL);
        }
    }
  else
    plot = model;

  //  Initiate .1aln file reading and read header information

  { char  *pwd, *root;
    FILE  *test;

    pwd   = PathTo(alnPath);
    root  = Root(alnPath,".1aln");
    input = open_Aln_Read(Catenate(pwd,"/",root,".1aln"),1,
                           &novl,&tspace,&src1_name,&src2_name,&cpath);
    if (input == NULL)
      { sprintf(EPLACE,"Could not open .1aln file %s/%s.1aln\n",pwd,root);
        free(root);
        free(pwd);
        if (model == NULL)
          free(plot);
        return (NULL);
      }
    free(root);
    free(pwd);

    test = fopen(src1_name,"r");
    if (test == NULL)
      { if (*src1_name != '/')
          test = fopen(Catenate(cpath,"/",src1_name,""),"r");
        if (test == NULL)
          { sprintf(EPLACE,"Could not find GDB %s\n",src1_name);
            goto error1;
          }
        pwd = Strdup(Catenate(cpath,"/",src1_name,""),"Allocating expanded name");
        if (pwd == NULL)
          { fclose(test);
            goto error1;
          }
        free(src1_name);
        src1_name = pwd;
      }
    else
      { if (*src1_name != '/')
          { free(cpath);
            cpath = getcwd(NULL,0);
            pwd = Strdup(Catenate(cpath,"/",src1_name,""),"Allocating expanded name");
            if (pwd == NULL)
              { fclose(test);
                goto error1;
              }
            free(src1_name);
            src1_name = pwd;
          }
      }
    fclose(test);

    if (src2_name != NULL)
      { test = fopen(src2_name,"r");
        if (test == NULL)
          { if (*src2_name != '/')
              test = fopen(Catenate(cpath,"/",src2_name,""),"r");
            if (test == NULL)
              { sprintf(EPLACE,"Could not find GDB %s\n",src2_name);
                goto error1;
              }
            pwd = Strdup(Catenate(cpath,"/",src2_name,""),"Allocating expanded name");
            if (pwd == NULL)
              { fclose(test);
                goto error1;
              }
            free(src2_name);
            src2_name = pwd;
          }
        else
          { if (*src2_name != '/')
              { free(cpath);
                cpath = getcwd(NULL,0);
                pwd = Strdup(Catenate(cpath,"/",src2_name,""),"Allocating expanded name");
                if (pwd == NULL)
                  { fclose(test);
                    goto error1;
                  }
                free(src2_name);
                src2_name = pwd;
              }
          }
        fclose(test);
      }
    free(cpath);
    cpath = NULL;
  }

  //  Create possibly temporary (if model != NULL) GDB's / DB records

  { char  *spath, *tpath;
    int    type;
    GDB   _gdb1, _gdb2;
    int    nogood;

    simplify_path(src1_name);
    if (src2_name != NULL)
      simplify_path(src2_name);

    gdb1 = &_gdb1;
    gdb2 = &_gdb2;

    type = Get_GDB_Paths(src1_name,NULL,&spath,&tpath,0);
    if (type < 0)
      goto error1;
    if (type != IS_GDB)
      nogood = (Create_GDB(gdb1,spath,type,1,tpath) == NULL);
    else
      nogood = Read_GDB(gdb1,tpath);
    free(spath);
    free(tpath);
    if (nogood)
      goto error1;
  
    if (src2_name != NULL)
      { type = Get_GDB_Paths(src2_name,NULL,&spath,&tpath,0);
        if (type < 0)
          { Close_GDB(gdb1);
            goto error1;
          }
        if (type != IS_GDB)
          nogood = (Create_GDB(gdb2,spath,type,1,tpath) == NULL);
        else
          nogood = Read_GDB(gdb2,tpath);
        free(spath);
        free(tpath);
        if (nogood)
          { Close_GDB(gdb1);
            goto error1;
          }
      }
    else
      gdb2 = gdb1;

    //  If first layer, make DB records and transfer GDB's, otherwise check the GDB's
    //     are equal and then close.

    if (model == NULL)
      { db1 = malloc(sizeof(DotGDB));
        db2 = NULL;
        if (db1 == NULL)
          { sprintf(EPLACE,"Cannot allocate GDB record\n");
            goto error2;
          }
        db1->nref = 1;
        db1->name = src1_name;
        db1->gdb  = _gdb1;
        if (src2_name == NULL)
          { db2 = db1;
            db2->nref += 1;
          }
        else
          { db2 = malloc(sizeof(DotGDB));
            if (db2 == NULL)
              { sprintf(EPLACE,"Cannot allocate GDB record\n");
                goto error2;
              }
            db2->nref = 1;
            db2->name = src2_name;
            db2->gdb  = _gdb2;
          }

        plot->db1 = db1;
        plot->db2 = db2;
        plot->dotmemory = dotplot_memory();
      }
    else
      { int comp1, comp2;

        db1 = model->db1;
        db2 = model->db2;

        free(src1_name);
        free(src2_name);
        src1_name = NULL;
        src2_name = NULL;

        comp1 = compare_GDB(&(db1->gdb),gdb1);
        comp2 = 0;
        if (db1 != db2)
          comp2 = compare_GDB(&(db2->gdb),gdb2);

        if (gdb2 != gdb1)
          Close_GDB(gdb2);
        Close_GDB(gdb1);

        if (comp1)
          { sprintf(EPLACE,"1st Genome is not the same as the base layer\n");
            goto error1;
          }
        if (comp2)
          { sprintf(EPLACE,"2nd Genome is not the same as the base layer\n");
            goto error1;
          }
      }

    gdb1 = &(db1->gdb);
    gdb2 = &(db2->gdb);
  }

  nscaff1  = gdb1->nscaff;
  nscaff2  = gdb2->nscaff;
  scaffs1  = gdb1->scaffolds;
  scaffs2  = gdb2->scaffolds;
  contigs1 = gdb1->contigs;
  contigs2 = gdb2->contigs;

  //  Set up scaffold name dictionary

  if (model == NULL)
    { int   s;
      char *head, *sptr, *eptr;
  
      hash1 = New_Hash_Table(nscaff1,0);
      if (hash1 == NULL)
        goto error2;
      head  = gdb1->headers;
      for (s = 0; s < nscaff1; s++)
        { sptr = head + scaffs1[s].hoff;
          for (eptr = sptr; *eptr != '\0'; eptr++)
            if (isspace(*eptr))
              break;
          *eptr = '\0';
          if (Hash_Lookup(hash1,sptr) < 0)
            { if (Hash_Add(hash1,sptr) < 0)
                goto error3;
            }
          else
            { sprintf(EPLACE,"Duplicate scaffold name: %s\n",sptr);
              goto error3;
            }
        }
  
      if (db1 != db2)
        { hash2 = New_Hash_Table(nscaff2,0);
          if (hash1 == NULL)
            goto error3;
          head  = gdb2->headers;
          for (s = 0; s < nscaff2; s++)
            { sptr = head + scaffs2[s].hoff;
              for (eptr = sptr; *eptr != '\0'; eptr++)
                if (isspace(*eptr))
                  break;
              *eptr = '\0';
              if (Hash_Lookup(hash2,sptr) < 0)
                { if (Hash_Add(hash2,sptr) < 0)
                    goto error4;
                }
              else
                { sprintf(EPLACE,"Duplicate scaffold name: %s\n",sptr);
                  goto error4;
                }
            }
        }
      else
        hash2 = hash1;
  
      db1->hash = hash1;
      db2->hash = hash2;
    }

  //  Adjust to global coords

  if (model == NULL)
    { int64 sum;
      int   s, c, b, e;
  
      sum = 0;
      for (s = 0; s < nscaff1; s++)
        { e = scaffs1[s].ectg;
          b = scaffs1[s].fctg;
          for (c = b; c < e; c++)
            contigs1[c].sbeg += sum;
          sum += scaffs1[s].slen;
        }
  
      if (db1 != db2)
        { sum = 0;
          for (s = 0; s < nscaff2; s++)
            { e = scaffs2[s].ectg;
              b = scaffs2[s].fctg;
              for (c = b; c < e; c++)
                contigs2[c].sbeg += sum;
              sum += scaffs2[s].slen;
            }
        }
  
      plot->alen = contigs1[scaffs1[nscaff1-1].fctg].sbeg + scaffs1[nscaff1-1].slen;
      plot->blen = contigs2[scaffs2[nscaff2-1].fctg].sbeg + scaffs2[nscaff2-1].slen;
    }

  //  Add layer

  { Overlap    _ovl, *ovl = &_ovl;
    DotSegment *segs;
    DotLayer   *layer;
    int64       aoff, boff;
    double      iid;
    int         j, k, nlay;

    if (model == NULL)
      { nlay = 1;
        plot->layers[0] = NULL;
      }
    else
      nlay = plot->nlays;
    if (nlay >= MAX_LAYERS)
      { sprintf(EPLACE,"Cannot have more than %d layers\n",MAX_LAYERS);
        goto error4;
      }
 
    segs = malloc(sizeof(DotSegment)*novl);
    if (segs == NULL)
      { sprintf(EPLACE,"Cannot allocate memory for %lld alignments\n",novl);
        goto error4;
      }

#ifdef DEBUG_LAYER
    printf("Initial ovls = %lld\n",novl); fflush(stdout);
#endif

    k = 0;
    for (j = 0; j < novl; j++)
      { Read_Aln_Overlap(input,ovl);
        Skip_Aln_Trace(input);

        if (ovl->path.aepos - ovl->path.abpos <= sCut)
          continue;

        iid  = 100. - (100. * ovl->path.diffs) / (ovl->path.aepos - ovl->path.abpos);

        if (iid <= iCut)
          continue;

        aoff = contigs1[ovl->aread].sbeg;
        segs[k].abeg = ovl->path.abpos + aoff;
        segs[k].aend = ovl->path.aepos + aoff;

        if (COMP(ovl->flags))
          { boff = contigs2[ovl->bread].sbeg + contigs2[ovl->bread].clen;
            segs[k].bbeg = boff - ovl->path.bbpos;
            segs[k].bend = boff - ovl->path.bepos;
          }
        else
          { boff = contigs2[ovl->bread].sbeg;
            segs[k].bbeg = ovl->path.bbpos + boff;
            segs[k].bend = ovl->path.bepos + boff;
          }

        segs[k].iid  = (int) iid;
        segs[k].idx  = j;
        segs[k].mark = 0;

        k += 1;
      }

    if (k < novl)
      { novl = k;
        segs = realloc(segs,sizeof(DotSegment)*novl);
      }

#ifdef DEBUG_LAYER
    printf("1st cull to = %lld\n",novl); fflush(stdout);
#endif

    if (lCut >= 0 && novl > lCut)
      { int *sarray;
        int  alen, digits;

        sarray = malloc(sizeof(int)*novl);
        if (sarray == NULL)
          { free(segs);
            goto error4;
          }

        for (j = 0; j < novl; j++)
          sarray[j] = segs[j].aend - segs[j].abeg;

        qsort(sarray, novl, sizeof(int), NSORT);

        alen = sarray[lCut-1];   // replace low digits with 0 up to loosing 10%

#ifdef DEBUG_LAYER
        printf("%d'th length = %d\n",lCut,alen); fflush(stdout);
#endif

        digits = 1;
        while (1)
          { if ((alen/digits)*digits < .9*alen)
              break;
            digits *= 10;
          }
        digits /= 10;
        alen    = (alen/digits)*digits;

#ifdef DEBUG_LAYER
        printf("Adjusted length = %d\n",alen); fflush(stdout);
#endif

        k = 0;
        for (j = 0; j < novl; j++)
          if (segs[j].aend - segs[j].abeg >= alen)
            segs[k++] = segs[j];

        free(sarray);

        novl = k;
        segs = realloc(segs,sizeof(DotSegment)*novl);

#ifdef DEBUG_LAYER
        printf("Culled to = %lld\n",novl); fflush(stdout);
#endif
      }

#ifdef DEBUG_LAYER
    printf("Start quad tree\n"); fflush(stdout);
#endif

    layer = malloc(sizeof(DotLayer));
    if (layer == NULL)
      { sprintf(EPLACE,"Could not allocate layer record\n");
        free(segs);
        goto error4;
      }

    layer->nref   = 1;
    layer->name   = Root(alnPath,NULL);
    layer->input  = input;
    layer->tspace = tspace;
    layer->novls  = novl;
    layer->segs   = segs;

    plot->layers[nlay] = layer;
    plot->nlays = nlay+1;

    Make_QuadTree(plot,nlay);

    (void) Show_QuadTree;
    (void) Stat_QuadTree;

    // Show_QuadTree(plot,nlay);
    // Stat_QuadTree(plot,nlay);
  }

  return (plot);

error4:
  if (model == NULL)
    Free_Hash_Table(hash2);
error3:
  if (model == NULL)
    Free_Hash_Table(hash1);
error2:
  if (model == NULL)
    { if (db1 != db2)
        free(db2);
      free(db1);
      if (gdb2 != gdb1)
        Close_GDB(gdb2);
      Close_GDB(gdb1);
    }
error1:
  free(cpath);
  free(src2_name);
  free(src1_name);
  oneFileClose(input);
  if (model == NULL)
    free(plot);
  return (NULL);
}

void Free_List(QuadLeaf *list)
{ QuadNode *block, *nlock;

  block = (QuadNode *) list;
  while (block != NULL)
    { nlock = block[BLK_SIZE].quads[0];
      free(block);
      block = nlock;
    }
}

static void Free_DotGDB(DotGDB *db)
{ if (db->nref-- > 1)
    return;
  Free_Hash_Table(db->hash);
  free(db->name);
  Close_GDB(&(db->gdb));
  free(db);
}

void Free_DotPlot(DotPlot *plot)
{ int i;

  for (i = 0; i < plot->nlays; i++)
    { if (plot->layers[i] == NULL || plot->layers[i]->nref-- > 1)
        continue;
      Free_List((QuadLeaf *) (plot->layers[i]->blocks));
      free(plot->layers[i]->name);
      free(plot->layers[i]->segs);
      oneFileClose(plot->layers[i]->input);
    }
  free(plot->dotmemory);
  Free_DotGDB(plot->db1);
  Free_DotGDB(plot->db2);
  free(plot);
}


/*******************************************************************************************
*
*   ALIGNMENT GENERATION
*
*******************************************************************************************/

static void print_seq(char *seq)
{ static char dna[4] = { 'a', 'c', 'g', 't' };
  int i;

  for (i = 0; seq[i] != 4; i++)
    printf("%c",dna[(int) seq[i]]);
  printf("\n");
  fflush(stdout);
}

static int   textlen;
static char *alnptr;

static void align_length(char *bit)
{ textlen += strlen(bit); }

static void align_cat(char *bit)
{ alnptr += sprintf(alnptr,"%s",bit); }

char *create_alignment(DotPlot *plot, DotLayer *layer, DotSegment *seg, char **title)
{ static Work_Data *work = NULL;
  static int        a_max = 0;
  static char      *aseq = NULL;
  static int        b_max = 0;
  static char      *bseq = NULL;
  static int        t_max = 0;
  static uint16    *trace = NULL;
  static int        l_max = 0;
  static char      *atext = NULL;
  static char      *tptr;
  static char       atitle[500];

  Overlap   _ovl, *ovl = &_ovl;
  Alignment _aln, *aln = &_aln;
  Path      *path = &(_ovl.path);

  GDB          *gdb1    = &(plot->db1->gdb);
  GDB          *gdb2    = &(plot->db2->gdb);
  OneFile      *in      = layer->input;
  GDB_SCAFFOLD *ascaffs = gdb1->scaffolds;
  GDB_SCAFFOLD *bscaffs = gdb2->scaffolds;
  GDB_CONTIG   *acontig = gdb1->contigs;
  GDB_CONTIG   *bcontig = gdb2->contigs;

  int acont, bcont;
  int ascaf, bscaf;
  int64 aclen, bclen;
  int64 aslen, bslen;
  int64 aoffs, boffs;
  int   amin, amax;
  int   bmin, bmax;

  (void) print_seq;

  if ( ! oneGoto(in,'A',seg->idx+1))
    return (NULL);

  if (work == NULL)
    { work = New_Work_Data();
      if (work == NULL)
        return (NULL);
    }

  if (in->info['T']->given.max > t_max)
    { t_max = in->info['T']->given.max;
      trace = (uint16 *) Realloc(trace,2*sizeof(uint16)*t_max,"Allocating trace vector");
      if (trace == NULL)
        { t_max = 0;
          return (NULL);
        }
    }

  oneReadLine(in);
  Read_Aln_Overlap(in,ovl);
  path->tlen  = Read_Aln_Trace(in,(uint8 *) trace);
  path->trace = trace;

  acont = ovl->aread;
  bcont = ovl->bread;
  ascaf = acontig[acont].scaf;
  bscaf = bcontig[bcont].scaf;
  aclen = acontig[acont].clen;
  bclen = bcontig[bcont].clen;
  aslen = ascaffs[ascaf].slen;
  bslen = bscaffs[bscaf].slen;
  aoffs = acontig[acont].sbeg - acontig[ascaffs[ascaf].fctg].sbeg;
  boffs = bcontig[bcont].sbeg - bcontig[bscaffs[bscaf].fctg].sbeg;

  aln->path  = path = &(ovl->path);
  aln->alen  = aclen;
  aln->blen  = bclen;
  aln->flags = ovl->flags;

  { int64 ab, ae;
    int64 bb, be;

    tptr = atitle + Number_To_String((int64) ascaf+1,0,atitle);
    tptr += sprintf(tptr,".%dn",(acont - ascaffs[ascaf].fctg)+1);
    ab = aoffs + path->abpos;
    ae = aoffs + path->aepos;
    if (ab == 0 || ab == aslen)
      tptr += sprintf(tptr,"<");
    else
      tptr += sprintf(tptr,"[");
    tptr += Number_To_String((int64) ab,0,tptr);
    tptr += sprintf(tptr,"..");
    tptr += Number_To_String((int64) ae,0,tptr);
    if (ae == 0 || ae == aslen)
      tptr += sprintf(tptr,"> x ");
    else
      tptr += sprintf(tptr,"] x ");

    tptr += Number_To_String((int64) bscaf+1,0,tptr);
    tptr += sprintf(tptr,".%d%c",(bcont - bscaffs[bscaf].fctg)+1,(COMP(ovl->flags) == 0)?'n':'c');
    if (COMP(ovl->flags))
      { bb = boffs+(bclen-path->bbpos);
        be = boffs+(bclen-path->bepos);
      }
    else
      { bb = boffs+path->bbpos;
        be = boffs+path->bepos;
      }
    if (bb == 0 || bb == bslen)
      tptr += sprintf(tptr,"<");
    else
      tptr += sprintf(tptr,"[");
    tptr += Number_To_String((int64) bb,0,tptr);
    tptr += sprintf(tptr,"..");
    tptr += Number_To_String((int64) be,0,tptr);
    if (be == 0 || be == bslen)
      tptr += sprintf(tptr,">");
    else
      tptr += sprintf(tptr,"]");
  }

  amin = path->abpos;
  amax = path->aepos;
  if (COMP(aln->flags))
    { bmin = bclen - path->bepos;
      bmax = bclen - path->bbpos;
    }
  else
    { bmin = path->bbpos;
      bmax = path->bepos;
    }

  if (amax-amin > a_max)
    { a_max = amax-amin;
      if (aseq == NULL)
        aseq = (char *) Malloc(a_max+12,"Allocating sequence buffer");
      else
        aseq = (char *) Realloc(aseq-1,a_max+12,"Allocating sequence buffer");
      if (aseq == NULL)
        { a_max = 0;
          return (NULL);
        }
      aseq += 1;
    }

  if (bmax-bmin > b_max)
    { b_max = bmax-bmin;
      if (bseq == NULL)
        bseq = (char *) Malloc(b_max+12,"Allocating sequence buffer");
      else
        bseq = (char *) Realloc(bseq-1,b_max+12,"Allocating sequence buffer");
      if (bseq == NULL)
        { b_max = 0;
          return (NULL);
        }
      bseq += 1;
    }

  Decompress_TraceTo16(ovl);

  aln->aseq = Get_Contig_Piece(gdb1,acont,amin,amax,NUMERIC,aseq);
  aln->bseq = Get_Contig_Piece(gdb2,bcont,bmin,bmax,NUMERIC,bseq);

  aln->aseq -= amin;
  if (COMP(aln->flags))
    { Complement_Seq(aln->bseq,bmax-bmin);
      aln->bseq -= (bclen-bmax);
    }
  else
    aln->bseq -= bmin;

  Compute_Trace_PTS(aln,work,layer->tspace,GREEDIEST);

  Gap_Improver(aln,work);

  { int   tlen   = path->tlen;
    int  *itrace = (int *) path->trace;
    int   i;

    path->abpos += aoffs;
    path->aepos += aoffs;
    aln->alen = 0;
    path->bbpos += boffs;
    path->bepos += boffs;
    if (COMP(ovl->flags))
      aln->blen += 2*boffs;
    else
      aln->blen = 0;

    aln->aseq -= aoffs;
    aln->bseq -= boffs;
    for (i = 0; i < tlen; i++)
      if (itrace[i] < 0)
        itrace[i] -= aoffs;
      else
        itrace[i] += boffs;
  }

  textlen = 0;
  Transmit_Alignment(align_length,aln,work,100,0,0,9,0);

  if (textlen > l_max)
    { l_max = 1.2*textlen + 1001; 
      atext = Realloc(atext,l_max+1,"Reallocating alignment text");
      if (atext == NULL)
        { l_max = 0;
          return (NULL);
        }
    }

  alnptr = atext;
  Transmit_Alignment(align_cat,aln,work,100,0,0,9,0);

  *title = atitle; 
  return (atext);
}
