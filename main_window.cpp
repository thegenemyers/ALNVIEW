#include <stdlib.h>
#include <math.h>

#include <QtGui>

extern "C" {
#include "GDB.h"
#include "float.h"
#include "gene_core.h"
#include "sticks.h"
#include "doter.h"
#include "select.h"
}

#undef DEBUG

#include "main_window.h"
#include "open_window.h"

#define CANVAS_MARGIN 20

#define PICK_LIMIT 1000000

#define LOCATOR_RECTANGLE_SIZE 100

QRect *DotWindow::screenGeometry = NULL;
int    DotWindow::windowWidth;
int    DotWindow::windowHeight;

int    DotCanvas::labelWidth;

/***************************************************************************************/
/*                                                                                     */
/*   CONSTANTS & CUSTOMIZED WIDGETS                                                    */
/*                                                                                     */
/***************************************************************************************/

DotLineEdit::DotLineEdit(QWidget *parent) : QLineEdit(parent) {}

void DotLineEdit::setChain(QWidget *p, QWidget *s)
{ mypred = p;
  mysucc = s;
}

void DotLineEdit::focusOutEvent(QFocusEvent *event)
{ if (event->reason() <= Qt::OtherFocusReason)
    emit focusOut();
  QLineEdit::focusOutEvent(event);
}

void DotLineEdit::keyPressEvent(QKeyEvent *event)
{ int key = event->key();
  if (key == Qt::Key_Return)  // || key == Qt::Key_Tab || key == Qt::Key_Backtab)
    { if ((event->modifiers() & Qt::ShiftModifier) != 0)
        mypred->setFocus();
      else
        mysucc->setFocus();
    }
  else
    QLineEdit::keyPressEvent(event);
}

bool DotLineEdit::focusNextPrevChild(bool next)
{ (void) next;
  return (false);
}

SeqLineEdit::SeqLineEdit(QWidget *parent, int id) : QLineEdit(parent)
{ setFixedHeight(20);
  setStyleSheet("border: 1px solid black");
  setReadOnly(true);
  setFont(QFont(tr("Monaco")));
  myid = id;
}

void SeqLineEdit::enterEvent(QEnterEvent *)
{ setFocus(); }

void SeqLineEdit::leaveEvent(QEvent *)
{ clearFocus(); }

void SeqLineEdit::keyPressEvent(QKeyEvent *event)
{ int key = event->key();
  if (key == Qt::Key_Left)
    emit focusShift(-1,myid);
  else if (key == Qt::Key_Right)
    emit focusShift(+1,myid);
  else
    return;
}

void SeqLineEdit::mousePressEvent(QMouseEvent *event)
{ (void) event; }

void SeqLineEdit::mouseMoveEvent(QMouseEvent *event)
{ (void) event; }

void SeqLineEdit::mouseReleaseEvent(QMouseEvent *event)
{ (void) event; }

MyMenu::MyMenu(QWidget *parent) : QMenu(parent) { }

void MyMenu::mouseReleaseEvent(QMouseEvent *ev)
{ QAction *act = actionAt(ev->pos());

  if (act != NULL)
    act->trigger();
  QCoreApplication::sendEvent(parent(),ev);
}

void MyMenu::mousePressEvent(QMouseEvent *ev)
{ QCoreApplication::sendEvent(parent(),ev); }

void MyMenu::keyReleaseEvent(QKeyEvent *ev)
{ QCoreApplication::sendEvent(parent(),ev); }


/***************************************************************************************/
/*                                                                                     */
/*   DOT CANVAS                                                                        */
/*                                                                                     */
/***************************************************************************************/

DotCanvas::~DotCanvas()
{ }

DotCanvas::DotCanvas(QWidget *parent) : QWidget(parent)
{ setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
  noFrame = true;
  pickedSeg = NULL;
  rubber  = new QRubberBand(QRubberBand::Rectangle, this);
  timer   = new QBasicTimer();

  image  = NULL;
  raster = new uchar *[screen()->availableGeometry().height()];

  setAttribute(Qt::WA_KeyCompression,false);

  aline = new QAction(tr(""),this);
    aline->setFont(QFont(tr("Monaco"),11));
    aline->setEnabled(false);

  bline = new QAction(tr(""),this);
    bline->setFont(QFont(tr("Monaco"),11));
    bline->setEnabled(false);

  alignAct = new QAction(tr("Alignment"),this);
    alignAct->setFont(QFont(tr("Monaco"),11));
    alignAct->setToolTip(tr("Show the alignment for this segment"));

  popup = new MyMenu(this);
    popup->addAction(aline);
    popup->addAction(bline);
    popup->addAction(alignAct);

  QPalette pal = popup->palette();
    pal.setBrush(QPalette::Disabled,QPalette::Text,QBrush(Qt::black));
    popup->setPalette(pal);

  QAction *mesgAct = new QAction(tr("Zoom in to width < 1Mbp for pick"),this);
    mesgAct->setFont(QFont(tr("Monaco"),11));
    mesgAct->setToolTip(tr("Show the alignment for this segment"));
    mesgAct->setEnabled(false);

  faraway = new MyMenu(this);
    faraway->addAction(mesgAct);

  QPalette paf = faraway->palette();
    paf.setBrush(QPalette::Disabled,QPalette::Text,QBrush(Qt::black));
    faraway->setPalette(paf);

  connect(alignAct,SIGNAL(triggered()),this,SLOT(showAlign()));
}

Frame *DotCanvas::shareData(DotState *s, DotPlot *p)
{ state = s;
  plot  = p;
  if (plot->db1->gdb.seqs == NULL || plot->db2->gdb.seqs == NULL)
    popup->actions().at(2)->setEnabled(false);
  return (&frame);
}

#define ASSIGN(frame,fx,fy,fw,fh)	\
{ frame.x = fx;				\
  frame.y = fy;				\
  frame.w = fw;				\
  frame.h = fh;				\
}

bool DotCanvas::zoomView(double zoomDel)
{ double vX = frame.x;
  double vY = frame.y;
  double vW = frame.w;
  double vH = frame.h;
  double nW = vW / zoomDel;
  double nH = vH / zoomDel;

  if (plot == NULL)
    return (false);

  if (nW > plot->alen && nH > plot->blen)
    { if (plot->alen/nW > plot->blen/nH)
        { nH = nH*(plot->alen/nW);
          nW = plot->alen;
        }
      else
        { nW = nW*(plot->blen/nH);
          nH = plot->blen;
        }
    }
  double nX = vX + (vW-nW)/2;
  double nY = vY + (vH-nH)/2;
  if (nW > plot->alen)
    nX = (plot->alen-nW)/2.;
  else
    { if (nX < 0.) nX = 0.;
      if (nX + nW > plot->alen) nX = plot->alen - nW; 
    }
  if (nH > plot->blen)
    nY = (plot->blen-nH)/2.;
  else
    { if (nY < 0.) nY = 0.;
      if (nY + nH > plot->blen) nY = plot->blen - nH; 
    }

  // if (nW < rectW)
    // return (false);

#ifdef DEBUG
  printf("Frame (zoom) = (%.1f,%.1f) %.1f vs %.1f\n",nX,nY,nW,nH);
#endif

  ASSIGN(frame,nX,nY,nW,nH);
  
  double zW = (1.*plot->alen) / nW;
  double zH = (1.*plot->blen) / nH;
  if (zW > zH)
    emit NewFrame(zW);
  else
    emit NewFrame(zH);

  return (true);
}

void DotCanvas::resizeEvent(QResizeEvent *event)
{ QSize oldS = event->oldSize();
  QSize newS = event->size();

  rectW = newS.width();
  rectH = newS.height();

  if (noFrame)
    { event->accept();
      return;
    }

  double vX = frame.x;
  double vY = frame.y;
  double vW = frame.w;
  double vH = frame.h;

  int oW = oldS.width();
  int nW = newS.width();
  int oH = oldS.height();
  int nH = newS.height();

  vW = (vW / oW) * nW;
  vH = (vH / oH) * nH;
  if (frame.w >= plot->alen && frame.h >= plot->blen)
    { if (plot->alen/vW > plot->blen/vH)
        { vH = vH*(plot->alen/vW);
          vW = plot->alen;
          vX = 0;
          vY = (plot->blen-vH)/2.;
        }
      else
        { vW = vW*(plot->blen/vH);
          vH = plot->blen;
          vX = (plot->alen-vW)/2.;
          vY = 0;
        }
    }
  else
    { if (vW > plot->alen)
        vX = (plot->alen-vW)/2.;
      else if (vX < 0)
        vX = 0.;
      else if (vX + vW > plot->alen)
        vX = plot->alen - vW;
      if (vH > plot->blen)
        vY = (plot->blen-vH)/2.;
      else if (vY < 0.)
        vY = 0.;
      else if (vY + vH > plot->blen)
        vY = plot->blen - vH;
    }

  ASSIGN(frame,vX,vY,vW,vH);

  event->accept();

  double zW = (1.*plot->alen) / vW;
  double zH = (1.*plot->blen) / vH;
  if (zW > zH)
    emit NewFrame(zW);
  else
    emit NewFrame(zH);
}

void DotCanvas::resetView()
{ double vX, vY, vW, vH;

  int rectW = width();
  int rectH = height();
  double h = (1.*rectW)/plot->alen;
  double v = (1.*rectH)/plot->blen;

  if (h > v)
    { vY = 0;
      vH = plot->blen;
      vW = rectW / v;
      vX = (plot->alen - vW) / 2.;
    }
  else
    { vX = 0;
      vW = plot->alen;
      vH = rectH / h;
      vY = (plot->blen - vH) / 2.;
    }

  ASSIGN(frame,vX,vY,vW,vH);

#ifdef DEBUG
  printf("Reset frame = (%.1f,%.1f) %.8f vs %.8f\n",vX,vY,vW,vH);
#endif
}

bool DotCanvas::viewToFrame()
{ double vX, vY, vW, vH;

  double h = (1.*rectW)/state->view.w;
  double v = (1.*rectH)/state->view.h;

  if (h > v)
    { vY = state->view.y;
      vH = state->view.h;
      vW = rectW / v;
      if (vW > plot->alen)
        vX = (plot->alen - vW) / 2.;
      else
        { vX = state->view.x + (state->view.w - vW) / 2.;
          if (vX < 0.)
            vX = 0.;
          if (vX + vW > plot->alen)
            vX = plot->alen - vW;
        }
    }
  else
    { vX = state->view.x;
      vW = state->view.w;
      vH = rectH / h;
      if (vH > plot->blen)
        vY = (plot->blen - vH) / 2.;
      else
        { vY = state->view.y + (state->view.h - vH) / 2.;
          if (vY < 0.)
            vY = 0.;
          if (vY + vH > plot->blen)
            vY = plot->blen - vH;
        }
    }

  // if (vW < rectW)
    // return (false);

  ASSIGN(frame,vX,vY,vW,vH);

#ifdef DEBUG
  printf("view 2 frame = (%.1f,%.1f) %.8f vs %.8f\n",vX,vY,vW,vH);
#endif

  double zW = (1.*plot->alen) / vW;
  double zH = (1.*plot->blen) / vH;
  if (zW > zH)
    emit NewFrame(zW);
  else
    emit NewFrame(zH);

  return (true);
}

void DotCanvas::enterEvent(QEnterEvent *)
{ QApplication::setOverrideCursor(Qt::ArrowCursor);
  setFocus();
#ifdef DEBUG
  printf("Enter\n");
#endif
}

void DotCanvas::leaveEvent(QEvent *)
{ QApplication::restoreOverrideCursor();
  clearFocus();
#ifdef DEBUG
  printf("Exit\n");
#endif
}

void DotCanvas::keyPressEvent(QKeyEvent *event)
{ if ((event->modifiers() & Qt::ShiftModifier) != 0)
    QApplication::changeOverrideCursor(Qt::CrossCursor);
  else
    QApplication::changeOverrideCursor(Qt::ArrowCursor);

  int key = event->key();
  if (key == lastKey)
    return;

  lastKey = key;
  if (key == Qt::Key_Z)
    zoomPop();
  else if (key == Qt::Key_Plus)
    { zoomView(1.41421356237);
      update();
    }
  else if (key == Qt::Key_Minus)
    { zoomView(0.70710678118);
      update();
    }
}

void DotCanvas::zoomPop()
{ double vX, vY, vW, vH;
  double mag, xC, yC;

  if (state->zMag.isEmpty())
    return;
  mag  = state->zMag.pop();
  xC   = state->zXct.pop();
  yC   = state->zYct.pop();
  if (mag == plot->alen)
    { resetView();
      emit NewFrame(1.0);
      update();
      return;
    }
  vW = mag*rectW;
  vH = mag*rectH;
  vX = xC - vW/2.;
  vY = yC - vH/2.;

  ASSIGN(frame,vX,vY,vW,vH);

  double zW = (1.*plot->alen) / vW;
  double zH = (1.*plot->blen) / vH;
  if (zW > zH)
    emit NewFrame(zW);
  else
    emit NewFrame(zH);

  update();
}

void DotCanvas::keyReleaseEvent(QKeyEvent *event)
{ if ((event->modifiers() & Qt::ShiftModifier) != 0)
    QApplication::changeOverrideCursor(Qt::CrossCursor);
  else
    QApplication::changeOverrideCursor(Qt::ArrowCursor);

  lastKey = -1;
#ifdef DEBUG
  printf("Key Release\n"); fflush(stdout);
#endif
}

DotSegment *DotCanvas::pick(int ex, int ey, DotLayer **pickedLayer)
{ QuadLeaf   *list;
  int         j, k;
  Frame       pframe;
  double      x, y;
  double      close;
  int         bestj, besti;

#ifdef DEBUG
  printf("Pick click\n");
#endif

  x = frame.x + ((ex-20.)/(rectW-40.))*frame.w;
  y = frame.y + ((ey-20.)/(rectH-40.))*frame.h;
  if (x < 0. || x < frame.x)
    return (NULL);
  if (x > plot->alen || x > frame.x+frame.w)
    return (NULL);
  if (y < 0. || y < frame.y)
    return (NULL);
  if (y > plot->blen || y > frame.y+frame.h)
    return (NULL);

  pframe.x = x-5000;
  pframe.y = y-5000;
  pframe.w = pframe.h = 10000;

  close = PICK_LIMIT+1;
  bestj = -1;
  besti = -1;
  for (j = state->nlays-1; j >= 0; j--)
    { k = state->order[j];
      if ( ! state->on[k] || k == 0)
        continue;

      list = Plot_Layer(plot,k,&pframe);

      { QuadLeaf   *curn, *next;
        DotSegment *segs, *line;
        int64       xbeg, xend;
        int64       ybeg, yend;
        int         i;
        double      d;

        segs = plot->layers[k]->segs;
        for (curn = list; curn != NULL; curn = next)
          { next = (QuadLeaf *) (((QuadNode *) curn)[BLK_SIZE].quads[0]);
            while (curn->length > 0)
              { for (i = 0; i < curn->length; i++)
                  { line = segs + curn->idx[i];
                    xbeg = line->abeg;
                    xend = line->aend;
                    ybeg = line->bbeg;
                    yend = line->bend;
                    if (abs(yend - ybeg) > abs(xend-xbeg))
                      { if (yend > ybeg)
                          { if (y >= yend)
                              d = (y-yend) + fabs(x-xend);
                            else if (y <= ybeg)
                              d = (ybeg-y) + fabs(x-xbeg);
                            else
                              d = fabs(x - (xbeg+(xend-xbeg)*((y-ybeg)/(yend-ybeg))));
                          }
                        else
                          { if (y >= ybeg)
                              d = (y-ybeg) + fabs(x-xbeg);
                            else if (y <= yend)
                              d = (yend-y) + fabs(x-xend);
                            else
                              d = fabs(x - (xend+(xbeg-xend)*((y-yend)/(ybeg-yend))));
                          }
                      }
                    else
                      { if (xend > xbeg)
                          { if (x >= xend)
                              d = (x-xend) + fabs(y-yend);
                            else if (x <= xbeg)
                              d = (xbeg-x) + fabs(y-ybeg);
                            else
                              d = fabs(y - (ybeg + (yend-ybeg)*((x-xbeg)/(xend-xbeg))));
                          }
                        else
                          { if (x >= xbeg)
                              d = (x-xbeg) + fabs(y-ybeg);
                            else if (x <= xend)
                              d = (xend-x) + fabs(y-yend);
                            else
                              d = fabs(y - (yend+(ybeg-yend)*((x-xend)/(xbeg-xend))));
                          }
                      }
                    if (d < close)
                      { bestj = k;
                        besti = curn->idx[i];
                        close = d;
                      }
                    line->mark = 0;
                  }
                curn += 1;
              }
          }
      }

      Free_List(list);
    }

  if (besti >= 0 && close < (50.*frame.w)/(rectW-40.))
    { *pickedLayer = plot->layers[bestj];
      return (plot->layers[bestj]->segs + besti);
    }
  else
    { *pickedLayer = NULL;
      return (NULL);
    }
}

void DotCanvas::mousePressEvent(QMouseEvent *event)
{ mouseX = event->position().toPoint().x();
  mouseY = event->position().toPoint().y();
#ifdef DEBUG
  printf("Press %d %d\n",mouseX,mouseY);
#endif
  if ((QApplication::keyboardModifiers() & Qt::ShiftModifier) != 0)
    if (lastKey == Qt::Key_X)
      { menuLock = true;
        goto pick_code;
      }
    else
      { select = true;
        rubber->setGeometry(QRect(mouseX,mouseY,0,0));
        rubber->show();
        QApplication::changeOverrideCursor(Qt::CrossCursor);
      }
  else if (lastKey == Qt::Key_X)
    { menuLock = false;
      goto pick_code;
    }
  else
    { select  = false;
      picking = false;
      nograb  = true;
      scaleX  = frame.w / width();
      scaleY  = frame.h / height();
    }
  return;

pick_code:
  select  = false;
  picking = true;
  pickedSeg = pick(event->position().toPoint().x(),
                   event->position().toPoint().y(),&pickedLayer);
  if (pickedSeg != NULL)
    { char  *s1, *s2;
      double d1, d2;
      int    len, prec;

      s1 = Map_Coord(&(plot->db1->gdb),pickedSeg->abeg,-1,
                     state->format,state->view.w);
      s2 = Map_Coord(&(plot->db2->gdb),-1,pickedSeg->bbeg,
                     state->format,state->view.h);
      aline->setText(tr("Beg: %1,%2").arg(s1).arg(s2));

      d1 = pickedSeg->aend - pickedSeg->abeg;
      d2 = llabs(pickedSeg->bend - pickedSeg->bbeg);
      len = (d1 + d2) / 2.;

      d1 = digits(len,&s1,&prec);

      bline->setText(tr("Len: %1%2   Id: %3\%").
                        arg(len/d1,4,'f',prec).arg(s1).arg(pickedSeg->iid));

      popup->popup(event->globalPosition().toPoint());

      return;
    }
  menuLock = false;
}

void DotCanvas::mouseMoveEvent(QMouseEvent *event)
{ xpos = event->position().toPoint().x();
  ypos = event->position().toPoint().y();
  if (select)
    rubber->setGeometry(QRect(mouseX,mouseY,xpos-mouseX,ypos-mouseY).normalized());
  else if (!picking)
    { if (nograb)
        { if (abs(xpos-mouseX) <= 1 && abs(ypos-mouseY) <= 1)
            return;
          QApplication::changeOverrideCursor(Qt::ClosedHandCursor);
          nograb = false;
        }
      if (frame.w < plot->alen)
        { frame.x += scaleX * (mouseX - xpos);
          if (frame.x < 0.)
            frame.x = 0.;
          if (frame.x + frame.w > plot->alen)
            frame.x = plot->alen - frame.w;
        }
      if (frame.h < plot->blen)
        { frame.y += scaleY * (mouseY - ypos);
          if (frame.y < 0.)
            frame.y = 0.;
          if (frame.y + frame.h > plot->blen)
            frame.y = plot->blen - frame.h;
        }
      mouseX = xpos;
      mouseY = ypos;

      update();
      emit NewFrame(0.);
    }
}

void DotCanvas::mouseReleaseEvent(QMouseEvent *event)
{ if (select)
    { rubber->hide();
#ifdef DEBUG
      printf("Selected\n");
#endif
      QRect reg  = rubber->geometry();
      View  undo = state->view;
      int64 xb = frame.x + ((reg.x()-20.)/(rectW-40.))*frame.w;
      int64 yb = frame.y + ((reg.y()-20.)/(rectH-40.))*frame.h;
      int64 xe = frame.x + (((reg.x()-20.)+reg.width())/(rectW-40.))*frame.w;
      int64 ye = frame.y + (((reg.y()-20.)+reg.height())/(rectH-40.))*frame.h;
      if (xb < 0) xb = 0;
      if (xb < frame.x) xb = frame.x;
      if (yb < 0) yb = 0;
      if (yb < frame.y) yb = frame.y;
      if (xe > plot->alen) xe = plot->alen;
      if (xe > frame.x+frame.w) xe = frame.x+frame.w;
      if (ye > plot->blen) ye = plot->blen;
      if (ye > frame.y+frame.h) ye = frame.y+frame.h;
      state->view.x = xb;
      state->view.y = yb;
      state->view.w = xe-xb;
      state->view.h = ye-yb;
#ifdef DEBUG
      printf("Region select (%d,%d) %d x %d\n",reg.x(),reg.y(),reg.width(),reg.height());
      printf("      = view  (%lld,%lld) (%lld,%lld)\n",xb,yb,xe,ye);
#endif
      if (viewToFrame())
        update();
      else
        state->view = undo;
    }
  else if (picking)
    { if (!menuLock)
        { popup->hide();
          faraway->hide();
          setFocus();
          pickedSeg = NULL;
        }
    }
  else if (nograb)
    { double ex = event->position().toPoint().x()-20.;
      double ey = event->position().toPoint().y()-20.;
#ifdef DEBUG
      printf("Focus release\n");
#endif
      int64 x = frame.x + (ex/(rectW-40))*frame.w;
      int64 y = frame.y + (ey/(rectH-40))*frame.h;
      if (x < 0) x = 0;
      if (x < frame.x) x = frame.x;
      if (x > plot->alen) x = plot->alen;
      if (x > frame.x+frame.w) x = frame.x+frame.w;
      if (y < 0) y = 0;
      if (y < frame.y) y = frame.y;
      if (y > plot->blen) y = plot->blen;
      if (y > frame.y+frame.h) y = frame.y+frame.h;
      state->focus.x = x;
      state->focus.y = y;
      update();
      emit NewFocus();
    }

  QApplication::changeOverrideCursor(Qt::ArrowCursor);
}

void DotCanvas::showAlign()
{ DotWindow   *par;
  AlignWindow *awin;
  char        *title, *align;

  align = create_alignment(plot,pickedLayer,pickedSeg,&title);

  par = (DotWindow *) parent();

  awin = new AlignWindow(tr(title),align,par,Qt::Tool);
  awin->show();
  awin->move(QPoint(state->wGeom.x()+DotWindow::windowHeight,state->wGeom.y()));
  awin->raise();

  // par->addTextBox(awin);
}
 
QVector<QRgb>  DotCanvas::ctable(2);
QVector<uchar> DotCanvas::imbit(8);

void DotCanvas::paintEvent(QPaintEvent *event)
{ QPainter  painter;
  double    cMag;

  (void) event;

  painter.begin(this);
  painter.setRenderHint(QPainter::Antialiasing,true);
  painter.setRenderHint(QPainter::SmoothPixmapTransform,true);

  if (noFrame)    //  First paint => set Frame (cannot do earlier as do not know size)
    { viewToFrame();
      noFrame = false;
      QFontMetrics met = QFontMetrics(painter.font());
      QRect br = met.boundingRect(tr("8.88M"));
      labelWidth = br.width();
    }

  double vX = frame.x;
  double vY = frame.y;
  double vW = frame.w;
  double vH = frame.h;

  double xa = (rectW-44)/vW;
  double xb = -vX*xa+22;

  double ya = (rectH-44)/vH;
  double yb = -vY*ya+22;

  if (vW >= plot->alen && vH >= plot->blen)
    cMag = plot->alen;
  else
    cMag = vW/rectW;
  if (cMag < state->lMag)
    { state->zMag.push(state->lMag);
      state->zXct.push(state->lXct);
      state->zYct.push(state->lYct);
    }
  else if (cMag > state->lMag)
    { while (!state->zMag.isEmpty())
        { double mag;

          mag = state->zMag.pop();
          if (mag > cMag)
            { state->zMag.push(mag);
              break;
            }
          // x = state->zXct.pop();
          // y = state->zYct.pop();
        }
    }
  state->lMag = cMag;
  state->lXct = vX+vW/2.;
  state->lYct = vY+vH/2.;

  painter.fillRect(0,0,rectW,rectH,QColor(0,0,0));

  { QPen dPen;
    int64  b, w;
    int64  c, h;
    int64  x, y, o;
    int64  v;
    int64  unit;
    double den;
    int    prec;
    char  *suf;

    dPen.setColor(QColor(255,255,255));
    dPen.setWidth(2);
    painter.setPen(dPen);
  
    if (xb < 22)
      b = 21;
    else
      b = xb-1;
    w = xb + plot->alen*xa + 1;
    if (w >= rectW-20)
      w = rectW-21;
    w = w-b;

    if (yb < 22)
      c = 21;
    else
      c = yb-1;
    h = yb + plot->blen*ya + 1;
    if (h >= rectH-20)
      h = rectH-21;
    h = h-c;

    painter.drawRect(b,c,w,h);

    den  = (double) digits((int64) state->view.w, &suf, &prec);

    unit = divide_bar((int) ((100./(rectW-40.))*vW));

    if (vX < 0)
      { y = plot->alen;
        o = 0;
      }
    else
      { y = vW;
        o = vX;
      }
    painter.drawLine(b,c,b,c-4);
    painter.drawLine(b,c,b-4,c);
    painter.drawText(QRect(b-21,c-8,20,4),Qt::AlignHCenter|Qt::AlignBottom|Qt::TextDontClip,
                     tr("0%2").arg(suf));
    for (x = unit; x < y; x += unit)
      { v = xb+xa*(x+o);
        painter.drawLine(v,c,v,c-4);
        if (v+40 > b+w)
          { if ((b+w)-v >= labelWidth/2 && unit*xa >= 1.5*labelWidth)
              painter.drawText(QRect(v-50,c-8,50,4),Qt::AlignRight|Qt::AlignBottom|Qt::TextDontClip,
                               tr("%1%2").arg(x/den,0,'f',prec).arg(suf));
          }
        else
          painter.drawText(QRect(v-10,c-8,20,4),Qt::AlignHCenter|Qt::AlignBottom|Qt::TextDontClip,
                           tr("%1%2").arg(x/den,0,'f',prec).arg(suf));
      }
    v = xb+xa*(y+o);
    painter.drawLine(v,c,v,c-4);
    painter.drawText(QRect(v-10,c-8,20,4),Qt::AlignHCenter|Qt::AlignBottom|Qt::TextDontClip,
                     tr("%1%2").arg(y/den,0,'f',prec).arg(suf));

    painter.save();
    painter.rotate(-90.);

    if (vY < 0)
      { y = plot->blen;
        o = 0;
      }
    else
      { y = vH;
        o = vY;
      }
    v = yb;
    painter.drawLine(-v,b,-v,b-4);
    for (x = unit; x < y; x += unit)
      { v = yb+ya*(x+o);
        painter.drawLine(-v,b,-v,b-4);
        if (v+40 > c+h)
          { if ((c+h)-v >= labelWidth/2 && unit*ya >= 1.5*labelWidth)
              painter.drawText(QRect(-v,b-9,50,4),Qt::AlignLeft|Qt::AlignBottom|Qt::TextDontClip,
                         tr("%1%2").arg(x/den,0,'f',prec).arg(suf));
          }
        else
          painter.drawText(QRect(-v-10,b-9,20,4),Qt::AlignHCenter|Qt::AlignBottom|Qt::TextDontClip,
                           tr("%1%2").arg(x/den,0,'f',prec).arg(suf));
      }
    v = yb+ya*(y+o);
    painter.drawLine(-v,b,-v,b-4);
    painter.drawText(QRect(-v-10,b-9,20,4),Qt::AlignHCenter|Qt::AlignBottom|Qt::TextDontClip,
                     tr("%1%2").arg(y/den,0,'f',prec).arg(suf));

    painter.restore();
  }

  int64 cxb, cxe;  // Clipping rectangle
  int64 cyb, cye;  
  
  if (xb < 22)
    cxb = 22;
  else
    cxb = (int) xb;
  cxe = xb + plot->alen*xa;
  if (cxe > rectW-22)
    cxe = rectW-22;
  
  if (yb < 22)
    cyb = 22;
  else
    cyb = (int) yb;
  cye = yb + plot->blen*ya;
  if (cye > rectH-22)
    cye = rectH-22;

  painter.setClipRegion(QRect(cxb,cyb,cxe-cxb,cye-cyb));

  { QPen  iPen;           //  Draw scaffold / contig lines
    int   i, j, x, y;
    int   alf, art;
    int   blf, brt;
    int64 x8, y8;
    int           nsf1,  nsf2;
    GDB_CONTIG   *ctg1, *ctg2;
    GDB_SCAFFOLD *scf1, *scf2;

    QVector<qreal> sdash;
    sdash << 8 << 8;

    QVector<qreal> cdash;
    cdash << 1 << 4;

    ctg1 = plot->db1->gdb.contigs;
    scf1 = plot->db1->gdb.scaffolds;
    nsf1 = plot->db1->gdb.nscaff-1;
    alf = art = 0;
    for (i = 1; i < nsf1; i++)
      { x8 = ctg1[scf1[i].fctg].sbeg;
        x  = (int) (x8*xa+xb);
        if (x >= 22)
          { if (x <= rectW-22)
              { art = i;
                // painter.drawLine(x,cyb,x,cye);
              }
          }
        else
          alf = i;
      }
    if (art < alf)
      art = alf;

    ctg2 = plot->db2->gdb.contigs;
    scf2 = plot->db2->gdb.scaffolds;
    nsf2 = plot->db2->gdb.nscaff-1;
    blf = brt = 0;
    for (i = 1; i < nsf2; i++)
      { y8 = ctg2[scf2[i].fctg].sbeg;
        y  = (int) (y8*ya+yb);
        if (y >= 22)
          { if (y <= rectH+22)
              { brt = i;
                // painter.drawLine(cxb,y,cxe,y);
              }
          }
        else
          blf = i;
      }
    if (brt < blf)
      brt = blf;

    if ((art-alf <= 2 || scf1[art].ectg - scf1[alf].fctg <= 20) &&
        (brt-blf <= 2 || scf2[brt].ectg - scf2[blf].fctg <= 20) )
      { iPen.setColor(QColor(200,200,255));
        iPen.setWidth(1);
        iPen.setDashPattern(cdash);
        painter.setPen(iPen);

        for (j = alf; j <= art; j++)
          for (i = scf1[j].fctg+1; i < scf1[j].ectg; i++)
            { x8 = ctg1[i].sbeg;
              x  = (int) (x8*xa+xb);
              if (x >= 22 && x <= rectW-22)
                painter.drawLine(x,cyb,x,cye);
            }
        for (j = blf; j <= brt; j++)
          for (i = scf2[j].fctg+1; i < scf2[j].ectg; i++)
            { y8 = ctg2[i].sbeg;
              y  = (int) (y8*ya+yb);
              if (y >= 22 && y <= rectH+22)
                painter.drawLine(cxb,y,cxe,y);
            }

        iPen.setColor(QColor(255,255,255));
        iPen.setWidth(1);
        iPen.setDashPattern(sdash);
        painter.setPen(iPen);
      }
    else
      { iPen.setColor(QColor(255,255,255));
        iPen.setWidth(1);
        iPen.setStyle(Qt::DashLine);
        painter.setPen(iPen);
      }

    for (i = 1; i < nsf1; i++)
      { x8 = ctg1[scf1[i].fctg].sbeg;
        x  = (int) (x8*xa+xb);
        if (x >= 22 && x <= rectW-22)
          painter.drawLine(x,cyb,x,cye);
      }

    for (i = 1; i < nsf2; i++)
      { y8 = ctg2[scf2[i].fctg].sbeg;
        y  = (int) (y8*ya+yb);
        if (y >= 22 && y <= rectH+22)
          painter.drawLine(cxb,y,cxe,y);
      }
  }

  { QuadLeaf *list;
    QPen      pPen;
    int       j, k;
    int       Thickint[5] = { 0, 1, 0, 2, 3 };
    qreal     Thickreal[5] = { .5, 0, 1.5, 0, 0 };

    pPen.setBrush(QBrush(QColor(255,255,255),Qt::Dense2Pattern));
    pPen.setWidth(5);

    for (j = 0; j < state->nlays; j++)
      { k = state->order[j];
        if ( ! state->on[k])
          continue;

        if (k == 0)
          { int   r, g, b;
            int   kmer, klen;
            Dots *dot;

            if (state->view.w > 1000000)
              continue;

            kmer = state->thick[0]+8;
            klen = kmer*xa;

            if (image == NULL || rectW != image->width() || rectH != image->height() ||
                (image->format() == QImage::Format_MonoLSB) == (klen >= 3))
              { if (image != NULL)
                  delete image;
                if (klen >= 3)
                  image = new QImage(rectW,rectH,QImage::Format_RGB32);
                else
                  { image = new QImage(rectW,rectH,QImage::Format_MonoLSB);
                    for (r = 0; r < rectH; r++)
                      raster[r] = image->scanLine(r); 
                  }
              }                  

            if (klen < 3)
              image->setColorTable(ctable);
            image->fill(0);

            dot  = dotplot(plot,kmer,&(state->view));

            { int  i, x;
              int *aplot = dot->aplot;

              for (i = 0; i < dot->ahit; i++)
                { x = aplot[i];
                  if (x >= 0)
                    aplot[i] = ((int) floor(xa*x+22.));
                }
            }

            r = state->colorF[0].red();
            g = state->colorF[0].green();
            b = state->colorF[0].blue();

            if (klen >= 3)
              { QPainter dotter(image);
                QPen     dPen;
                int     *aplot = dot->aplot;
                Tuple   *blist = dot->blist;
                int      i, x, y, k;

                dPen.setColor(QColor(r,g,b));
                dPen.setWidth(1);

                dotter.setRenderHint(QPainter::Antialiasing,true);
                dotter.setClipRegion(QRect(22,22,rectW-22,rectH-22));
                dotter.setPen(dPen);

                for (i = 0; i < dot->brun; i++)
                  { y = (int) floor(ya*blist[i].pos+22.);
                    k = blist[i].code;
                    while (1)
                      { x = aplot[k++];
                        if (x < 0)
                          break;
                        dotter.drawLine(x,y,x+klen,y+klen);
                      }
                  }
              }
            else
              { int i, x, k; 
                uint8 *ras;
                int   *aplot = dot->aplot;
                Tuple *blist = dot->blist;

                ctable[0] = qRgb(0,0,0);
                ctable[1] = qRgb(r,g,b);
                for (r = 0; r < 8; r++)
                  imbit[r] = (1<<r); 

                for (i = 0; i < dot->brun; i++)
                  { ras = raster[((int) floor(ya*blist[i].pos+22.))];
                    k = blist[i].code;
                    while (1)
                      { x = aplot[k++];
                        if (x < 0)
                          break;
                        ras[x>>3] |= imbit[x&0x7];
                      }
                  }
              }

            painter.drawImage(QPoint(0,0),*image);

            continue;
          }

        list = Plot_Layer(plot,k,&frame);

        { QPen        fPen, rPen;
          QuadLeaf   *curn, *next;
          DotSegment *segs, *line;
          int         xbeg, xend;
          int         ybeg, yend;
          int         i;
    
          fPen.setColor(state->colorF[k]);
          rPen.setColor(state->colorR[k]);
          i = state->thick[k];
          if (i == 0 || i == 2)
            { fPen.setWidthF(Thickreal[i]);
              rPen.setWidthF(Thickreal[i]);
            }
          else
            { fPen.setWidth(Thickint[i]);
              rPen.setWidth(Thickint[i]);
            }

          segs = plot->layers[k]->segs;
          for (curn = list; curn != NULL; curn = next)
            { next = (QuadLeaf *) (((QuadNode *) curn)[BLK_SIZE].quads[0]);
              while (curn->length > 0)
                { for (i = 0; i < curn->length; i++)
                    { line = segs + curn->idx[i];
                      xbeg = (int) (line->abeg*xa+xb);
                      ybeg = (int) (line->bbeg*ya+yb);
                      xend = (int) (line->aend*xa+xb);
                      yend = (int) (line->bend*ya+yb);
                      if (line == pickedSeg)
                        { painter.setPen(pPen);
                          painter.drawLine(xbeg,ybeg,xend,yend);
                        }
                      if (line->bbeg < line->bend)
                        painter.setPen(fPen);
                      else
                        painter.setPen(rPen);
                      painter.drawLine(xbeg,ybeg,xend,yend);
                      line->mark = 0;
                    }
                  curn += 1;
                }
            }
        }

        Free_List(list);
      }
  }

  if (state->lViz)
    if (vW < plot->alen*.7 || vH < plot->blen*.7)
      { int64 lX, lY;
        int64 lW, lH;
	int64 x, y, w, h;

        painter.setPen(state->lColor);

	if (plot->alen > plot->blen)
	  { lW = LOCATOR_RECTANGLE_SIZE;
            lH = LOCATOR_RECTANGLE_SIZE*((1.*plot->blen)/plot->alen);
          }
        else
          { lH = LOCATOR_RECTANGLE_SIZE;
            lW = LOCATOR_RECTANGLE_SIZE*((1.*plot->alen)/plot->blen);
          }
            
        if (state->lQuad%2)
          lX = (rectW-20) - (lW + 6);
        else
          lX = 26;
        if (state->lQuad < 2)
          lY = 26;
        else
          lY = (rectH-20) - (lH + 6);

        painter.drawRect(lX,lY,lW,lH);

        if (vX < 0)
          { x = 0;
            w = lW;
          }
        else
          { x = lW*(vX/plot->alen);
            w = lW*(vW/plot->alen); 
          }
        if (vY < 0)
          { y = 0;
            h = lH;
          }
        else
          { y = lH*(vY/plot->blen);
            h = lH*(vH/plot->blen); 
          }

        if (w >= 2 && h >= 2)
          painter.drawRect(lX+x,lY+y,w,h);
        else
          { x += lX+w/2;
            y += lY+h/2;
            painter.drawLine(x-2,y,x+2,y);
            painter.drawLine(x,y-2,x,y+2);
          }
      }

  if (state->fOn)
    { painter.setPen(state->fColor);
      int x = xa*state->focus.x + xb;
      int y = ya*state->focus.y + yb;
      if (state->fViz)
        { if (xb-2 < 22)
            if (xa*plot->alen+xb+2 > rectW-22)
              painter.drawLine(22,y,rectW-22,y);
            else
              painter.drawLine(22,y,xa*plot->alen+xb+2,y);
          else
            if (xa*plot->alen+xb+2 > rectW-22)
              painter.drawLine(xb-2,y,rectW-22,y);
            else
              painter.drawLine(xb-2,y,xa*plot->alen+xb+2,y);
          if (yb-2 < 22)
            if (ya*plot->blen+yb+2 > rectH-22)
              painter.drawLine(x,22,x,rectH-22);
            else
              painter.drawLine(x,22,x,ya*plot->blen+yb+2);
          else
            if (ya*plot->blen+yb+2 > rectH-22)
              painter.drawLine(x,yb-2,x,rectH-22);
            else
              painter.drawLine(x,yb-2,x,ya*plot->blen+yb+2);
        }
      else
        { painter.drawLine(x-10,y,x+10,y);
          painter.drawLine(x,y-10,x,y+10);
        }
    }

  painter.end();
}



/*****************************************************************************\
*
*  DRAG AND DROP LAYER WIDGET
*
\*****************************************************************************/

LayerWidget::LayerWidget(DotState *istate, DotCanvas *icanvas, QWidget *parent) : QWidget(parent)
{ setAcceptDrops(true);
  state  = istate;
  canvas = icanvas;
}

void LayerWidget::mousePressEvent(QMouseEvent *ev)
{
  QWidget *child = static_cast<QWidget *>(childAt(ev->pos()));

  if (child == NULL)
    return;
  if (children().contains(child))
    return;

  QWidget *list = static_cast<QWidget *>(child->parent());
  if (child != list->children().at(1))
    return;

  QPixmap line(100,5);
    QPainter pnt(&line);
    QPen     pen(Qt::darkBlue,5);
    pnt.setPen(pen);
    pnt.drawLine(0,2,100,2);

  QDrag *drag = new QDrag(this);

  QMimeData *mimeData = new QMimeData;
    drag->setMimeData(mimeData);

  drag->setPixmap(line);  // list->grab()

  QPoint hotSpot = ev->pos() - list->pos();
    hotSpot.setY(hotSpot.y()-list->height()/2);
    drag->setHotSpot(hotSpot);

  bool wasOn = list->isEnabled();
  list->setEnabled(false);

  if (drag->exec() == Qt::IgnoreAction)
    { list->setEnabled(wasOn);
      return;
   }

  { int w, f, t;
    int x, j;
    QVBoxLayout *lman = static_cast<QVBoxLayout *>(layout());

    w = lman->itemAt(1)->widget()->pos().y();
    f = list->pos().y()/w;
    t = (dropY - hotSpot.y() + w/2) / w;

    if (t < f || t > f+1)
      { lman->removeWidget(list);
        if (t < f)
          lman->insertWidget(t,list);
        else
          lman->insertWidget(t-1,list);
        x = state->order[f];
        if (t < f)
          { for (j = f; j > t; j--)
              state->order[j] = state->order[j-1];
            state->order[t] = x;
          }
        else
          { for (j = f; j < t-1; j++)
              state->order[j] = state->order[j+1];
            state->order[t-1] = x;
          }
      }
  }

  list->setEnabled(wasOn);

  canvas->update();
}

void LayerWidget::dragEnterEvent(QDragEnterEvent *ev)
{ ev->setDropAction(Qt::MoveAction);
  ev->accept();
}

void LayerWidget::dragMoveEvent(QDragMoveEvent *ev)
{ ev->setDropAction(Qt::MoveAction);
  ev->accept();
}

void LayerWidget::dropEvent(QDropEvent *ev)
{ dropY = ev->position().toPoint().y();
  ev->accept();
}

void DotWindow::activateLayer(int act)
{ bool on;
  int  j;

  (void) act;

  for (j = 0; j < state.nlays; j++)
    { state.on[j] = on = layerOn[j]->isChecked();
      layerTitle[j]->setEnabled(on);
      layerFBox[j]->setEnabled(on);
      layerFText[j]->setEnabled(on);
      layerRBox[j]->setEnabled(on);
      layerRText[j]->setEnabled(on);
      layerThick[j]->setEnabled(on);
    }

  update();
}

void DotWindow::layerFChange()
{ int j;
  for (j = 0; j < state.nlays; j++)
    if (layerFBox[j]->isDown())
      break;
  QColor newColor = QColorDialog::getColor(state.colorF[j],this);
  layerFBox[j]->setDown(false);
  if ( ! newColor.isValid()) return;

  state.colorF[j] = newColor;
  QPixmap blob = QPixmap(16,16);
    blob.fill(newColor);
  layerFBox[j]->setIcon(QIcon(blob));
}

void DotWindow::layerRChange()
{ int j;

  for (j = 0; j < state.nlays; j++)
    if (layerRBox[j]->isDown())
      break;
  QColor newColor = QColorDialog::getColor(state.colorR[j],this);
  layerRBox[j]->setDown(false);
  if ( ! newColor.isValid()) return;

  state.colorR[j] = newColor;
  QPixmap blob = QPixmap(16,16);
    blob.fill(newColor);
  layerRBox[j]->setIcon(QIcon(blob));
}


/***************************************************************************************/
/*                                                                                     */
/*   ALIGNMENT VIEWER                                                                  */
/*                                                                                     */
/***************************************************************************************/

AlignWindow::AlignWindow(QString title, char *align, QWidget *parent, Qt::WindowFlags flags) : QMainWindow(parent,flags)
{
  edit = new QTextEdit();
  edit->setCurrentFont(QFont(tr("Monaco"),11));
  edit->setPlainText(tr(align));
  edit->setReadOnly(true);
  edit->setWordWrapMode(QTextOption::NoWrap);
  edit->setLineWrapMode(QTextEdit::NoWrap);

  setCentralWidget(edit);
  setWindowTitle(title);

  setMinimumWidth(820);
  setMinimumHeight(500);

  raise();
}


/***************************************************************************************/
/*                                                                                     */
/*   DOT WINDOW DATA TYPE                                                              */
/*                                                                                     */
/***************************************************************************************/

OpenDialog *DotWindow::openDialog;
Open_State  DotWindow::dataset;

QList<DotWindow *> DotWindow::dotwindows;

void DotWindow::openCopy()
{ DotPlot *nplot;

  this->state.wGeom = geometry();

  nplot = copyPlot(this->plot);

  DotWindow *dot = new DotWindow(nplot,&this->state,1);
  dot->show();
  dot->move(QPoint(state.wGeom.x()+windowHeight,state.wGeom.y()));
  dot->raise();

  dotwindows += dot;
}

void DotWindow::openFile()    // static
{ Open_State ostate;
  DotState  *sptr;
  DotPlot   *plot;

tryagain:

  openDialog->getState(ostate);
  if (openDialog->exec() == QDialog::Accepted)
    openDialog->getState(dataset);
  else
    { if (dotwindows.length() == 0)
        { QApplication::quit();
          exit (0);
        }
      openDialog->putState(ostate);
      return;
    }

  plot = createPlot(dataset.alnInfo->absoluteFilePath().toLatin1().data(),
                    dataset.longCut,dataset.idCut,dataset.sizeCut,NULL);
  if (plot == NULL)
    { DotWindow::warning(tr(Ebuffer),NULL,DotWindow::ERROR,tr("OK"));
      if (dotwindows.length() == 0)
        goto tryagain;
      else
        return;
    }
  if (plot->db1->gdb.seqs == NULL)
    DotWindow::warning(tr("Could not find source for %1, alignments & dot-plot disabled").arg(plot->db1->name),NULL,DotWindow::ERROR,tr("OK"));
  else if (plot->db1 != plot->db2 && plot->db2->gdb.seqs == NULL)
    DotWindow::warning(tr("Could not find source for %1, alignments & dot-plot disabled").arg(plot->db2->name),NULL,DotWindow::ERROR,tr("OK"));

  if (dotwindows.length() == 0)
    sptr = NULL;
  else
    sptr = & dotwindows.at(dotwindows.length()-1)->state;
  DotWindow *dot = new DotWindow(plot,sptr,0);
  dot->raise();
  dot->show();

  dotwindows += dot;
}

void DotWindow::openOverlay()
{ Open_State ostate;
  DotPlot   *nplot;

  openDialog->getState(ostate);
  if (openDialog->exec() == QDialog::Accepted)
    openDialog->getState(dataset);
  else
    { openDialog->putState(ostate);
      return;
    }

  nplot = createPlot(dataset.alnInfo->absoluteFilePath().toLatin1().data(),
                     dataset.longCut,dataset.idCut,dataset.sizeCut,plot);
  if (nplot == NULL)
    { DotWindow::warning(tr(Ebuffer),this,DotWindow::ERROR,tr("OK"));
      return;
    }

  state.nlays = plot->nlays;
  layerTitle[state.nlays-1]->setText(tr(plot->layers[state.nlays-1]->name));
  layerWidget[state.nlays-1]->setVisible(true);
}

DotWindow::~DotWindow()
{ if (plot != NULL)
    Free_DotPlot(plot);
  printf("Deleting %p\n",plot); fflush(stdout);
}

DotWindow::DotWindow(DotPlot *model, DotState *startState, bool isCopy) : QMainWindow()
{ int j;

  // setAttribute(Qt::WA_DeleteOnClose);

  if (screenGeometry == NULL)
    screenGeometry = new QRect(screen()->availableGeometry());

  plot = model;

  canvas = new DotCanvas(this);

  frame = canvas->shareData(&state,plot);

      QLabel *zoomLabel = new QLabel(tr("Zoom: "));

      QToolButton *zoomIn = new QToolButton();
        zoomIn->setToolButtonStyle(Qt::ToolButtonTextOnly);
        zoomIn->setText(tr("+"));
        zoomIn->setFixedSize(26,26);
        zoomIn->setToolTip(tr("Click to zoom in a notch (sqrt(2))"));

      QDoubleValidator *zVal = new QDoubleValidator(0.,DBL_MAX,2,this);

      zoomEdit = new DotLineEdit(this);
        zoomEdit->setFixedHeight(24);
        zoomEdit->setFixedWidth(73);
        zoomEdit->setFrame(false);
        zoomEdit->setFont(QFont(tr("Monaco")));
        zoomEdit->setAlignment(Qt::AlignCenter);
        zoomEdit->setValidator(zVal);

      QToolButton *zoomOut = new QToolButton();
        zoomOut->setToolButtonStyle(Qt::ToolButtonTextOnly);
        zoomOut->setText(tr("-"));
        zoomOut->setFixedSize(26,26);
        zoomOut->setToolTip(tr("Click to zoom out a notch (sqrt(2))"));

      QToolButton *zoom0 = new QToolButton();
        zoom0->setIconSize(QSize(16,16));
        zoom0->setFixedSize(26,26);
        zoom0->setIcon(QIcon(":/images/backarrow.png"));
        zoom0->setToolTip(tr("Click to return to previous zoom level"));

  QHBoxLayout *zoomLayout = new QHBoxLayout;
    zoomLayout->addWidget(zoomLabel);
    zoomLayout->addSpacing(10);
    zoomLayout->addWidget(zoomOut);
    zoomLayout->addWidget(zoomEdit);
    zoomLayout->addWidget(zoomIn);
    zoomLayout->addSpacing(10);
    zoomLayout->addWidget(zoom0);
    zoomLayout->addStretch(1);
    zoomLayout->setSpacing(0);

      QLabel *Alabel = new QLabel(tr("X: "));

      Arng = new DotLineEdit(this);
        Arng->setFixedHeight(24);
        Arng->setFrame(false);
        Arng->setFont(QFont(tr("Monaco")));
        Arng->setAlignment(Qt::AlignLeft);

  QHBoxLayout *viewLayoutA = new QHBoxLayout;
    viewLayoutA->addWidget(Alabel,0);
    viewLayoutA->addWidget(Arng,1);
    viewLayoutA->setSpacing(0);

      QLabel *Blabel = new QLabel(tr("Y: "));

      Brng = new DotLineEdit(this);
        Brng->setFixedHeight(24);
        Brng->setFrame(false);
        Brng->setFont(QFont(tr("Monaco")));
        Brng->setAlignment(Qt::AlignLeft);

  QHBoxLayout *viewLayoutB = new QHBoxLayout;
    viewLayoutB->addWidget(Blabel,0);
    viewLayoutB->addWidget(Brng,1);
    viewLayoutB->setSpacing(0);

  QToolButton *ABpush = new QToolButton();
    ABpush->setToolButtonStyle(Qt::ToolButtonTextOnly);
    ABpush->setText(tr("="));
    ABpush->setFixedSize(26,26);
    ABpush->setFocusPolicy(Qt::NoFocus);
    ABpush->setToolTip(tr("Click to set view to A- & B-ranges"));

      QLabel *posLabel = new QLabel(tr("Focus:  "));

      Fpnt = new DotLineEdit(this);
        Fpnt->setFixedHeight(24);
	Fpnt->setFrame(false);
        Fpnt->setFont(QFont(tr("Monaco")));
        Fpnt->setAlignment(Qt::AlignLeft);

  QHBoxLayout *cursorLayout = new QHBoxLayout;
    cursorLayout->addWidget(posLabel,0);
    cursorLayout->addWidget(Fpnt,1);
    cursorLayout->setSpacing(0);

  QToolButton *Fpush = new QToolButton();
    Fpush->setToolButtonStyle(Qt::ToolButtonTextOnly);
    Fpush->setText(tr("="));
    Fpush->setFixedSize(26,26);
    Fpush->setToolTip(tr("Click to set center of view to focus coords"));

          cFormat = new QComboBox();
            cFormat->addItem(tr("#"));
            cFormat->addItem(tr(".c:#"));
            cFormat->addItem(tr("@s:#"));
            cFormat->addItem(tr("@s.c:#"));
            cFormat->addItem(tr("@id:#"));
            cFormat->addItem(tr("@id.c:#"));

          QLabel *formatLabel = new QLabel(tr("Format: "));

  QHBoxLayout *formatLayout = new QHBoxLayout;
    formatLayout->addWidget(formatLabel);
    formatLayout->addSpacing(6);
    formatLayout->addWidget(cFormat);
    formatLayout->addStretch(1);
    formatLayout->setSpacing(0);

      QLabel *focusLabel = new QLabel(tr("Focus:  "));

      focusOn = new QCheckBox(tr("On"));

      focusBox = new QToolButton();
        focusBox->setIconSize(QSize(16,16));
        focusBox->setFixedSize(20,20);

      focusCheck = new QCheckBox(tr("Cross Hairs"));
    
  QHBoxLayout *focusLayout = new QHBoxLayout;
    focusLayout->addWidget(focusLabel);
    focusLayout->addWidget(focusOn);
    focusLayout->addSpacing(12);
    focusLayout->addWidget(focusBox);
    focusLayout->addSpacing(10);
    focusLayout->addWidget(focusCheck);
    focusLayout->addStretch(1);
    focusLayout->setSpacing(0);

  QPixmap upd = QPixmap(tr(":/images/UpDown.png")).
                    scaled(16,16,Qt::IgnoreAspectRatio,Qt::SmoothTransformation);

  for (j = 0; j < MAX_LAYERS; j++)
    { layerOn[j] = new QCheckBox();

      if (j == 0)
        layerTitle[j] = new QLabel(tr("Dot Plot"));
      else
        layerTitle[j] = new QLabel(tr("Name%1").arg(j));

      layerFBox[j] = new QToolButton();
        layerFBox[j]->setIconSize(QSize(16,16));
        layerFBox[j]->setFixedSize(20,20);

      if (j == 0)
        layerFText[j] = new QLabel(tr("     K-mer"));
      else
        layerFText[j] = new QLabel(tr(" F"));

      layerRBox[j] = new QToolButton();
        layerRBox[j]->setIconSize(QSize(16,16));
        layerRBox[j]->setFixedSize(20,20);

      layerRText[j] = new QLabel(tr(" R"));

      layerThick[j] = new QComboBox();
      if (j == 0)
        { for (int k = 8; k <= 32; k++)
            layerThick[j]->addItem(tr("%1").arg(k));
        }
      else
        { layerThick[j]->addItem(tr(".5"));
          layerThick[j]->addItem(tr("1"));
          layerThick[j]->addItem(tr("1.5"));
          layerThick[j]->addItem(tr("2"));
          layerThick[j]->addItem(tr("3"));
        }

      QLabel *layb = new QLabel();
        layb->setFixedSize(16,16);
        layb->setPixmap(upd);

      QHBoxLayout *layerLayout1 = new QHBoxLayout();
        layerLayout1->setContentsMargins(0,0,0,0);
        layerLayout1->addSpacing(6);
        layerLayout1->addWidget(layerOn[j]);
        layerLayout1->addSpacing(9);
        layerLayout1->addWidget(layerTitle[j]);
        layerLayout1->addStretch(1);

      QHBoxLayout *layerLayout2 = new QHBoxLayout();
        layerLayout2->setContentsMargins(0,0,0,0);
        layerLayout2->setSpacing(0);
        layerLayout2->addSpacing(24);
        layerLayout2->addWidget(layerFBox[j]);
        layerLayout2->addWidget(layerFText[j]);
        if (j == 0)
          layerLayout2->addSpacing(2);
        else
          { layerLayout2->addSpacing(6);
            layerLayout2->addWidget(layerRBox[j]);
            layerLayout2->addWidget(layerRText[j]);
            layerLayout2->addSpacing(6);
          }
	layerLayout2->addWidget(layerThick[j]);
	layerLayout2->addStretch(1);

      QVBoxLayout *layerLayout3 = new QVBoxLayout();
        layerLayout3->setContentsMargins(0,0,0,0);
        layerLayout3->setSpacing(0);
        layerLayout3->addLayout(layerLayout1);
        layerLayout3->addLayout(layerLayout2);

      QHBoxLayout *layerLayout = new QHBoxLayout();
        layerLayout->setContentsMargins(0,0,0,0);
        layerLayout->setSpacing(0);
        layerLayout->addWidget(layb);
        layerLayout->addLayout(layerLayout3);

      layerWidget[j] = new QWidget();
        layerWidget[j]->setLayout(layerLayout);
    }

  QVBoxLayout *layerLayout = new QVBoxLayout();
    layerLayout->setContentsMargins(10,0,0,0);
    layerLayout->setSpacing(0);
    layerLayout->setSizeConstraint(QLayout::SetFixedSize);
    for (j = 0; j < MAX_LAYERS; j++)
      layerLayout->addWidget(layerWidget[j]);

  layerPanel = new LayerWidget(&state,canvas);
    layerPanel->setLayout(layerLayout);

  QScrollArea *layerArea = new QScrollArea();
    layerArea->setWidget(layerPanel);
    layerArea->setAlignment(Qt::AlignLeft|Qt::AlignTop);
    layerArea->setWidgetResizable(false);

  QLabel *layerMaster = new QLabel(tr("Layers:"));

  QVBoxLayout *layerMargin = new QVBoxLayout();
    layerMargin->setContentsMargins(2,15,2,2);
    layerMargin->addWidget(layerMaster);
    layerMargin->addWidget(layerArea);

      QLabel *locatorLabel = new QLabel(tr("Navi: "));

      locatorBox = new QToolButton();
        locatorBox->setIconSize(QSize(16,16));
        locatorBox->setFixedSize(20,20);

      QToolButton *locatorTL = new QToolButton();
        locatorTL->setIconSize(QSize(6,6));
        locatorTL->setFixedSize(10,10);
        locatorTL->setCheckable(true);
        locatorTL->setIcon(QIcon(":/images/topleft.png"));

      QToolButton *locatorTR = new QToolButton();
        locatorTR->setIconSize(QSize(6,6));
        locatorTR->setFixedSize(10,10);
        locatorTR->setCheckable(true);
        locatorTR->setIcon(QIcon(":/images/topright.png"));

      QToolButton *locatorBL = new QToolButton();
        locatorBL->setIconSize(QSize(6,6));
        locatorBL->setFixedSize(10,10);
        locatorBL->setCheckable(true);
        locatorBL->setIcon(QIcon(":/images/bottomleft.png"));

      QToolButton *locatorBR = new QToolButton();
        locatorBR->setIconSize(QSize(6,6));
        locatorBR->setFixedSize(10,10);
        locatorBR->setCheckable(true);
        locatorBR->setIcon(QIcon(":/images/bottomright.png"));

      locatorQuad = new QButtonGroup;
        locatorQuad->addButton(locatorTL,0);
        locatorQuad->addButton(locatorTR,1);
        locatorQuad->addButton(locatorBL,2);
        locatorQuad->addButton(locatorBR,3);
        locatorQuad->setExclusive(true);

      QGridLayout *directLayout = new QGridLayout;
        directLayout->addWidget(locatorTL,0,0);
        directLayout->addWidget(locatorTR,0,1);
        directLayout->addWidget(locatorBL,1,0);
        directLayout->addWidget(locatorBR,1,1);
        directLayout->setSpacing(0);
        directLayout->setContentsMargins(0,0,0,0);

      locatorCheck = new QCheckBox(tr("On"));

  QHBoxLayout *locatorLayout = new QHBoxLayout;
    locatorLayout->addWidget(locatorLabel);
    locatorLayout->addSpacing(6);
    locatorLayout->addWidget(locatorCheck);
    locatorLayout->addWidget(locatorBox);
    locatorLayout->addLayout(directLayout);
    locatorLayout->addStretch(1);
    locatorLayout->setContentsMargins(0,0,0,0);

  QLabel *panel = new QLabel();
    panel->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    panel->setFrameStyle(QFrame::StyledPanel|QFrame::Plain);
    panel->setLineWidth(1);
    panel->setAutoFillBackground(true);
    panel->setText(tr("%1<br>%2<br>%3<br>%4<br>%5<br>%6<br>%7")
                         .arg("<b>Over plot area:</b>")
                         .arg("&nbsp;&nbsp;Click to set focus")
                         .arg("&nbsp;&nbsp;Press-Drag to pan")
                         .arg("&nbsp;&nbsp;Shift-Press-Drag for zoom area")
                         .arg("&nbsp;&nbsp;x-Press for temporary pick")
                         .arg("&nbsp;&nbsp;X-Press for fixed pick")
                         .arg("&nbsp;&nbsp;z, +, - for zoom control"));
    panel->setAlignment(Qt::AlignLeft|Qt::AlignTop);

      QPalette panelColor = panel->palette();
        panelColor.setColor(QPalette::Active,QPalette::Window,QColor(200,200,220));
        panel->setPalette(panelColor);

  QVBoxLayout *controlLayout = new QVBoxLayout;
    controlLayout->addStrut(120);
    controlLayout->addLayout(zoomLayout);
    controlLayout->addSpacing(30);
    controlLayout->addLayout(focusLayout);
    controlLayout->addSpacing(23);
    controlLayout->addLayout(formatLayout);
    controlLayout->addSpacing(15);
    controlLayout->addLayout(layerMargin);
    controlLayout->addSpacing(5);
    controlLayout->addLayout(locatorLayout);
    controlLayout->addWidget(panel);

  QWidget *controlPart = new QWidget;
    controlPart->setLayout(controlLayout);
    controlPart->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Expanding);

  QHBoxLayout *toolLayout = new QHBoxLayout;
    toolLayout->addWidget(canvas,1);
    toolLayout->addWidget(controlPart,0);
    toolLayout->setSpacing(0);
    toolLayout->setContentsMargins(0,0,0,0);

  QHBoxLayout *viewLayout = new QHBoxLayout;
    viewLayout->setSpacing(0);
    viewLayout->addSpacing(20);
    viewLayout->addLayout(viewLayoutA,1);
    viewLayout->addSpacing(10);
    viewLayout->addLayout(viewLayoutB,1);
    viewLayout->addSpacing(10);
    viewLayout->addWidget(ABpush);
    viewLayout->addSpacing(20);
    viewLayout->addLayout(cursorLayout,1);
    viewLayout->addSpacing(10);
    viewLayout->addWidget(Fpush);
    viewLayout->addSpacing(20);
    viewLayout->setContentsMargins(0,5,0,5);

  QVBoxLayout *rangeLayout = new QVBoxLayout;
    rangeLayout->addLayout(toolLayout,1);
    rangeLayout->addLayout(viewLayout,0);
    rangeLayout->setSpacing(0);
    rangeLayout->setContentsMargins(0,0,0,0);

  QWidget *centralFrame = new QWidget();
    centralFrame->setLayout(rangeLayout);

  setCentralWidget(centralFrame);

  setWindowTitle(tr("ALNview: %1 vs %2").arg(plot->db1->name).arg(plot->db2->name));
  setMinimumWidth(800);
  setMinimumHeight(600);

  createActions();
  createMenus();
  createToolBars();

  if (startState == NULL)
    { readAndApplySettings();
      state.lMag = plot->alen;
      state.lXct = plot->alen/2.;
      state.lYct = plot->blen/2.;
      state.zMag.empty();
      state.zXct.empty();
      state.zYct.empty();
      for (j = 0; j < MAX_LAYERS; j++)
        { state.order[j] = j;
          state.on[j] = true;
        }
      if (plot->db1->gdb.seqs == NULL || plot->db2->gdb.seqs == NULL)
        { state.on[0] = false;
          layerOn[0]->setEnabled(false);
        }
    }
  else
    { state.wGeom  = startState->wGeom;
      state.view   = startState->view;
      state.zoom   = startState->zoom;
      state.format = startState->format;
      state.fOn    = startState->fOn;
      state.focus  = startState->focus;
      state.fColor = startState->fColor;
      state.fViz   = startState->fViz;
      state.lColor = startState->lColor;
      state.lViz   = startState->lViz;
      state.lQuad  = startState->lQuad;
      state.lMag   = startState->lMag;
      state.lXct   = startState->lXct;
      state.lYct   = startState->lYct;
      state.zMag   = startState->zMag;
      state.zXct   = startState->zXct;
      state.zYct   = startState->zYct;
      state.nlays  = startState->nlays;
      for (j = 0; j < MAX_LAYERS; j++)
        { state.order[j]  = startState->order[j];
          state.on[j]     = startState->on[j];
          state.colorF[j] = startState->colorF[j];
          state.colorR[j] = startState->colorR[j];
          state.thick[j] = startState->thick[j];
        }
      if (plot->db1->gdb.seqs == NULL || plot->db2->gdb.seqs == NULL)
        layerOn[0]->setEnabled(false);
    }
  if ( ! isCopy)
    { state.view.x = 0;
      state.view.y = 0;
      state.view.w = plot->alen;
      state.view.h = plot->blen;
      state.focus.x = 0;
      state.focus.y = 0;
    }
  state.nlays = plot->nlays;
  pushState();

  connect(canvas,SIGNAL(NewFrame(double)),this,SLOT(frameToView(double)));
  connect(canvas,SIGNAL(NewFocus()),this,SLOT(clickToFocus()));

  connect(zoomIn,SIGNAL(clicked()),this,SLOT(zoomUp()));
  connect(zoomOut,SIGNAL(clicked()),this,SLOT(zoomDown()));
  connect(zoom0,SIGNAL(clicked()),canvas,SLOT(zoomPop()));
  connect(zoomEdit,SIGNAL(focusOut()),this,SLOT(zoomTo()));

  connect(Arng,SIGNAL(focusOut()),this,SLOT(checkRangeA()));
  connect(Brng,SIGNAL(focusOut()),this,SLOT(checkRangeB()));
  connect(Fpnt,SIGNAL(focusOut()),this,SLOT(checkFocus()));

  connect(ABpush,SIGNAL(clicked()),this,SLOT(viewChange()));
  connect(Fpush,SIGNAL(clicked()),this,SLOT(focusChange()));

  connect(cFormat,SIGNAL(currentIndexChanged(int)),this,SLOT(formatChange(int)));

  connect(focusOn,SIGNAL(stateChanged(int)),this,SLOT(focusOnChange()));
  connect(focusBox,SIGNAL(clicked()),this,SLOT(focusColorChange()));
  connect(focusCheck,SIGNAL(stateChanged(int)),this,SLOT(hairsChange()));

  for (j = 0; j < MAX_LAYERS; j++)
    { connect(layerOn[j],SIGNAL(stateChanged(int)),this,SLOT(activateLayer(int)));
      connect(layerFBox[j],SIGNAL(pressed()),this,SLOT(layerFChange()));
      connect(layerRBox[j],SIGNAL(pressed()),this,SLOT(layerRChange()));
      connect(layerThick[j],SIGNAL(currentIndexChanged(int)),this,SLOT(thickChange(int)));
    }

  connect(locatorCheck,SIGNAL(stateChanged(int)),this,SLOT(locatorChange()));
  connect(locatorBox,SIGNAL(clicked()),this,SLOT(locatorColorChange()));
  connect(locatorTL,SIGNAL(clicked()),this,SLOT(locatorChange()));
  connect(locatorTR,SIGNAL(clicked()),this,SLOT(locatorChange()));
  connect(locatorBL,SIGNAL(clicked()),this,SLOT(locatorChange()));
  connect(locatorBR,SIGNAL(clicked()),this,SLOT(locatorChange()));

  Arng->setFocus();
  Arng->setChain(zoomEdit,Brng);
  Brng->setChain(Arng,Fpnt);
  Fpnt->setChain(Brng,zoomEdit);
  zoomEdit->setChain(Fpnt,Arng);

  windowWidth     = frameGeometry().width() - geometry().width();
  windowHeight    = frameGeometry().height() - geometry().height();
}

void DotWindow::toggleToolBar()
{ int i;

  if (fileToolBar->isVisible())
    { for (i = 0; i < dotwindows.length(); i++)
        { dotwindows[i]->toolArea = dotwindows[i]->toolBarArea(dotwindows[i]->fileToolBar);
          dotwindows[i]->removeToolBar(dotwindows[i]->fileToolBar);
          dotwindows[i]->toolAct->setText(tr("Show Toolbar"));
        }
    }
  else
    { for (i = 0; i < dotwindows.length(); i++)
        { dotwindows[i]->addToolBar(dotwindows[i]->toolArea,dotwindows[i]->fileToolBar);
          dotwindows[i]->fileToolBar->setVisible(true);
          dotwindows[i]->toolAct->setText(tr("Hide Toolbar"));
        }
    }
}

void DotWindow::createActions()
{ openAct = new QAction(QIcon(":/images/open.png"), tr("&Open"), this);
    openAct->setShortcut(tr("Ctrl+O"));
    openAct->setToolTip(tr("Open an overlap data set"));

  overAct = new QAction(QIcon(":/images/overlay.png"), tr("&Overlay"), this);
    overAct->setShortcut(tr("Ctrl+L"));
    overAct->setToolTip(tr("Overlay on this window's view"));

  copyAct = new QAction(QIcon(":/images/copy.png"), tr("&Duplicate"), this);
    copyAct->setShortcut(tr("Ctrl+D"));
    copyAct->setToolTip(tr("Duplicate this window's view"));

  exitAct = new QAction(tr("Exit"), this);
    exitAct->setShortcut(tr("Ctrl+Q"));
    exitAct->setToolTip(tr("Exit the application"));
    exitAct->setMenuRole(QAction::QuitRole);

  unminimizeAllAct = new QAction(tr("&Unminimize All"), this);
    unminimizeAllAct->setShortcut(tr("Ctrl+U"));
    unminimizeAllAct->setToolTip(tr("Un-minimize all windows"));

  raiseAllAct = new QAction(tr("Raise all windows"), this);
    raiseAllAct->setShortcut(tr("Ctrl+R"));
    raiseAllAct->setToolTip(tr("Bring all windows to front"));

  connect(openAct, SIGNAL(triggered()), this, SLOT(openFile()));
  connect(overAct, SIGNAL(triggered()), this, SLOT(openOverlay()));
  connect(copyAct, SIGNAL(triggered()), this, SLOT(openCopy()));
  connect(exitAct, SIGNAL(triggered()), this, SLOT(closeAll()));

  connect(unminimizeAllAct, SIGNAL(triggered()), this, SLOT(unminimizeAll()));
  connect(raiseAllAct, SIGNAL(triggered()), this, SLOT(raiseAll()));

  toolAct = new QAction(tr("Hide Toolbar"), this);
    toolAct->setToolTip(tr("Show/Hide toolbars in all pile windows"));

  tileAct = new QAction(QIcon(":/images/tile.png"), tr("&Tile"), this);
    tileAct->setShortcut(tr("Ctrl+T"));
    tileAct->setToolTip(tr("Tile all open image windows"));

  cascadeAct = new QAction(QIcon(":/images/cascade.png"), tr("&Cascade"), this);
    cascadeAct->setShortcut(tr("Ctrl+C"));
    cascadeAct->setToolTip(tr("Cascade all open image windows"));

  connect(toolAct, SIGNAL(triggered()), this, SLOT(toggleToolBar()));
  connect(tileAct,SIGNAL(triggered()),this,SLOT(tileImages()));
  connect(cascadeAct,SIGNAL(triggered()),this,SLOT(cascadeImages()));
}

void DotWindow::createToolBars()
{ fileToolBar = addToolBar(tr("File ToolBar"));
    fileToolBar->addAction(overAct);
    fileToolBar->addAction(copyAct);
    fileToolBar->addAction(tileAct);
    fileToolBar->addAction(cascadeAct);
  toolArea = Qt::TopToolBarArea;
}

void DotWindow::createMenus()
{ QMenuBar *bar = menuBar();

  QMenu *fileMenu = bar->addMenu(tr("&File"));
    fileMenu->addAction(openAct);
    fileMenu->addSeparator();
    fileMenu->addAction(exitAct);

  QMenu *imageMenu = bar->addMenu(tr("&Image"));
    imageMenu->addAction(tileAct);
    imageMenu->addAction(cascadeAct);
    imageMenu->addAction(toolAct);
 
  QMenu *windowMenu = bar->addMenu(tr("&Window"));
    windowMenu->addAction(unminimizeAllAct);
    windowMenu->addAction(raiseAllAct);
}

void DotWindow::frameToView(double newZ)
{ char *s;

  int64 ab = frame->x;
  int64 ae = frame->x + frame->w;
  int64 bb = frame->y;
  int64 be = frame->y + frame->h;
  if (ab < 0)
    { ab = 0;
      ae = plot->alen;
    }
  if (bb < 0)
    { bb = 0;
      be = plot->blen;
    }

  state.view.x = ab;
  state.view.y = bb;
  state.view.w = ae-ab;
  state.view.h = be-bb;
  if (newZ > 0.)
    state.zoom = newZ;

  s = Map_Coord(&(plot->db1->gdb),state.view.x,state.view.x+state.view.w,state.format,state.view.w);
  Arng->setText(tr("%1").arg(s));

  s = Map_Coord(&(plot->db2->gdb),state.view.y,state.view.y+state.view.h,state.format,state.view.h);
  Brng->setText(tr("%1").arg(s));

  clickToFocus();

  zoomEdit->setText(tr("%1").arg(state.zoom,0,'f',2));
}

void DotWindow::clickToFocus()
{ char *s1 = Map_Coord(& (plot->db1->gdb),state.focus.x,-1,state.format,state.view.w);
  char *s2 = Map_Coord(& (plot->db2->gdb),-1,state.focus.y,state.format,state.view.h);
  Fpnt->setText(tr("%1,%2").arg(s1).arg(s2));
}

void DotWindow::seqMove(int dir, int id)
{ if (id <= 2)
    { state.focus.x = state.focus.x-dir;
      if (state.focus.x < 0)
        state.focus.x = 0;
      else if (state.focus.x > plot->alen)
        state.focus.x = plot->alen;
    }
  else
    { state.focus.y = state.focus.y-dir;
      if (state.focus.y < 0)
        state.focus.y = 0;
      else if (state.focus.y > plot->blen)
        state.focus.y = plot->blen;
    }
  clickToFocus();
  update();
}

void DotWindow::zoomDown()
{ if (canvas->zoomView(0.70710678118))
    update();
}

void DotWindow::zoomUp()
{ if (canvas->zoomView(1.41421356237))
    update();
}

void DotWindow::zoomTo()
{ double nZoom = zoomEdit->text().toDouble();
  if (canvas->zoomView(nZoom/state.zoom))
    update();
  else
    { zoomEdit->setText(tr("%1").arg(state.zoom));
      DotWindow::warning(tr("Beyond maximal zoom of 1bp per pixel"),this,DotWindow::ERROR,tr("OK"));
    }
}

void DotWindow::checkRangeA()
{ Selection sel;

  if (plot != NULL && interpret_range(&sel,Arng->text().toLatin1().data(),&(plot->db1->gdb),plot->db1->hash))
    DotWindow::warning(tr(EPLACE),this,DotWindow::ERROR,tr("OK"));
}

void DotWindow::checkRangeB()
{ Selection sel;

  if (plot != NULL && interpret_range(&sel,Brng->text().toLatin1().data(),&(plot->db2->gdb),plot->db2->hash))
    DotWindow::warning(tr(EPLACE),this,DotWindow::ERROR,tr("OK"));
}

void DotWindow::viewChange()
{ Selection asel;
  Selection bsel;
  char     *s;
  int64     ab, ae;
  int64     bb, be;

  if (interpret_range(&asel,Arng->text().toLatin1().data(),&(plot->db1->gdb),plot->db1->hash))
    { DotWindow::warning(tr(EPLACE),this,DotWindow::ERROR,tr("OK"));
      return;
    }

  if (interpret_range(&bsel,Brng->text().toLatin1().data(),&(plot->db2->gdb),plot->db2->hash))
    { DotWindow::warning(tr(EPLACE),this,DotWindow::ERROR,tr("OK"));
      return;
    }

  ab = plot->db1->gdb.contigs[asel.c1].sbeg + asel.p1;
  ae = plot->db1->gdb.contigs[asel.c2].sbeg + asel.p2;
  bb = plot->db2->gdb.contigs[bsel.c1].sbeg + bsel.p1;
  be = plot->db2->gdb.contigs[bsel.c2].sbeg + bsel.p2;

  View undo = state.view;
  state.view.x = ab;
  state.view.y = bb;
  state.view.w = ae-ab;
  state.view.h = be-bb;
#ifdef DEBUG
  printf("View Change %d - %d %d - %d\n",ab,ae,bb,be); fflush(stdout);
#endif
  if (canvas->viewToFrame())
    update();
  else
    { state.view = undo;
      DotWindow::warning(tr("Beyond maximal zoom of 1bp per pixel"),this,DotWindow::ERROR,tr("OK"));

      s = Map_Coord(&(plot->db1->gdb),state.view.x,state.view.x+state.view.w,
                   state.format,state.view.w);
      Arng->setText(tr("%1").arg(s));

      s = Map_Coord(&(plot->db2->gdb),state.view.y,state.view.y+state.view.h,
                    state.format,state.view.h);
      Brng->setText(tr("%1").arg(s));
    }
}

void DotWindow::checkFocus()
{ Selection sel;

  if (plot != NULL && interpret_point(&sel,Fpnt->text().toLatin1().data(),&(plot->db1->gdb),plot->db1->hash,&(plot->db2->gdb),plot->db2->hash))
    DotWindow::warning(tr(EPLACE),this,DotWindow::ERROR,tr("OK"));
}

void DotWindow::focusChange()
{ Selection sel;

  if (interpret_point(&sel,Fpnt->text().toLatin1().data(),&(plot->db1->gdb),plot->db1->hash,&(plot->db2->gdb),plot->db2->hash))
    { DotWindow::warning(tr(EPLACE),this,DotWindow::ERROR,tr("OK"));
      return;
    }

  state.focus.x = plot->db1->gdb.contigs[sel.c1].sbeg + sel.p1;
  state.focus.y = plot->db2->gdb.contigs[sel.c2].sbeg + sel.p2;
  update();
}

void DotWindow::formatChange(int index)
{ state.format = index;
  clickToFocus();
  frameToView(state.zoom);
}

void DotWindow::thickChange(int index)
{ int j;

  (void) index;

  for (j = 0; j < state.nlays; j++)
    state.thick[j] = layerThick[j]->currentIndex();
  update();
}

void DotWindow::focusColorChange()
{ state.fColor = QColorDialog::getColor(state.fColor);
  focusBox->setDown(false);

  QPixmap blob = QPixmap(16,16);
     blob.fill(state.fColor);
  focusBox->setIcon(QIcon(blob));

  update();
}

void DotWindow::focusOnChange()
{ state.fOn = (focusOn->checkState() == Qt::Checked);

  bool en = state.fOn;

  Fpnt->setEnabled(en);
  focusBox->setEnabled(en);
  focusCheck->setEnabled(en);

  update();
}

void DotWindow::hairsChange()
{ state.fViz = (focusCheck->checkState() == Qt::Checked);
  update();
}

void DotWindow::locatorChange()
{ state.lViz   = (locatorCheck->checkState() == Qt::Checked);
  state.lQuad  = (LocatorQuad) locatorQuad->checkedId();
  update();
}

void DotWindow::locatorColorChange()
{ state.lColor = QColorDialog::getColor(state.lColor);
  locatorBox->setDown(false);

  QPixmap blob = QPixmap(16,16);
     blob.fill(state.lColor);
  locatorBox->setIcon(QIcon(blob));

  update();
}

void DotWindow::addTextBox(AlignWindow *awin)
{ alignWindows += awin; }

void DotWindow::closeEvent(QCloseEvent *event)
{ int i;

  if (dotwindows.length() == 1)
    writeSettings();

  for (i = 0; i < dotwindows.length(); i++)
    if (dotwindows[i] == this)
      { dotwindows.removeAt(i);
        break;
      }

  for (i = 0; i < this->alignWindows.length(); i++)
    alignWindows[i]->close();

  if (this->plot != NULL)
    Free_DotPlot(this->plot);
  this->plot = NULL;
       
  event->accept();
}

#define SCROLL_PADDING 4

void DotWindow::tileImages()
{ int tileH, tileW, tileSize;

  int screenH = screenGeometry->height();
  int screenW = screenGeometry->width();
  int maximH  = screenH;
  int maximW  = screenW;

  if (dotwindows.length() != 0)
    { tileSize   = dotwindows.length();
      for (int i = 0; i < dotwindows.length(); i++)
        { DotWindow *win = dotwindows[i];
          int h = win->maximumHeight();
          if (h < maximH)
            maximH = h;
          int w = win->maximumWidth();
          if (w < maximW)
            maximW = w;
        }
    }
  else
    return;

  maximH += windowHeight;
  maximW += windowWidth;

  int deltaW  = SCROLL_PADDING + windowWidth;
  int deltaH  = SCROLL_PADDING + windowHeight;
  int minimH  = minimumHeight() + windowHeight;
  int minimW  = minimumWidth() + deltaW;

  { int wmax, hmax, wmin;

    wmax = screenW / minimW;
    hmax = screenH / minimH;
    wmin = (tileSize-1) / hmax + 1;

    if (wmin > wmax)
      { tileW = wmax;
        tileH = hmax;
      }
    else
      { for (int h = 1; h <= hmax; h++)
          { tileH = h;
            tileW = (tileSize-1) / h + 1;
            if ( screenW / tileW - deltaW >= screenH / tileH - deltaH)
              break;
          }
      }
  }

  { int w, h;

    int screenX = screenGeometry->x();
    int screenY = screenGeometry->y();
    int stepY   = screenH / tileH;
    int stepX   = screenW / tileW;
    if (stepX > maximW)
      stepX = maximW;
    if (stepY > maximH)
      stepY = maximH;

    QSize winSize = QSize(stepX-windowWidth,stepY-windowHeight);

    w = 0;
    h = 0;
    for (int t = 0; t < tileSize; t++)
      { dotwindows[t]->move(QPoint(screenX+w*stepX,screenY+h*stepY));
        dotwindows[t]->resize(winSize);
        dotwindows[t]->raise();
        h += 1;
        if (h >= tileH)
          { h = 0;
            w += 1;
            if (w >= tileW)
              w = 0;
          }
      }
    raise();
  }
}

void DotWindow::cascadeImages()
{ int x = screenGeometry->x();
  int y = screenGeometry->y();

  for (int i = 0; i < dotwindows.length(); i++)
    { dotwindows[i]->move(QPoint(x,y));
      dotwindows[i]->raise();
      x += windowHeight;
      y += windowHeight;
      dotwindows[i]->raise();
    }
}

void DotWindow::unminimizeAll()
{ int i;

  for (i = 0; i < dotwindows.length(); i++)
    dotwindows[i]->setWindowState( windowState() & ~ Qt::WindowMinimized);
}

void DotWindow::raiseAll()
{ int i;

  for (i = 0; i < dotwindows.length(); i++)
    dotwindows[i]->raise();
}

void DotWindow::closeAll()
{ int i;

  for (i = dotwindows.length()-1; i >= 0; i--)
    dotwindows[i]->close();
}

void DotWindow::pushState()
{ QVBoxLayout *lman = static_cast<QVBoxLayout *>(layerPanel->layout());
  QWidget     *widget[MAX_LAYERS];
  int          j, k;
  char        *s1, *s2;

  setGeometry(state.wGeom);

  cFormat->setCurrentIndex(state.format);

  s1 = Map_Coord(&(plot->db1->gdb),state.view.x,state.view.x+state.view.w,
                 state.format,state.view.w);
  Arng->setText(tr("%1").arg(s1));

  s1 = Map_Coord(&(plot->db2->gdb),state.view.y,state.view.y+state.view.h,
                 state.format,state.view.h);
  Brng->setText(tr("%1").arg(s1));

  zoomEdit->setText(tr("%1").arg(state.zoom,0,'f',2));

  canvas->viewToFrame();

  focusOn->setCheckState(state.fOn?Qt::Checked:Qt::Unchecked);
  Fpnt->setEnabled(state.fOn);
  focusBox->setEnabled(state.fOn);
  focusCheck->setEnabled(state.fOn);
  s1 = Map_Coord(& (plot->db2->gdb),state.focus.x,-1,state.format,state.view.w);
  s2 = Map_Coord(& (plot->db2->gdb),-1,state.focus.y,state.format,state.view.h);
  Fpnt->setText(tr("%1,%2").arg(s1).arg(s2));

  QPixmap blob = QPixmap(16,16);

  blob.fill(state.fColor);
  focusBox->setIcon(QIcon(blob));
  focusCheck->setCheckState(state.fViz?Qt::Checked:Qt::Unchecked);

  blob.fill(state.lColor);
  locatorBox->setIcon(QIcon(blob));
  locatorQuad->button((int) state.lQuad)->setChecked(true);
  locatorCheck->setCheckState(state.lViz?(Qt::Checked):(Qt::Unchecked));

  for (j = 0; j < MAX_LAYERS; j++)
    { if (j >= state.nlays)
        layerWidget[j]->setVisible(false);
      else if (plot->layers[j] != NULL)
        layerTitle[j]->setText(tr(plot->layers[j]->name));
      layerOn[j]->setCheckState(state.on[j]?(Qt::Checked):(Qt::Unchecked));
      QPixmap blob1 = QPixmap(16,16);
        blob1.fill(state.colorF[j]);
      layerFBox[j]->setIcon(QIcon(blob1));
      QPixmap blob2 = QPixmap(16,16);
        blob2.fill(state.colorR[j]);
      layerRBox[j]->setIcon(QIcon(blob2));
      layerThick[j]->setCurrentIndex(state.thick[j]);
    }
  activateLayer(0);

  lman = static_cast<QVBoxLayout *>(layerPanel->layout());
  for (j = 0; j < state.nlays; j++)
    widget[j] = lman->itemAt(j)->widget();
  for (j = state.nlays-1; j >= 0; j--)
    { k = state.order[j];
      lman->insertWidget(0,widget[k]);
    }
}

void DotWindow::readAndApplySettings()
{ QRgb colF[MAX_LAYERS];
  QRgb colR[MAX_LAYERS];
  int  j;

  QSettings settings("FASTGA", "ALNview");

  if ( ! QFile(settings.fileName()).exists())
    settings.clear();

  settings.beginGroup("window");
    bool   tbarVisible = settings.value("tbarVisible", true).toBool();
    Qt::ToolBarArea tbarArea =
        (Qt::ToolBarArea) settings.value("tbarArea", Qt::TopToolBarArea).toInt();

    state.wGeom  = settings.value("geom", QRect(0,0,800,600)).toRect();
    state.zoom   = settings.value("zoom", 1.0).toFloat();
    state.fOn    = settings.value("fOn", false).toBool();
    state.focus.x = settings.value("focusX", 0 ).toLongLong();
    state.focus.y = settings.value("focusY", 0 ).toLongLong();
    state.format  = settings.value("format", 0 ).toInt();
    QRgb fRGB    = settings.value("fColor", QColor(255,255,0).rgb()).toUInt();
    state.fViz   = settings.value("fViz", false).toBool();
    QRgb lRGB    = settings.value("lColor", QColor(255,0,255).rgb()).toUInt();
    state.lViz   = settings.value("lViz", true).toBool();
    state.lQuad  = (LocatorQuad) settings.value("lQuad", 0).toInt();
    for (j = 0; j < MAX_LAYERS; j++)
      { colF[j] = settings.value(tr("colF%1").arg(j), QColor(0,255,0).rgb()).toUInt();
        colR[j] = settings.value(tr("colR%1").arg(j), QColor(255,0,0).rgb()).toUInt();
        state.thick[j] = settings.value(tr("thick%1").arg(j), 1).toInt();
      }
  settings.endGroup();

  state.fColor.setRgb(fRGB);
  state.lColor.setRgb(lRGB);
  for (j = 0; j < MAX_LAYERS; j++)
    { state.colorF[j].setRgb(colF[j]);
      state.colorR[j].setRgb(colR[j]);
    }

  if (tbarVisible)
    { removeToolBar(fileToolBar);
      addToolBar(tbarArea,fileToolBar);
      fileToolBar->setVisible(true);
      toolAct->setText(tr("Hide Toolbar"));
    }
  else
    { removeToolBar(fileToolBar);
      toolAct->setText(tr("Show Toolbar"));
    }
}

void DotWindow::writeSettings()
{ int j;

  QSettings settings("FASTGA", "ALNview");

  state.wGeom = geometry();

  settings.beginGroup("window");
    if (fileToolBar->isVisible())
      settings.setValue("tbarArea", toolBarArea(fileToolBar));
    else
      settings.setValue("tbarArea", toolArea);
    settings.setValue("tbarVisible", fileToolBar->isVisible());

    settings.setValue("geom", state.wGeom);
    settings.setValue("zoom", state.zoom);
    settings.setValue("fOn", state.fOn);
    settings.setValue("focusX", state.focus.x);
    settings.setValue("focusY", state.focus.y);
    settings.setValue("format", state.format);
    settings.setValue("fColor", state.fColor.rgb());
    settings.setValue("fViz", state.fViz);
    settings.setValue("lColor", state.lColor.rgb());
    settings.setValue("lViz", state.lViz);
    settings.setValue("lQuad", state.lQuad);
    for (j = 0; j < MAX_LAYERS; j++)
      { settings.setValue(tr("colF%1").arg(j), state.colorF[j].rgb());
        settings.setValue(tr("colR%1").arg(j), state.colorR[j].rgb());
        settings.setValue(tr("thick%1").arg(j), state.thick[j]);
      }
  settings.endGroup();

  openDialog->writeSettings(settings);
}

int DotWindow::warning(const QString& message, QWidget *parent, MessageKind kind,
                        const QString& label1, const QString& label2, const QString& label3)
{ QPixmap     *symbol;
  QPushButton *button1 = 0, *button2 = 0, *button3 = 0;

  QMessageBox msg(QMessageBox::Critical, QObject::tr("DaViewer"),
                  message, QMessageBox::NoButton, parent,
                  Qt::Dialog | Qt::FramelessWindowHint);

  msg.setAutoFillBackground(true);
  msg.setPalette(QPalette(QColor(255,255,180,255)));
  if (kind == INFORM)
    symbol = new QPixmap(":/images/inform.png","PNG");
  else if (kind == WARNING)
    symbol = new QPixmap(":/images/warning.png","PNG");
  else //  kind == ERROR
    symbol = new QPixmap(":/images/error.png","PNG");
  msg.setIconPixmap(*symbol);

  if (label1.isEmpty())
    button1 = msg.addButton(QObject::tr("OK"),QMessageBox::AcceptRole);
  else
    { button1 = msg.addButton(label1,QMessageBox::AcceptRole);
      if ( ! label2.isEmpty())
        { button2 = msg.addButton(label2,QMessageBox::RejectRole);
          button2->setFocusPolicy(Qt::StrongFocus);
          msg.setTabOrder(button1,button2);
        }
      if ( ! label3.isEmpty())
        { button3 = msg.addButton(label3,QMessageBox::RejectRole);
          button3->setFocusPolicy(Qt::StrongFocus);
          msg.setTabOrder(button2,button3);
        }
    }
  button1->setFocusPolicy(Qt::StrongFocus);   //  Mac buttons are not tabable by default
  button1->setFocus(Qt::OtherFocusReason);    //    so ensure that it is the same on all gui's.

  msg.exec();

  if (msg.clickedButton() == button1)
    return (0);
  else if (msg.clickedButton() == button2)
    return (1);
  else
    return (2);
}
