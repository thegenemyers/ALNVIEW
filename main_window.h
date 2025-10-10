#ifndef IMAGE_WINDOW_H
#define IMAGE_WINDOW_H

#include <QtGui>
#include <QtWidgets>
#include "open_window.h"

extern "C" {
#include "GDB.h"
#include "sticks.h"
}

class DotCanvas;
class DotWindow;
class OpenDialog;

/***************************************************************************************/
/*                                                                                     */
/*   VIEW CONTROL STATE & DOT PLOT MODEL                                               */
/*                                                                                     */
/***************************************************************************************/

typedef enum { TOP_LEFT = 0, TOP_RIGHT = 1, BOTTOM_LEFT = 2, BOTTOM_RIGHT = 3 } LocatorQuad;

typedef struct
{ QRect           wGeom;

  View            view;
  double          zoom;
  int             format;

  bool            fOn;
  Focus           focus;
  QColor          fColor;
  bool            fViz;

  int             nlays;
  int             order[MAX_LAYERS];
  bool            on[MAX_LAYERS];
  QColor          colorF[MAX_LAYERS];
  QColor          colorR[MAX_LAYERS];
  int             thick[MAX_LAYERS];
  
  QColor          lColor;
  bool            lViz;
  LocatorQuad     lQuad;

  double          lMag;
  double          lXct;
  double          lYct;

  QStack<double>  zMag;
  QStack<double>  zXct;
  QStack<double>  zYct;

} DotState;

class DotLineEdit : public QLineEdit
{
  Q_OBJECT

public:
  DotLineEdit(QWidget *parent = 0);
  void setChain(QWidget *p, QWidget *s);

signals:
  void focusOut();

protected:
  void focusOutEvent(QFocusEvent *event);
  void keyPressEvent(QKeyEvent *event);
  bool focusNextPrevChild(bool next);

private:
  QWidget *mypred;
  QWidget *mysucc;
};

class SeqLineEdit : public QLineEdit
{
  Q_OBJECT

public:
  SeqLineEdit(QWidget *parent = 0, int id = 0);

signals:
  void focusShift(int, int);

protected:
  void keyPressEvent(QKeyEvent *);
  void mousePressEvent(QMouseEvent *);
  void mouseMoveEvent(QMouseEvent *);
  void mouseReleaseEvent(QMouseEvent *);
  void enterEvent(QEnterEvent *);
  void leaveEvent(QEvent *);

private:
  int  myid;
};

class MyMenu : public QMenu
{
  Q_OBJECT

public:
  MyMenu(QWidget *parent = 0);

protected:
  void mouseReleaseEvent(QMouseEvent *ev);
  void mousePressEvent(QMouseEvent *ev);
  void keyReleaseEvent(QKeyEvent *ev);
};


/***************************************************************************************/
/*                                                                                     */
/*   DOT CANVAS                                                                        */
/*                                                                                     */
/***************************************************************************************/

class DotCanvas : public QWidget
{
  Q_OBJECT

public:
  DotCanvas(QWidget *parent = 0);
  ~DotCanvas();

  Frame *shareData(DotState *state, DotPlot *plot);
  bool   zoomView(double zoomDel);
  void   resetView();
  bool   viewToFrame();

  static int labelWidth;

  static QVector<QRgb> ctable;
  static QVector<uchar> imbit;

public slots:
  void zoomPop();

signals:
  void NewFrame(double newZ);
  void NewFocus();
  void Moved();

protected:
  void resizeEvent(QResizeEvent *event);

  void mousePressEvent(QMouseEvent *);
  void mouseMoveEvent(QMouseEvent *);
  void mouseReleaseEvent(QMouseEvent *);

  void keyPressEvent(QKeyEvent *);
  void keyReleaseEvent(QKeyEvent *);

  void paintEvent(QPaintEvent *event);

  void enterEvent(QEnterEvent *);
  void leaveEvent(QEvent *);

private slots:
  void showAlign();

private:
  DotSegment *pick(int x, int y, DotLayer **layer);
  DotSegment *pickedSeg;
  DotLayer   *pickedLayer;

  bool         noFrame;
  Frame        frame;

  int          rectW;
  int          rectH;

  QImage      *image;
  uchar      **raster;

  DotPlot     *plot;
  DotState    *state;

  int          mouseX;
  int          mouseY;
  double       scaleX;
  double       scaleY;

  int          lastKey;
  bool         select;
  bool         nograb;
  bool         picking;
  int          menuLock;

  QBasicTimer *timer;
  int          xpos;
  int          ypos;
  int          vwid;
  int          vhgt;

  QRubberBand *rubber;

  MyMenu  *popup;
  QAction *aline;
  QAction *bline;
  QAction *alignAct;
  MyMenu  *faraway;
};

class LayerWidget : public QWidget
{
  Q_OBJECT

public:
  LayerWidget(DotState *istate, DotCanvas *icanvas, QWidget *parent = 0);

protected:
  void mousePressEvent(QMouseEvent *ev);
  void dragEnterEvent(QDragEnterEvent *ev);
  void dragMoveEvent(QDragMoveEvent *ev);
  void dropEvent(QDropEvent *ev);

private:
  int        dropY;
  DotState  *state;
  DotCanvas *canvas;
};

/***************************************************************************************/
/*                                                                                     */
/*   ALIGNMENT VIEWER                                                                  */
/*                                                                                     */
/***************************************************************************************/

class AlignWindow : public QMainWindow
{
  Q_OBJECT

public:
  AlignWindow(QString title, char *atext, QWidget *parent = NULL, Qt::WindowFlags flags = Qt::WindowFlags());

private:
  QTextEdit *edit;
};


/***************************************************************************************/
/*                                                                                     */
/*   DOT WINDOW DATA TYPE                                                              */
/*                                                                                     */
/***************************************************************************************/

class DotWindow : public QMainWindow
{
  Q_OBJECT

public:
  DotWindow(DotPlot *plot, DotState *state, bool isCopy);
  ~DotWindow();

  typedef enum { INFORM, WARNING, ERROR } MessageKind;

  static int warning(const QString& message, QWidget *parent, MessageKind,
                     const QString& label1 = QString(),
                     const QString& label2 = QString(),
                     const QString& label3 = QString() );

  static OpenDialog *openDialog;
  static Open_State  dataset;

  static QRect *screenGeometry;
  static int    windowWidth;
  static int    windowHeight;

  static QList<DotWindow *> dotwindows;

public slots:
  static void openFile();

  void zoomUp();
  void zoomDown();
  void zoomTo();

  void viewChange();

  void formatChange(int index);
  void thickChange(int index);
 
  void focusOnChange();
  void focusChange();
  void focusColorChange();
  void hairsChange();

  void locatorChange();
  void locatorColorChange();

  void openOverlay();
  void openCopy();
  void closeAll();

  void toggleToolBar();
  void tileImages();
  void cascadeImages();
  void unminimizeAll();
  void raiseAll();

  void frameToView(double newZ);
  void clickToFocus();
  void seqMove(int,int);

   void addTextBox(AlignWindow *align);

protected:
  void closeEvent(QCloseEvent *event);

private slots:
  void checkRangeA();
  void checkRangeB();
  void checkFocus();
  void activateLayer(int);
  void layerFChange();
  void layerRChange();

private:
  void readAndApplySettings();
  void writeSettings();

  void pushState();

  DotPlot            *plot;
  Frame              *frame;
  QImage             *pixels;
  DotCanvas          *canvas;

  void                createMenus();
  void                createActions();
  void                createToolBars();

  QAction *exitAct;

  QAction *openAct;
  QAction *overAct;
  QAction *copyAct;
  QAction *toolAct;

  QAction *tileAct;
  QAction *cascadeAct;

  QAction *unminimizeAllAct;
  QAction *raiseAllAct;

  QToolBar        *fileToolBar;
  Qt::ToolBarArea  toolArea;

  DotState            state;

  DotLineEdit        *zoomEdit;

  DotLineEdit        *Arng;
  DotLineEdit        *Brng;
  DotLineEdit        *Fpnt;

  QCheckBox          *focusOn;
  QToolButton        *focusBox;
  QCheckBox          *focusCheck;

  QComboBox          *cFormat;

  LayerWidget *layerPanel;
  int          nmasks;
  QWidget   *layerWidget[MAX_LAYERS];
    QCheckBox   *layerOn[MAX_LAYERS];
    QLabel      *layerTitle[MAX_LAYERS];
    QToolButton *layerFBox[MAX_LAYERS];
    QLabel      *layerFText[MAX_LAYERS];
    QToolButton *layerRBox[MAX_LAYERS];
    QLabel      *layerRText[MAX_LAYERS];
    QComboBox   *layerThick[MAX_LAYERS];

  QToolButton        *locatorBox;
  QCheckBox          *locatorCheck;
  QButtonGroup       *locatorQuad;

  QList<AlignWindow *> alignWindows;
};

#endif
