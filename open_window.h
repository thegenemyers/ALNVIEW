#ifndef OPEN_DIALOG_H
#define OPEN_DIALOG_H

#include <QtGui>
#include <QtWidgets>

extern "C" {
#include "GDB.h"
#include "ONElib.h"
#include "gene_core.h"
}

class DotWindow;

class MyLineEdit : public QLineEdit
{
  Q_OBJECT

public:
  MyLineEdit(QWidget *parent = 0);

  void processed(bool on);

signals:
  void touched();
  void focusOut();

protected:
  void keyPressEvent(QKeyEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void focusOutEvent(QFocusEvent *event);

private:
  bool process;
};

typedef struct
  { QFileInfo *alnInfo;
    QString    alnText;
    bool       longOn;
    int        longCut;
    bool       idOn;
    int        idCut;
    bool       sizeOn;
    int        sizeCut;
  } Open_State;

class OpenDialog : public QDialog
{
  Q_OBJECT

public:
  OpenDialog(QWidget *parent = 0);

  void getState(Open_State &state);
  void putState(Open_State &state);

  void readAndApplySettings();
  void writeSettings(QSettings &);

private slots:
  void activateLongest(int);
  void activateIdentity(int);
  void activateSize(int);

  void openALN();
  void aboutTo();

private:
  QCheckBox   *lBox;
    QLabel      *lLabel;
    QLineEdit   *lCutoff;
    QLabel      *lSuffix;

  QCheckBox   *iBox;
    QLineEdit   *iCutoff;
    QLabel      *iLabel;

  QCheckBox   *sBox;
    QLineEdit   *sCutoff;
    QLabel      *sLabel;

  QPushButton *open;
  QPushButton *cancel;

  QLineEdit *alnFile;
  QFileInfo *alnInfo;
  QString    alnText;

  int        lCut, iCut, sCut;
};

#endif
