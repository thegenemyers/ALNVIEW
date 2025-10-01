#include <stdio.h>
#include <math.h>

#include <QtGui>

extern "C" {
#include "GDB.h"
#include "align.h"
#include "alncode.h"
#include "ONElib.h"
#include "gene_core.h"
}

#include "open_window.h"
#include "main_window.h"


/*****************************************************************************\
*
*  MY LINE-EDIT
*
\*****************************************************************************/

MyLineEdit::MyLineEdit(QWidget *parent) : QLineEdit(parent) { process = false; }

void MyLineEdit::keyPressEvent(QKeyEvent *event)
{ if (process)
    { process = false;
      emit touched();
    }
  QLineEdit::keyPressEvent(event);
}

void MyLineEdit::mousePressEvent(QMouseEvent *event)
{ if (process)
    { process = false;
      emit touched();
    }
  QLineEdit::mousePressEvent(event);
}

void MyLineEdit::focusOutEvent(QFocusEvent *event)
{ emit focusOut();
  QLineEdit::focusOutEvent(event);
}

void MyLineEdit::processed(bool on)
{ process = on; }


/*****************************************************************************\
*
*  OPEN DIALOG
*
\*****************************************************************************/

void OpenDialog::openALN()
{ QString dir;

  if (alnInfo == NULL)
    dir = tr(".");
  else
    dir = alnInfo->absolutePath();
  QString name = QFileDialog::getOpenFileName(this,
                    tr("Open .1aln file"),dir,tr("Alignment file (*.1aln)"));
  if ( ! name.isNull())
    { delete alnInfo;
      alnInfo = new QFileInfo(name);
      alnText = alnInfo->fileName();
      alnFile->setText(alnText);
    }
}

void OpenDialog::activateLongest(int state)
{ bool on;

  on = (state == Qt::Checked);

  lLabel->setEnabled(on);
  lSuffix->setEnabled(on);
  lCutoff->setEnabled(on);
  if (on == false)
    lCutoff->setText(tr(""));
  else
    { if (lCut >= 0)
        lCutoff->setText(tr("%1").arg(lCut));
    }
}

void OpenDialog::activateIdentity(int state)
{ bool on;

  on = (state == Qt::Checked);

  iLabel->setEnabled(on);
  iCutoff->setEnabled(on);
  if (on == false)
    iCutoff->setText(tr(""));
  else
    { if (iCut >= 0)
        iCutoff->setText(tr("%1").arg(iCut));
    }
}

void OpenDialog::activateSize(int state)
{ bool on;

  on = (state == Qt::Checked);

  sLabel->setEnabled(on);
  sCutoff->setEnabled(on);
  if (on == false)
    sCutoff->setText(tr(""));
  else
    { if (sCut >= 0)
        sCutoff->setText(tr("%1").arg(sCut));
    }
}

void OpenDialog::aboutTo()
{ QString    dir;
  QFileInfo *ninfo;

  alnText = alnFile->text();
  if (alnText.isEmpty())
    { DotWindow::warning(
             tr("ALN-file name is emtpy"),
             this,DotWindow::WARNING,tr("OK"));
      return;
    }
  if (alnInfo == NULL)
    dir = tr(".");
  else
    dir = alnInfo->absolutePath();
  ninfo = new QFileInfo(dir + tr("/") + alnText);
  if ( ! ninfo->exists())
    { DotWindow::warning(
             tr("File ")+ninfo->absoluteFilePath()+tr(" does not exist"),
             this,DotWindow::WARNING,tr("OK"));
      delete ninfo;
      return;
    }
  delete alnInfo;
  alnInfo = ninfo;

  if (lBox->isChecked())
    { if (lCutoff->text().isEmpty())
        { DotWindow::warning(tr("Number of alignments not specified"),
                     this,DotWindow::ERROR,tr("OK"));
          return;
        }
      lCut = lCutoff->text().toInt();
      if (lCut <= 0)
        { DotWindow::warning(tr("Number of alignments (%1) is not a positive integer").arg(lCut),
                     this,DotWindow::ERROR,tr("OK"));
          return;
        }
    }
  else
    lCut = -1;

  if (iBox->isChecked())
    { if (iCutoff->text().isEmpty())
        { DotWindow::warning(tr("percent identity not specified"),
                     this,DotWindow::ERROR,tr("OK"));
          return;
        }
      iCut = iCutoff->text().toInt();
      if (iCut <= 0 || iCut > 100)
        { DotWindow::warning(tr("Percentage (%1) is not in [0,100]").arg(iCut),
                     this,DotWindow::ERROR,tr("OK"));
          return;
        }
    }
  else
    iCut = -1;

  if (sBox->isChecked())
    { if (sCutoff->text().isEmpty())
        { DotWindow::warning(tr("Alignment length not specified"),
                     this,DotWindow::ERROR,tr("OK"));
          return;
        }
      sCut = sCutoff->text().toInt();
      if (sCut <= 0)
        { DotWindow::warning(tr("Length (%1) is not a positive integer").arg(sCut),
                     this,DotWindow::ERROR,tr("OK"));
          return;
        }
    }
  else
    sCut = -1;

  accept();
}

OpenDialog::OpenDialog(QWidget *parent) : QDialog(parent)
{
  alnInfo = NULL;
  lCut = iCut = sCut = -1;

  QIntValidator *validInt = new QIntValidator(1,INT32_MAX,this);

  QLabel *alnLabel = new QLabel(tr("ALN file:"));
  QPushButton *alnSelect = new QPushButton("Pick");

  alnFile = new MyLineEdit();

  lLabel  = new QLabel(tr("Longest"));
  lSuffix = new QLabel(tr("ALNs"));
  lCutoff = new QLineEdit();
    lCutoff->setFixedWidth(80);
    lCutoff->setValidator(validInt);
  lBox = new QCheckBox();
    lBox->setCheckState(Qt::Unchecked);
    lLabel->setEnabled(false);
    lCutoff->setEnabled(false);
    lSuffix->setEnabled(false);

  iLabel  = new QLabel(tr("% id or better"));
  iCutoff = new QLineEdit();
    iCutoff->setFixedWidth(80);
    iCutoff->setValidator(validInt);
    iCutoff->setAlignment(Qt::AlignRight);
  iBox = new QCheckBox();
    iBox->setCheckState(Qt::Unchecked);
    iLabel->setEnabled(false);
    iCutoff->setEnabled(false);

  sLabel  = new QLabel(tr("bp or longer"));
  sCutoff = new QLineEdit();
    sCutoff->setFixedWidth(80);
    sCutoff->setValidator(validInt);
    sCutoff->setAlignment(Qt::AlignRight);
  sBox = new QCheckBox();
    sBox->setCheckState(Qt::Unchecked);
    sLabel->setEnabled(false);
    sCutoff->setEnabled(false);

  QHBoxLayout *file = new QHBoxLayout();
    file->addWidget(alnLabel);
    file->addWidget(alnFile,1);
    file->addWidget(alnSelect);

  QHBoxLayout *leng = new QHBoxLayout();
    leng->addSpacing(20);
    leng->addWidget(lBox);
    leng->addWidget(lLabel);
    leng->addWidget(lCutoff);
    leng->addWidget(lSuffix);
    leng->addStretch(1);

  QHBoxLayout *id = new QHBoxLayout();
    id->addSpacing(20);
    id->addWidget(iBox);
    id->addWidget(iCutoff);
    id->addWidget(iLabel);
    id->addStretch(1);

  QHBoxLayout *size = new QHBoxLayout();
    size->addSpacing(20);
    size->addWidget(sBox);
    size->addWidget(sCutoff);
    size->addWidget(sLabel);
    size->addStretch(1);

  QVBoxLayout *select = new QVBoxLayout();
    select->addLayout(file);
    select->addSpacing(5);
    select->addLayout(leng);
    select->addLayout(id);
    select->addLayout(size);

  cancel = new QPushButton("Cancel");
  open = new QPushButton("Open");

  QHBoxLayout *decision = new QHBoxLayout();
    decision->addStretch(1);
    decision->addWidget(cancel);
    decision->addSpacing(5);
    decision->addWidget(open);

  QVBoxLayout *central = new QVBoxLayout();
    central->addLayout(select);
    central->addSpacing(15);
    central->addStretch(1);
    central->addLayout(decision);

  setLayout(central);

  open->setDefault(true);
  setWindowTitle(tr("Open ALN Dataset"));
  setModal(true);
  setSizeGripEnabled(true);

  readAndApplySettings();

  connect(open,SIGNAL(clicked()),this,SLOT(aboutTo()));
  connect(cancel,SIGNAL(clicked()),this,SLOT(reject()));

  connect(alnSelect,SIGNAL(clicked()),this,SLOT(openALN()));

  connect(lBox,SIGNAL(stateChanged(int)),this,SLOT(activateLongest(int)));
  connect(iBox,SIGNAL(stateChanged(int)),this,SLOT(activateIdentity(int)));
  connect(sBox,SIGNAL(stateChanged(int)),this,SLOT(activateSize(int)));
}

void OpenDialog::getState(Open_State &state)
{ state.alnInfo       = alnInfo;
  state.alnText       = alnText;
  state.longOn        = lBox->isChecked();
  state.longCut       = lCut;
  state.idOn          = iBox->isChecked();
  state.idCut         = iCut;
  state.sizeOn        = sBox->isChecked();
  state.sizeCut       = sCut;
}

void OpenDialog::putState(Open_State &state)
{ QPixmap blob = QPixmap(16,16);

  alnInfo       = state.alnInfo;
  alnText       = state.alnText;

  lCut          = state.longCut;
  iCut          = state.idCut;
  sCut          = state.sizeCut;

  alnFile->setText(alnText);
  
  if (state.longOn)
    activateLongest(Qt::Checked);
  else
    activateLongest(Qt::Unchecked);
  lBox->setChecked(state.longOn);
  
  if (state.idOn)
    activateIdentity(Qt::Checked);
  else
    activateIdentity(Qt::Unchecked);
  iBox->setChecked(state.idOn);
  
  if (state.sizeOn)
    activateSize(Qt::Checked);
  else
    activateSize(Qt::Unchecked);
  sBox->setChecked(state.sizeOn);
}

void OpenDialog::readAndApplySettings()
{ Open_State state;

  QSettings settings("FASTGA", "ALNview");

  settings.beginGroup("open");
    QString filepath = settings.value("aln",tr("")).toString();
    state.alnText  = settings.value("aln",tr("")).toString();
    state.longCut  = settings.value("longCut",-1).toInt();
    state.longOn   = settings.value("longOn",false).toBool();
    state.idCut    = settings.value("idCut",-1).toInt();
    state.idOn     = settings.value("idOn",false).toBool();
    state.sizeCut  = settings.value("sizeCut",-1).toInt();
    state.sizeOn   = settings.value("sizeOn",false).toBool();
  settings.endGroup();

  if (filepath.isEmpty())
    { state.alnInfo = NULL;
      state.alnText = tr("");
    }
  else
    { state.alnInfo = new QFileInfo(filepath);
      state.alnText = state.alnInfo->fileName();
    }

  putState(state);
}

void OpenDialog::writeSettings(QSettings &settings)
{ Open_State state;

  getState(state);

  settings.beginGroup("open");
    if (state.alnInfo != NULL)
      settings.setValue("aln", state.alnInfo->absoluteFilePath());
    else
      settings.setValue("aln", tr(""));
    settings.setValue("longCut", state.longCut);
    settings.setValue(  "idCut", state.idCut);
    settings.setValue("sizeCut", state.sizeCut);
    settings.setValue("longOn", state.longOn);
    settings.setValue(  "idOn", state.idOn);
    settings.setValue("sizeOn", state.sizeOn);
  settings.endGroup();
}
