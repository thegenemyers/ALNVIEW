#include <unistd.h>
#include <sys/types.h>

#include <QtGui>

#ifdef Q_WS_MAC
#include <QMacStyle>
#endif

#include "main_window.h"

int main(int argc, char *argv[])
{
  QApplication app(argc, argv);

  DotWindow::openDialog = new OpenDialog(NULL);

  DotWindow::openFile();

  return app.exec();
}
