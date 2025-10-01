QMAKE_SPEC = macx-clang
CONFIG += release
QMAKE_CXXFLAGS += -DINTERACTIVE -fvisibility=hidden -Wall -Wno-unused-result
QMAKE_CFLAGS += -DINTERACTIVE -Wall -Wno-unused-result
QMAKE_LIBS += -lz

MOC_DIR = BUILD
OBJECTS_DIR = BUILD
RCC_DIR = BUILD

QT += widgets

HEADERS       = main_window.h open_window.h sticks.h doter.h alncode.h align.h gene_core.h ONElib.h GDB.h hash.h select.h
SOURCES       = main.cpp main_window.cpp open_window.cpp sticks.c doter.c alncode.c align.c gene_core.c ONElib.c GDB.c hash.c select.c
TARGET        = ALNview
RESOURCES     = viewer.qrc
