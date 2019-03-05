BUILDDIR = $$PWD/../build/mist

include(../common/common.pri)

TEMPLATE = app
CONFIG -= app_bundle qt

HEADERS += ../builddate.h

SOURCES += mist.cpp

builddate.commands = touch builddate.h
QMAKE_EXTRA_TARGETS = builddate
PRE_TARGETDEPS = builddate

