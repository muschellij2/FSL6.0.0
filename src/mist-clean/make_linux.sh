#!/bin/sh

export QMAKESPEC=/usr/lib64/qt4/mkspecs/linux-g++-64
/usr/bin/qmake-qt4 CONFIG+=release

cd build
make
cd ..

