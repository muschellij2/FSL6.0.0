#!/bin/sh

export QTDIR=/opt/Qt/5.5/clang_64

$QTDIR/bin/qmake -spec mkspecs/macx-clang-libc++ CONFIG+=release

cd build
make
cd ..

