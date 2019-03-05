QMAKE_CXXFLAGS += -std=c++11
# -DBOOST_ALL_DYN_LINK

# This is a synonym for -Wextra, which is enabled by qmake if warnings are enabled
QMAKE_CXXFLAGS_WARN_ON -= -W

# GCC produces these at -Wall, which makes the important warnings hard to spot
QMAKE_CXXFLAGS_WARN_ON += -Wno-sign-compare -Wno-deprecated-declarations

QMAKE_MAC_SDK = macosx10.11

DESTDIR = $$BUILDDIR
OBJECTS_DIR = $$BUILDDIR
MOC_DIR = $$BUILDDIR
RCC_DIR = $$BUILDDIR
UI_DIR = $$BUILDDIR
MAKEFILE = $$BUILDDIR/Makefile

include(../external_dependencies.pri)

INCLUDEPATH += $$PWD


LIBS += -lwarpfns -lbasisfield -lnewimage -lmiscmaths -lnewmat -lgdc -lgd -lpng12 -lfslio -lniftiio -lutils -lznz -lm -lz
LIBS += -lnlopt
LIBS += -lvtkIOLegacy$$VTKSUFFIX -lvtkIOCore$$VTKSUFFIX -lvtkFiltersModeling$$VTKSUFFIX -lvtkFiltersCore$$VTKSUFFIX \
	-lvtkCommonExecutionModel$$VTKSUFFIX -lvtkCommonDataModel$$VTKSUFFIX -lvtkCommonMisc$$VTKSUFFIX \
	-lvtkCommonSystem$$VTKSUFFIX -lvtkCommonTransforms$$VTKSUFFIX -lvtkCommonMath$$VTKSUFFIX \
	-lvtkCommonCore$$VTKSUFFIX -lvtksys$$VTKSUFFIX
LIBS += -lboost_log$$BOOSTSUFFIX -lboost_log_setup$$BOOSTSUFFIX -lboost_thread$$BOOSTSUFFIX \
	-lboost_filesystem$$BOOSTSUFFIX -lboost_date_time$$BOOSTSUFFIX -lboost_chrono$$BOOSTSUFFIX \
	-lboost_system$$BOOSTSUFFIX -lboost_serialization$$BOOSTSUFFIX -lboost_regex$$BOOSTSUFFIX
LIBS += -lrt -lpthread -ldl $$SQLITEDIR/lib/libsqlite3.a

# NOTE: This is not compiled as a library because boost::serialization seems to require that all serialisation code
# is moved out of the public header files in that case ...

SOURCES += $$PWD/shape.cpp \
    $$PWD/shapemodel.cpp \
    $$PWD/mvnshapemodel.cpp \
    $$PWD/profilemodel.cpp \
    $$PWD/profilepriors.cpp \
    $$PWD/serialisation.cpp \
    $$PWD/transformation.cpp \
    $$PWD/profilefilters.cpp \
    $$PWD/stats.cpp \
    $$PWD/profilemixtures.cpp \
    $$PWD/plotting.cpp \
    $$PWD/mrfshapemodel.cpp \
    $$PWD/gibbsshapemodel.cpp

HEADERS += \
    $$PWD/stats.h \
    $$PWD/profilemodel.h \
    $$PWD/profilemixtures.h \
    $$PWD/shapemodel.h \
    $$PWD/mvnshapemodel.h \
    $$PWD/serialisation.h \
    $$PWD/shape.h \
    $$PWD/profilepriors.h \
    $$PWD/transformation.h \
    $$PWD/profilefilters.h \
    $$PWD/plotting.h \
    $$PWD/mrfshapemodel.h \
    $$PWD/gibbsshapemodel.h
