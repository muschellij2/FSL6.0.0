BUILDDIR=$$PWD/../build/mistview

include(../common/common.pri)

TEMPLATE = app
CONFIG += qt
QT += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

LIBS = -lvtkGUISupportQt$$VTKSUFFIX -lvtkRenderingContextOpenGL$$VTKSUFFIX -lvtkRenderingOpenGL$$VTKSUFFIX \
    -lvtkInteractionStyle$$VTKSUFFIX \
    -lvtkViewsContext2D$$VTKSUFFIX -lvtkViewsCore$$VTKSUFFIX \
	-lvtkCommonExecutionModel$$VTKSUFFIX -lvtkRenderingImage$$VTKSUFFIX \
	-lvtkCommonDataModel$$VTKSUFFIX -lvtkIOImage$$VTKSUFFIX -lvtkIOLegacy$$VTKSUFFIX -lvtkFiltersCore$$VTKSUFFIX \
	-lvtkInteractionWidgets$$VTKSUFFIX -lvtkFiltersGeneral$$VTKSUFFIX \
	-lvtkImagingCore$$VTKSUFFIX -lvtkRenderingAnnotation$$VTKSUFFIX \
	-lvtkRenderingFreeType$$VTKSUFFIX -lvtkChartsCore$$VTKSUFFIX \
    -lvtkRenderingCore$$VTKSUFFIX -lvtkFiltersExtraction$$VTKSUFFIX -lvtkFiltersSources$$VTKSUFFIX \ 
	-lvtkRenderingContext2D$$VTKSUFFIX -lvtkCommonColor$$VTKSUFFIX -lvtkCommonSystem$$VTKSUFFIX \
	-lvtkCommonMisc$$VTKSUFFIX -lvtkCommonTransforms$$VTKSUFFIX -lvtkCommonMath$$VTKSUFFIX \
	-lvtkCommonCore$$VTKSUFFIX -lvtksys$$VTKSUFFIX -lvtkzlib$$VTKSUFFIX \
	-lvtkfreetype$$VTKSUFFIX -lvtkpng$$VTKSUFFIX \
    -lvtkzlib$$VTKSUFFIX \
    $$LIBS

LIBS += -lGL -lXt -lX11

macx {
	LIBS += -framework Cocoa -framework IOKit -framework OpenGL
}

FORMS    += mainwindow.ui \
    setupfilenameswindow.ui \
    plotwindow.ui \
    multishapewindow.ui

SOURCES += main.cpp\
        mainwindow.cpp \
    orthoviews.cpp \
    viewbase.cpp \
    viewdata.cpp \
    viewshape.cpp \
    normalview.cpp \
    custominteractorstyle.cpp \
    setupfilenameswindow.cpp \
    plotwindow.cpp \
    multishapewindow.cpp \
    nokeysqvtkwidget.cpp

HEADERS  += mainwindow.h \
    orthoviews.h \
    viewbase.h \
    viewdata.h \
    viewshape.h \
    normalview.h \
    custominteractorstyle.h \
    setupfilenameswindow.h \
    plotwindow.h \
    multishapewindow.h \
    nokeysqvtkwidget.h

