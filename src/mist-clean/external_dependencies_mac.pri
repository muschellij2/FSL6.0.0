# NOTE: Use -DCMAKE_INSTALL_NAME_DIR=<libdir> when compiling VTK on OSX

BOOSTSUFFIX = -mt
NLOPTDIR = $$(HOME)/src/nlopt-install
SQLITEDIR = $$(HOME)/src/sqlite-install
VTKDIR = /opt
VTKSUFFIX = -6.3
FSLDIR = /opt/fsl

INCLUDEPATH += $$NLOPTDIR/include $$SQLITEDIR/include $$VTKDIR/include/vtk$$VTKSUFFIX
INCLUDEPATH += $$FSLDIR/include $$FSLDIR/include/newimage $$FSLDIR/extras/include $$FSLDIR/extras/include/newmat $$FSLDIR/extras/include/libgdc
LIBS += -L$$NLOPTDIR/lib -L$$VTKDIR/lib
LIBS += -L$$FSLDIR/lib -L$$FSLDIR/extras/lib

