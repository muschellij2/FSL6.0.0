BOOSTSUFFIX =
NLOPTDIR = $$(HOME)/src/nlopt-install
SQLITEDIR = $$(HOME)/src/sqlite-install
VTKDIR = $$(HOME)/src/VTK-7.0.0-static-install
VTKSUFFIX = -7.0
FSLDIR = /opt/fmrib/fsl

INCLUDEPATH += $$NLOPTDIR/include $$SQLITEDIR/include $$VTKDIR/include/vtk$$VTKSUFFIX
INCLUDEPATH += $$FSLDIR/include $$FSLDIR/include/newimage $$FSLDIR/extras/include $$FSLDIR/extras/include/newmat $$FSLDIR/extras/include/libgdc $$FSLDIR/extras/include/boost
LIBS += -L$$NLOPTDIR/lib -L$$VTKDIR/lib
LIBS += -L$$FSLDIR/lib -L$$FSLDIR/extras/lib

