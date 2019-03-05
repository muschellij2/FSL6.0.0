# $Id: externallibs.mk,v 1.10 2018/08/31 14:54:40 mwebster Exp $

# External Library and Include Paths

FSLEXTLIB=${FSLDIR}/extras/lib
FSLEXTINC=${FSLDIR}/extras/include
FSLEXTBIN=${FSLDIR}/extras/bin

# GD library
LIB_GD = ${FSLEXTLIB}
INC_GD = ${FSLEXTINC}

# GDC library
LIB_GDC = ${FSLEXTLIB}
INC_GDC = ${FSLEXTINC}/libgdc

# LIBXML2 library
INC_XML2 = ${FSLEXTINC}/libxml2

# LIBXML++ library
INC_XML++ = ${FSLEXTINC}/libxml++-2.6
INC_XML++CONF = ${FSLEXTLIB}/libxml++-2.6/include
# GSL library
LIB_GSL = ${FSLEXTLIB}
INC_GSL = ${FSLEXTINC}/gsl

# PNG library
LIB_PNG = ${FSLEXTLIB}
INC_PNG = ${FSLEXTINC}

# PROB library
LIB_PROB = ${FSLEXTLIB}
INC_PROB = ${FSLEXTINC}/libprob

# CPROB library
LIB_CPROB = ${FSLEXTLIB}
INC_CPROB = ${FSLEXTINC}/libcprob

# NEWMAT library
#LIB_NEWMAT = ${FSLEXTLIB} -llapack -lblas or just -lopenblas
LIB_NEWMAT = ${FSLEXTLIB} -lopenblas
INC_NEWMAT = ${FSLDIR}/extras/include/armawrap/armawrap -DARMA_USE_LAPACK -DARMA_USE_BLAS -DARMA_64BIT_WORD

# NEWRAN library
LIB_NEWRAN = ${FSLEXTLIB}
INC_NEWRAN = ${FSLEXTINC}/newran

# ZLIB library
LIB_ZLIB = /lib64
INC_ZLIB = /usr/include

# BOOST library
BOOSTDIR = ${FSLEXTINC}/boost
LIB_BOOST = ${BOOSTDIR}
INC_BOOST = ${BOOSTDIR}

# QT library
QTDIR = /usr/lib/qt3
LIB_QT = ${QTDIR}/lib
INC_QT = ${QTDIR}/include

# QWT library
QWTDIR = /usr/local/qwt
LIB_QWT = ${QWTDIR}/lib
INC_QWT = ${QWTDIR}/include
 
# FFTW3 library
LIB_FFTW3 = ${FSLEXTLIB}
INC_FFTW3 = ${FSLEXTINC}/fftw3

# VTK library 
VTKDIR_INC = /home/fs0/cowboy/var/caper_linux_64-gcc4.4/VTK7/include/vtk-7.0
VTKDIR_LIB = /home/fs0/cowboy/var/caper_linux_64-gcc4.4/VTK7/lib
VTKSUFFIX = -7.0
