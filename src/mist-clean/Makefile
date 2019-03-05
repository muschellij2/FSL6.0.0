include ${FSLCONFDIR}/default.mk

NLOPT_INC = ${FSLEXTINC}
NLOPT_LIB = ${FSLEXTLIB}
SQLITE_INC = ${FSLEXTINC}/libsqlite
SQLITE_LIB = ${FSLEXTLIB}

PROJNAME = mist

XFILES = mist/mist
SCRIPTS = bin/mist_1_train bin/mist_2_fit bin/mist_FA_reg bin/mist_display bin/mist_mesh_utils

USRCXXFLAGS = -std=c++11
USRINCFLAGS = -I${FSLDIR}/include/newimage -I${INC_NEWMAT} -I${INC_ZLIB} -I${INC_GDC} -I${INC_GD} -I${SQLITE_INC} -I${NLOPT_INC} -I${VTKDIR_INC} -Icommon
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB} -L${LIB_BOOST} -L${LIB_GDC} -L${LIB_GD} -L${NLOPT_LIB} -L${VTKDIR_LIB}

UNAME := $(shell uname)
ifneq (${UNAME},Darwin)
        LIBRT = -lrt
endif

LDLIBS =	-lvtkIOLegacy${VTKSUFFIX} -lvtkIOCore${VTKSUFFIX} -lvtkFiltersModeling${VTKSUFFIX} -lvtkFiltersCore${VTKSUFFIX} \
			-lvtkCommonExecutionModel${VTKSUFFIX} -lvtkCommonDataModel${VTKSUFFIX} -lvtkCommonMisc${VTKSUFFIX} \
			-lvtkCommonSystem${VTKSUFFIX} -lvtkCommonTransforms${VTKSUFFIX} -lvtkCommonMath${VTKSUFFIX} \
			-lvtkCommonCore${VTKSUFFIX} -lvtksys${VTKSUFFIX} \
			-lboost_log -lboost_log_setup -lboost_thread -lboost_filesystem -lboost_date_time -lboost_chrono \
			-lboost_system -lboost_serialization -lboost_regex -lwarpfns -lbasisfield -lnewimage -lmiscmaths \
			-lnewmat -lgdc -lgd -lpng -lNewNifti -lutils -lznz -lm -lz -lnlopt -lpthread -ldl \
			${SQLITE_LIB}/libsqlite3.a ${LIBRT}

HDRS = common/gibbsshapemodel.h common/mrfshapemodel.h common/mvnshapemodel.h common/plotting.h common/profilefilters.h common/profilemixtures.h common/profilemodel.h common/profilepriors.h common/serialisation.h common/shape.h common/shapemodel.h common/stats.h common/transformation.h mist/builddate.h

OBJS = common/gibbsshapemodel.o common/mrfshapemodel.o common/mvnshapemodel.o common/plotting.o common/profilefilters.o common/profilemixtures.o common/profilemodel.o common/profilepriors.o common/serialisation.o common/shape.o common/shapemodel.o common/stats.o common/transformation.o

all:		${XFILES} installdata installpython

mist/mist:	${OBJS} mist/mist.o
			${CXX11} ${CXXFLAGS} ${LDFLAGS} -o $@ $^ ${LDLIBS}

# Always rebuild if something changed (to get the build date right)
mist/mist.o:	${OBJS}

%.o:		%.cpp ${HDRS}
			${CXX11} -c ${CXXFLAGS} -o $@ $<

installdata:
			mkdir -p ${DESTDIR}/data/mist
			cp -r data/masks data/meshes ${DESTDIR}/data/mist/

installpython:
		mkdir -p ${DESTDIR}/python/mist
		cp -r python/* ${DESTDIR}/python/mist

clean:
			rm -f ${OBJS} mist/mist.o mist/mist

.PHONY:		all clean installdata
