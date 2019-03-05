include ${FSLCONFDIR}/default.mk

PROJNAME = flame

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB} -DCIFTILIB_USE_XMLPP -I${FSLEXTINC} -I${INC_XML2} -I${INC_XML++} -I${INC_XML++CONF} -I${INC_BOOST} -I${FSLDIR}/include/ciftiio
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

UNAME := $(shell uname)
ifeq (${UNAME},Darwin)
        LIBS = -liconv
endif

LIBS += -lutils -lnewimage -lmiscmaths -lm -lnewmat -lNewNifti -lcifti -lxml++-2.6 -lxml2 -lboost_filesystem -lboost_system -lznz -lprob -lz 

XFILES = flameo

OBJS = design.o mcmc_mh.o gsoptions.o gsmanager.o gsa.o

all:	${XFILES}

flameo:${OBJS} flameo.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} flameo.o ${LIBS}
