include ${FSLCONFDIR}/default.mk

PROJNAME = mcflirt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

MOBJS = mcflirt.o Globaloptions.o Log.o
M2OBJS = mcflirt2.o

LIBS = -lnewimage -lmiscmaths -lutils -lNewNifti -lznz -lnewmat -lprob -lm -lz

XFILES = mcflirt 

all:	${XFILES}

mcflirt:${MOBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${MOBJS} ${LIBS}

mcflirt2:${M2OBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${M2OBJS} ${LIBS}
