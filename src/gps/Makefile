include ${FSLCONFDIR}/default.mk

PROJNAME=gps

USRINCFLAGS=-I${INC_NEWMAT} -I${INC_BOOST} 
USRLDFLAGS=-L${LIB_NEWMAT} 

OBJS=gps.o
LIBS=-lmiscmaths -lNewNifti -lznz -lutils -lnewmat -lm -lz

XFILES=gps

all: ${XFILES}

gps: ${OBJS}
	${CXX} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}
