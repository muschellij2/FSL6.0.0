include ${FSLCONFDIR}/default.mk

PROJNAME = NewNifti

USRINCFLAGS = -I${INC_ZLIB} -DHAVE_ZLIB

LIBS =-lm -lznz -lz

OBJS=legacyFunctions.o NewNifti.o

XFILES=

all: libNewNifti.a ${XFILES}

libNewNifti.a: ${OBJS}
	${AR} -r libNewNifti.a ${OBJS}

test: test.o ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ test.o ${OBJS} ${LIBS}
