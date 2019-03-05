# $Id: Makefile,v 1.31 2015/07/01 16:05:31 mwebster Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = first

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} -I${INC_PROB}

USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB} -L${LIB_PROB}



LIBS= -lshapeModel -lfirst_lib -lfslvtkio -lmeshclass -lnewimage -lprob -lmiscmaths -lNewNifti -lznz -lnewmat  -lutils -lz 

FMBLIBS= -lnewimage -lmiscmaths -lNewNifti -lznz -lnewmat  -lutils -lz 

XFILES= first first_utils first_mult_bcorr

SCRIPTS = run_first run_first_all first_flirt concat_bvars first_roi_slicesdir first_boundary_corr surface_fdr first3Dview

all: ${XFILES}

first: first.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ first.o ${LIBS}

first_utils: first_utils.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ first_utils.o ${LIBS}

first_mult_bcorr: first_mult_bcorr.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ first_mult_bcorr.o ${FMBLIBS}
