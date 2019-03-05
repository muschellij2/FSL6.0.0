# $id: Makefile,v 1.1.1.1 2013/10/08 14:06:18 emmar Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = MSM
compilewithhocr=false
USRCXXFLAGS=-g
DISCRETELIB=-ldiscrete
TBBLIBS=
ifdef HOCR
	USRCXXFLAGS:= ${USRCXXFLAGS} -std=c++11  -DHAS_HOCR
	DISCRETELIB=-ldiscreteHOCR 
endif
ifdef TBB
	USRCXXFLAGS:= ${USRCXXFLAGS} -DHAS_TBB
	TBBLIBS= -ltbb -lANN
endif 

#Matlab version not for release
#USRINCFLAGS =-I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_ZLIB} -I${INC_CPROB} -I${INC_PROB}  -I${INC_GD}  -I${INC_GDC} -I${INC_PNG} -I${INC_BOOSTNEW} -I./include -I${HOME}/fsldev/include -I/usr/lib/matlab/R2015a/extern/include/  -I/vol/vipdata/packages/tbb41_20130116oss/include
#USRLDFLAGS =-Wl,-rpath  /usr/lib/matlab/R2015a/bin/glnxa64/  -L/usr/lib/matlab/R2015a/bin/glnxa64/ -L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_ZLIB}  -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_GD}   -L. -L${HOME}/fsldev/lib -L/vol/vipdata/packages/tbb41_20130116oss/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/


USRINCFLAGS =-I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_ZLIB} -I${INC_CPROB} -I${INC_PROB}  -I${INC_GD}  -I${INC_GDC} -I${INC_PNG} -I${INC_BOOST} -I./include -I/vol/vipdata/packages/tbb41_20130116oss/include
USRLDFLAGS =-L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${LIB_ZLIB}  -L${LIB_CPROB} -L${LIB_PROB} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_GD} -L. -L/vol/vipdata/packages/tbb41_20130116oss/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21/
LIBS= -lmeshreg -lfastpd ${DISCRETELIB} ${TBBLIBS} -lnewmeshclass -lfslsurface -lgiftiio -lmm  -lexpat  -lfslvtkio -lnewimage  -lmiscmaths -lNewNifti -lmiscplot -lmiscpic -lnewmat -lutils  -lbint -lnewran -lprob -lznz -lgdc -lgd -lpng  -lz

XFILES=msm_metricmath msmresample  msmapplywarp estimate_strains estimate_curvature average_surfaces msm estimate_metric_distortion surfconvert \
calculatevertexwisecorr msmresamplemetric #simulatetask groupsegmentation #  Msmapplywarp msmresamplemetric estimate_metric_distortion average_surfaces

ifdef GROUP
        XFILES:= msm_group ${XFILES}           
endif


all: ${XFILES}

libmsm.a: ${MSMLIBOBJS}
	${AR} -r libmsm.a ${MSMLIBOBJS}


msm: msm.o
	${CXX} ${LDFLAGS} -o $@ msm.o ${LIBS} 
msm_group: msm_group.o
	${CXX} ${LDFLAGS} -o $@ msm_group.o ${LIBS} 
msm_makeicos: msm_makeicos.o
	${CXX} ${LDFLAGS} -o $@ msm_makeicos.o ${LIBS} 
intensity_normalise: intensity_normalise.o
	${CXX} ${LDFLAGS} -o $@ intensity_normalise.o ${LIBS} 
msm_metricmath: msm_metricmath.o
	${CXX} ${LDFLAGS} -o $@ msm_metricmath.o ${LIBS} 
msm_surface_average: msm_surface_average.o
	${CXX} ${LDFLAGS} -o $@ msm_surface_average.o ${LIBS} 
transform: transform.o
	${CXX} ${LDFLAGS} -o $@ transform.o ${LIBS} 
createartefacts: createartefacts.o
	${CXX} ${LDFLAGS} -o $@ createartefacts.o ${LIBS} 
calculatevertexwisecorr: calculatevertexwisecorr.o
	${CXX} ${LDFLAGS} -o $@ calculatevertexwisecorr.o ${LIBS} 
estimate_strains: estimate_strains.o
	${CXX} ${LDFLAGS} -o $@ estimate_strains.o ${LIBS} 
estimate_curvature: estimate_curvature.o
	${CXX} ${LDFLAGS} -o $@ estimate_curvature.o ${LIBS} 

###SparseFeaturesTest: SparseFeaturesTest.o
	#${CXX} ${LDFLAGS} -o $@ SparseFeaturesTest.o ${LIBS} 

engdemo: engdemo.o
	${CXX} ${LDFLAGS} -o $@ engdemo.o ${LIBS} 

simulatetask: simulatetask.o
	${CXX} ${LDFLAGS} -o $@ simulatetask.o ${LIBS} 

groupsegmentation: groupsegmentation.o
	${CXX} ${LDFLAGS} -o $@ groupsegmentation.o ${LIBS} 

msmapplywarp: msmapplywarp.o
	${CXX}  ${LDFLAGS} -o $@  msmapplywarp.o ${LIBS} 

msmresample: msmresample.o
	${CXX}  ${LDFLAGS} -o $@  msmresample.o ${LIBS} 

msmresamplemetric: msmresamplemetric.o
	${CXX} ${LDFLAGS} -o $@ msmresamplemetric.o ${LIBS} 

estimate_metric_distortion: estimate_metric_distortion.o
	${CXX} ${LDFLAGS} -o $@ estimate_metric_distortion.o ${LIBS} 

minimise_anat-to-sphere_distortions: minimise_anat-to-sphere_distortions.o
	${CXX} ${LDFLAGS} -o $@ minimise_anat-to-sphere_distortions.o ${LIBS} 

average_surfaces: average_surfaces.o
	${CXX}  ${LDFLAGS} -o $@  average_surfaces.o ${LIBS} 

average_surfaces_TBB: average_surfaces_TBB.o
	${CXX}  ${LDFLAGS} -o $@  average_surfaces_TBB.o ${LIBS} 

surfconvert: surfconvert.o
	${CXX}  ${LDFLAGS} -o $@  surfconvert.o ${LIBS} 


