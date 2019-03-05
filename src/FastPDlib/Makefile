include ${FSLCONFDIR}/default.mk

PROJNAME=FastPD

USRCXXFLAGS=-g 

USRINCFLAGS =-I${INC_NEWMAT} -I${INC_NEWRAN} -I${INC_BOOST}  -I${HOME}/fsldev/include 
USRLDFLAGS =-L${LIB_NEWMAT} -L${LIB_NEWRAN} -L${HOME}/fsldev/lib 

LIBS=-ldiscrete

ifdef HOCR
	USRCXXFLAGS:= -std=c++11
endif
#------------------------------------------------------------------------------		
MYPROGRAMOBJ=FastPD.o graph.o
MYPROGRAM=FastPD
MYlib=libfastPD	


#------------------------------------------------------------------------------

all: ${MYlib}.a

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

${MYlib}.a: ${MYPROGRAMOBJ}
	${AR} -r libfastpd.a ${MYPROGRAMOBJ}

