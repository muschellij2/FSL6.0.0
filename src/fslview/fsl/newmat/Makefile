include ${FSLCONFDIR}/default.mk

PROJNAME = dummy_newmat

all:
	g++ -c -o newmat.o newmat.cpp
	ar -cr libnewmat.a newmat.o
	ranlib libnewmat.a
