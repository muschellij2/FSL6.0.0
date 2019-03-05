# $Id: systemvars.mk,v 1.13 2017/01/20 15:47:17 mwebster Exp $

# System dependent paths

RM = /bin/rm
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
CP = /bin/cp
MV = /bin/mv
INSTALL = install -p
TCLSH = ${FSLDIR}/bin/fsltclsh
RANLIB = echo

FSLML = ${FSLDIR}/bin/fslml

# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# Compiler dependent variables

CC = gcc
CXX = c++
CXX11 = scl enable devtoolset-2 -- c++
CSTATICFLAGS = -static
CXXSTATICFLAGS = -static

ARCHFLAGS = -m64 

PARALLELFLAGS = -fopenmp

DEPENDFLAGS = -MM

OPTFLAGS = -g -O3 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS = -g
GNU_ANSI_FLAGS = -Wall -ansi -pedantic -Wno-long-long
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}

# CUDA development environment
CUDA_INSTALLATION = /opt/cuda-7.5
GENCODE_FLAGS = $(shell ${FSLDIR}/config/common/supportedGencodes.sh ${CUDA_INSTALLATION})
LIB_CUDA = ${CUDA_INSTALLATION}/lib64
INC_CUDA = ${CUDA_INSTALLATION}/include
NVCC = ${CUDA_INSTALLATION}/bin/nvcc
NVCC11=scl enable devtoolset-2 -- ${CUDA_INSTALLATION}/bin/nvcc
