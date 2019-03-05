
# available values are 'armawrap', 'newmat', and 'armadillo'
LIB = armawrap

CFLAGS   = -O3
LDFLAGS  =
CFLAGS  += -Wall -pedantic -ansi
CFLAGS  += -Wno-unknown-warning-option
CFLAGS  += -Wno-misleading-indentation
CFLAGS  += -Wno-c++11-compat
CFLAGS  += -Wno-terminate
CFLAGS  += --coverage -fno-inline

CFLAGS += -DARMAWRAP_API=0
CFLAGS += -DARMADILLO_API=1
CFLAGS += -DNEWMAT_API=2

ifeq ($(LIB),armawrap)
CFLAGS += -DLIB=0
else ifeq ($(LIB),armadillo)
CFLAGS += -DLIB=1
else ifeq ($(LIB),newmat)
CFLAGS += -DLIB=2
endif

ifeq ($(LIB),newmat)
CFLAGS  += -Inewmat/
LDFLAGS += -Lnewmat/
LDFLAGS += -lnewmat
CFLAGS  += -DUSING_NEWMAT_API

else
LDFLAGS += -llapack
LDFLAGS += -lblas
CFLAGS  += -DARMA_USE_LAPACK
CFLAGS  += -DARMA_USE_BLAS
CFLAGS  += -I../armawrap/
endif

ifeq ($(LIB),armawrap)
CFLAGS  += -DUSING_NEWMAT_API
endif

testsrcs = $(wildcard test*.cpp)
testexcs = $(basename $(testsrcs))

default: $(testexcs)

%: %.cpp
	@echo "Compiling $< -> $@ ..."
	@g++ $(CFLAGS) -o $@ $< $(LDFLAGS)

clean:
	@rm -f $(testexcs) *.gcda *.gcno
