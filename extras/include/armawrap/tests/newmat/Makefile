all: libnewmat.a

clean:
	@rm -f *.o *.a

%.o: %.cpp
	@g++  -c $*.cpp

OBJ_LIB = bandmat.o cholesky.o evalue.o fft.o hholder.o     \
  jacobi.o myexcept.o newmat1.o newmat2.o newmat3.o         \
  newmat4.o newmat5.o newmat6.o newmat7.o newmat8.o         \
  newmat9.o newmatex.o newmatnl.o newmatrm.o solution.o     \
  sort.o submat.o svd.o newfft.o

libnewmat.a: $(OBJ_LIB)
	@ar cr $@ $(OBJ_LIB)
	@ranlib $@
