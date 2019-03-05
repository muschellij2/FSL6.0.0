//#define ARMA_EXTRA_DEBUG
#include "tests.hpp"

int main(int argc, char *argv[]) {

  Matrix m(8, 8);
  randu(m);
  Matrix b = m.Rows(3, 5).SubMatrix(1, 2, 3, 6);
  Matrix c = m.SubMatrix(1, 6, 3, 6).Rows(4, 6);
  Matrix d = (m.SymSubMatrix(3, 7) * 100).Rows(1, 3) / 50;

  printmat("m: ", m);
  printmat("b: ", b);
  printmat("c: ", c);
  printmat("d: ", d);
 
  
  return 0;
}
