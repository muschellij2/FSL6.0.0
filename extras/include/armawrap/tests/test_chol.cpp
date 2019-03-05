#include "tests.hpp"

int main(int argc, char *argv[]) {

  DiagonalMatrix A(5);

  randu(A);

  SymmetricMatrix       PD = A;
  LowerTriangularMatrix L  = Cholesky(PD);

  cout << "PD: " << endl;
  cout <<  PD    << endl;
  cout << "L: "  << endl;
  cout <<  L     << endl;

  return 0;
}
