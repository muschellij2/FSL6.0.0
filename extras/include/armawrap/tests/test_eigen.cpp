#include "tests.hpp"

int main(int argc, char *argv[]) {

  SymmetricMatrix  s(5);
  DiagonalMatrix vals;
  Matrix         vecs;

  randu(s);

  EigenValues(s, vals, vecs);

  printmat("s: ", s);
  cout << s << endl;
  printmat("vals: ", vals);
  cout << vals << endl;
  printmat("vecs: ", vecs);
  cout << vecs << endl;


  return 0;
}
