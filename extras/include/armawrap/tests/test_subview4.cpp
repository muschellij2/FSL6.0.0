#include "tests.hpp"

int main(int argc, char *argv[]) {

  ColumnVector vrow( 6);
  Matrix       acEst(6, 6);

  vrow = 0;
  randu(acEst);

  int          size = 3;
  int          i    = 5;

  printmat("\nvrow:  \n", vrow);
  printmat("\nacEst: \n", acEst);

  vrow.Rows(1, size) = acEst.Column(i).Rows(1, size);
  
  printmat("\nvrow:  \n", vrow);
  cout << "The test case worked" << endl;
  return 0;
}
