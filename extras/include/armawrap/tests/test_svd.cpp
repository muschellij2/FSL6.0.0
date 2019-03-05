#include "tests.hpp"

void printAll(const char *msg, Matrix &A, DiagonalMatrix &D, Matrix &U, Matrix &V) {
  cout << msg   << endl;
  cout << "D: " << endl;
  cout <<  D    << endl;
  cout << "U: " << endl;
  cout <<  U    << endl;
  cout << "V: " << endl;
  cout <<  V    << endl;

  D = 0;
  U = 0;
  V = 0;
}

int main(int argc, char *argv[]) {

  Matrix         A(5, 5);
  DiagonalMatrix D(5);
  Matrix         U(5, 5);
  Matrix         V(5, 5);

  randu(A);

  cout << "A: " << endl;
  cout <<  A    << endl << endl;;

  cout << "NOTE: For SVD(..., false) and SVD(..., false, X), "    <<
          "Newmat uses the U matrix as a workspace. Results may " <<
          "therefore differ between armawrap and newmat."         << endl;

  SVD(A, D);                     printAll("SVD(A, D)",                     A, D, U, V); 
  SVD(A, D, U);                  printAll("SVD(A, D, U)",                  A, D, U, V);
  SVD(A, D, U,    false);        printAll("SVD(A, D, U, false)",           A, D, U, V);
  SVD(A, D, U,    true);         printAll("SVD(A, D, U, true)",            A, D, U, V);
  SVD(A, D, U, V);               printAll("SVD(A, D, U, V)",               A, D, U, V);
  SVD(A, D, U, V, false);        printAll("SVD(A, D, U, V, false)",        A, D, U, V);
  SVD(A, D, U, V, true);         printAll("SVD(A, D, U, V, true)",         A, D, U, V);
  SVD(A, D, U, V, false, false); printAll("SVD(A, D, U, V, false, false)", A, D, U, V);
  SVD(A, D, U, V, false, true);  printAll("SVD(A, D, U, V, false, true)",  A, D, U, V);
  SVD(A, D, U, V, true,  false); printAll("SVD(A, D, U, V, true,  false)", A, D, U, V);
  SVD(A, D, U, V, true,  true);  printAll("SVD(A, D, U, V, true,  true)",  A, D, U, V);


  cout << "Passing uninitialised matrices .. " << endl;
  Matrix         AA(5, 5);
  DiagonalMatrix DD;
  Matrix         UU;
  Matrix         VV;

  randu(AA);

  SVD(AA, DD, UU, VV);

  cout << "AA" << endl;
  cout <<  AA  << endl;
  cout << "DD" << endl;
  cout <<  DD  << endl;
  cout << "UU" << endl;
  cout <<  UU  << endl;
  cout << "VV" << endl;
  cout <<  VV  << endl; 

  return 0;
}
