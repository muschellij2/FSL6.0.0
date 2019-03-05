#include "tests.hpp"


int main(int argc, char *argv[]) {

  Matrix                Ao(8, 5);
  Matrix                A(8, 5);
  UpperTriangularMatrix U(8);

  randu(Ao);

  A = Ao;

  cout << "A: " << endl;
  cout <<  A    << endl;

  cout << "QRZ(A, U): "         << endl;
  QRZ(A, U);

  abs(A);
  abs(U);
  cout << "A (after): "         << endl;
  cout <<  A                    << endl;
  cout << "U: "                 << endl;
  cout <<  U                    << endl;

  return 0;
}
