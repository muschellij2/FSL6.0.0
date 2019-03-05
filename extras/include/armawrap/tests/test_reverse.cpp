//#define ARMA_EXTRA_DEBUG
#include "tests.hpp"

int main(int argc, char *argv[]) {

  Matrix m(4, 3);
  randu(m);
  Matrix b = m.Reverse().Reverse();


  cout << SP(((b.Reverse() * 54) + 34), m) << endl;
  
  //  cout << "b" << endl;
  
  cout << b      << endl;

  //  cout << "in line" << endl;

  cout << m.Reverse().Reverse() << endl;

  return 0;
}
