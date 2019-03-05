#include "tests.hpp"

int main(int argc, char *argv[]) {

  BandMatrix          b(5, 1, 1);
  LowerBandMatrix     l(5, 2);
  UpperBandMatrix     u(5, 2);
  SymmetricBandMatrix s(5, 1);

  randu(b);
  randu(l);
  randu(u);
  randu(s);

  cout << "b: " << endl;
  cout <<  b    << endl;
  cout << "l: " << endl;
  cout <<  l    << endl;
  cout << "u: " << endl;
  cout <<  u    << endl;
  cout << "s: " << endl;
  cout <<  s    << endl; 

  return 0;
}
