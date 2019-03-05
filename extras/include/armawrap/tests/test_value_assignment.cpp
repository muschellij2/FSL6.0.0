
#include "tests.hpp"

int main(int argc, char *argv[]) {

  Matrix                a(5, 5);
  RowVector             r(5);
  ColumnVector          c(5);
  IdentityMatrix        i(5);
  DiagonalMatrix        d(5);
  UpperTriangularMatrix ut(5);
  LowerTriangularMatrix lt(5);
  SymmetricMatrix       s(5);
  BandMatrix            b(5, 2, 2);
  UpperBandMatrix      ub(5, 2);
  LowerBandMatrix      lb(5, 2);

  randu(a);
  randu(r);
  randu(c);
  randu(d);
  randu(i);
  randu(ut);
  randu(lt);
  randu(s);
  randu(b);
  randu(ub);
  randu(lb);

  a(3, 2) = 1843;
  r(4)    = 12346;
  c(5)    = 235364;
  d(3)    = 324890;
  ut(3,3) = 123;
  ut(3,4) = 456;

  lt(5,4) = 888;
  lt(4,4) = 999;
  
  s(3, 3) = 33;
  s(3, 4) = 34;
  s(5, 2) = 52;

  b(4, 4) = 99;
  b(1, 1) = 99;
  b(1, 3) = 99;
  b(3, 1) = 99;
  b(5, 3) = 99;

  ub(3, 3) = 99;
  ub(3, 4) = 99;
  ub(3, 5) = 99;
  ub(1, 3) = 99;

  lb(3, 3) = 99;
  lb(4, 3) = 99;
  lb(5, 3) = 99;
  lb(3, 1) = 99;

  printmat("a:  ", a);
  printmat("r:  ", r);
  printmat("c:  ", c);
  printmat("i:  ", i);
  printmat("d:  ", d);
  printmat("ut: ", ut);
  printmat("lt: ", lt);
  printmat("s:  ", s);
  printmat("b:  ", b);
  printmat("ub: ", ub);
  printmat("lb: ", lb);

  cout << "a: "  << endl;
  cout <<  a     << endl;
  cout << "r: "  << endl;
  cout <<  r     << endl;
  cout << "c: "  << endl;
  cout <<  c     << endl;
  cout << "i: "  << endl;
  cout <<  i     << endl;
  cout << "d: "  << endl;
  cout <<  d     << endl;
  cout << "ut: " << endl;
  cout <<  ut    << endl;
  cout << "lt: " << endl;
  cout <<  lt    << endl;
  cout << "s: "  << endl;
  cout <<  s     << endl;
  cout << "b: "  << endl;
  cout <<  b     << endl;
  cout << "ub: " << endl;
  cout <<  ub    << endl;
  cout << "lb: " << endl;
  cout <<  lb    << endl;

  cout << "a(3, 2): " << a(3, 2)  << endl;
  cout << "r(4):    " << r(4)     << endl;
  cout << "c(5):    " << c(5)     << endl;
  cout << "d(3):    " << d(3)     << endl;
  cout << "d(2, 2): " << d(2, 2)  << endl;
  cout << "ut(3, 3):" << ut(3, 3) << endl;
  cout << "ut(3, 4):" << ut(3, 4) << endl;
  cout << "lt(4, 4):" << lt(4, 4) << endl;
  cout << "lt(5, 4):" << lt(5, 4) << endl;
  cout << "s(3, 3): " << s(3, 3)  << endl;
  cout << "s(3, 4): " << s(3, 4)  << endl;
  cout << "s(4, 3): " << s(4, 3)  << endl;
  cout << "s(5, 2): " << s(5, 2)  << endl;
  cout << "s(2, 5): " << s(2, 5)  << endl;


  cout << "max(s(2, 5), lt(5, 4)): " << std::max<double>(s(2, 5), lt(5, 4)) << endl;
}
