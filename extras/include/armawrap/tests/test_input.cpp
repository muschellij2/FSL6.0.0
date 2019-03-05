#include "tests.hpp"

int main(int argc, char *argv[]) {
  
  Matrix       a( 4, 3);
  Matrix       as(9, 9);
  RowVector    r( 8);
  ColumnVector c( 4);

  c << 0  << 0  << 0  << 1;
  //c << 20 << 21 << 22 << 23;

  

  a << 1  << 2  << 3 
    << 4  << 5  << 6
    << 7  << 8  << 9
    << 10 << 11 << 12;

  r << 10 << 11 << 12 << 13 << 14 << 15 << 16 << 17;
  randu(as);

  cout << "as.SubMatrix(3, 5, 3, 5): " << endl;
  cout << setw(10) << setprecision(5) << as.SubMatrix(3, 5, 3, 5) << endl;

  // the following line is identical to as.Row(3) << ...;
  as.SubMatrix(3, 3, 1, 9) << 11 << 12 << 13 << 14 << 15 << 16 << 17 << 18 << 19;
  as.Row(4)                << 21 << 22 << 23 << 24 << 25 << 26 << 27 << 28 << 29;
  as.Row(5)                << 31 << 32 << 33 << 34 << 35 << 36 << 37 << 38 << 39;

  printmat("a:  ", a);
  printmat("as: ", as);
  printmat("r:  ", r);
  printmat("c:  ", c);

  cout << "a: " << endl;
  cout << setw(10) << setprecision(5) << a << endl;
  cout << "as: " << endl;
  cout << setw(10) << setprecision(5) << as << endl;
  cout << "r: " << endl;
  cout << setw(10) << setprecision(5) << r << endl;
  cout << "c: " << endl;
  cout << setw(10) << setprecision(5) << c << endl;

  Matrix                 aa(3, 3);
  RowVector              ar(9);
  ColumnVector           ac(9);
  IdentityMatrix         ai(3);
  DiagonalMatrix         ad(3);
  UpperTriangularMatrix aut(3);
  LowerTriangularMatrix alt(3);

  double vals1[] = {9,8,7,6,5,4,3,2,1};
  double vals2[] = {99};
  double vals3[] = {99, 98, 97};
  double vals4[] = {99, 98, 97, 96, 95, 94};
  double vals5[] = {94, 95, 96, 97, 98, 99};
  aa << vals1;
  ar << vals1;
  ac << vals1;

  ai  << vals2;
  ad  << vals3;
  aut << vals4;
  alt << vals5;

  printmat("aa:  ", aa);
  printmat("ar:  ", ar);
  printmat("ac:  ", ac);
  printmat("ai:  ", ai);
  printmat("ad:  ", ad);
  printmat("aut: ", aut);
  printmat("alt: ", alt);

  cout << "aa: " << endl;
  cout << setw(10) << setprecision(5) << aa << endl;
  cout << "ar: " << endl;
  cout << setw(10) << setprecision(5) << ar << endl;
  cout << "ac: " << endl;
  cout << setw(10) << setprecision(5) << ac << endl;
  cout << "ai: " << endl;
  cout << setw(10) << setprecision(5) << ai << endl;
  cout << "ad: " << endl;
  cout << setw(10) << setprecision(5) << ad << endl;
  cout << "aut: " << endl;
  cout << setw(10) << setprecision(5) << aut << endl;
  cout << "alt: " << endl;
  cout << setw(10) << setprecision(5) << alt << endl;


  IdentityMatrix        i(4);
  DiagonalMatrix        d(4);
  UpperTriangularMatrix ut(4);
  LowerTriangularMatrix lt(4);

  i  << 2;
  d  << 1  << 2 << 3 << 4;
  ut << 1  << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9 << 10;
  lt << 10 << 9 << 8 << 7 << 6 << 5 << 4 << 3 << 2 << 1;

  printmat("i:  ", i);
  printmat("d:  ", d);
  printmat("ut: ", ut);
  printmat("lt: ", lt);

  cout << "i:  " << i  << endl;
  cout << "d:  " << d  << endl;
  cout << "ut: " << ut << endl;
  cout << "lt: " << lt << endl;


  Matrix       em(4, 4);
  RowVector    er(4);
  ColumnVector ec(4);

  randu(er);
  randu(ec);

  em << ec * er;

  cout << "er: " << endl;
  cout <<  er    << endl;
  cout << "ec: " << endl;
  cout <<  ec    << endl;
  cout << "em: " << endl;
  cout <<  em    << endl;

  return 0;
}
