#include "tests.hpp"

int main(int argc, char *argv[]) {

  Matrix                m(   10, 10);
  IdentityMatrix        i(   10);
  DiagonalMatrix        d(   10);
  RowVector             r(   10);
  ColumnVector          c(   10);
  UpperTriangularMatrix ut(  10);
  LowerTriangularMatrix lt(  10);
  SymmetricMatrix       s(   10);
  Matrix                tmp1(4, 4);
  Matrix                tmp2(3, 6);
  RowVector             tmp3(3);
  ColumnVector          tmp4(4);
  Matrix                tmp5(2, 2);

  randu(m);
  randu(d);
  randu(r);
  randu(c);
  randu(s);
  randu(ut);
  randu(lt);
  randu(tmp1);
  randu(tmp2);
  randu(tmp3);
  randu(tmp4);

  tmp5 = 850;
  
  m.SubMatrix(2, 5,  4, 7) = tmp1 + 50;
  r.SubMatrix(1, 1,  1, 3) = tmp3 + 250;
  c.SubMatrix(7, 10, 1, 1) = tmp4 + 250;

  m.Row(    9)     = 9;
  r.Columns(9, 10) = 9;
  c.Rows(   2, 3)  = 9;

  // armawrap supports subview assignment of 
  // special matrix, but newmat doesn't.
  #if LIB == ARMAWRAP_API
  d .SubMatrix(2, 4, 4, 9) = tmp2 - 100;
  ut.SubMatrix(2, 5, 4, 7) = tmp1 - 50;
  lt.SubMatrix(4, 7, 3, 6) = tmp1 + 200;
  s .SubMatrix(4, 6, 3, 8) = tmp2 + 999;

  d .Rows(   1, 2)  = 5;
  ut.Rows(   8, 10) = 5;
  lt.Columns(6, 7)  = 5;
  s .Rows(   4, 5)  = 5;
  #else
  for (int row = 2; row <= 4; row ++) {
    for (int col = 4; col <= 9; col++) { if (row == col) d( row, col) = tmp2(row-1, col-3) - 100; } }
  for (int row = 2; row <= 5; row ++) {
    for (int col = 4; col <= 7; col++) { if (row <= col) ut(row, col) = tmp1(row-1, col-3) - 50; } }
  for (int row = 4; row <= 7; row ++) {
    for (int col = 3; col <= 6; col++) { if (row >= col) lt(row, col) = tmp1(row-3, col-2) + 200; } }
  for (int row = 4; row <= 6; row ++) {
    for (int col = 3; col <= 8; col++) { if (row >= col) s( row, col) = tmp2(row-3, col-2) + 999; } }

  d( 1,  1)  = 5;
  d( 2,  2)  = 5;
  ut(8,  8)  = 5;
  ut(8,  9)  = 5;
  ut(8,  10) = 5;
  ut(9,  9)  = 5;
  ut(9,  10) = 5;
  ut(10, 10) = 5;
  for (int row = 6; row <= 10; row++) {
    for (int col = 6; col <= 7; col++) { if (row >= col) lt(row, col) = 5; } }
  for (int row = 4; row <= 5; row++) {
    for (int col = 1; col <= 10; col++) { if (row >= col) s(row, col) = 5; } } 
  #endif

  m .Rows(        1, 2)        *= 20;
  r .Columns(     4, 6)        -= 50;
  c .Rows(        7, 8)        += 50;
  #if LIB == ARMAWRAP_API  
  d .SymSubMatrix(7, 8)        /= 100000;
  ut.SubMatrix(   1, 2, 9, 10) += tmp5;
  lt.SubMatrix(   9, 10, 1, 2) -= tmp5;
  #else
  d(7, 7) /= 100000;
  d(8, 8) /= 100000;
  ut.SubMatrix(   1, 2, 9, 10) += tmp5;
  lt.SubMatrix(   9, 10, 1, 2) -= tmp5; 
  #endif

  m.Column(1) = m.Row(5).t();

  cout << "i: "                    << endl;
  cout <<  i                       << endl;
  cout << "i.SymSubMatrix(3, 6): " << endl;
  cout <<  i.SymSubMatrix(3, 6)    << endl;


  cout << "m: "  << endl;
  cout <<  m     << endl; 
  cout << "d: "  << endl;
  cout <<  d     << endl;
  cout << "r: "  << endl;
  cout <<  r     << endl; 
  cout << "c: "  << endl;
  cout <<  c     << endl; 
  cout << "ut: " << endl;
  cout <<  ut    << endl; 
  cout << "lt: " << endl;
  cout <<  lt    << endl; 
  cout << "m.SubMatrix( 1, 10, 1, 10): " << endl;
  cout <<  m.SubMatrix( 1, 10, 1, 10)    << endl;
  cout << "d.SubMatrix( 1, 10, 1, 10): " << endl;
  cout <<  d.SubMatrix( 1, 10, 1, 10)    << endl;
  cout << "ut.SubMatrix(1, 10, 1, 10): " << endl;
  cout <<  ut.SubMatrix(1, 10, 1, 10)    << endl;
  cout << "lt.SubMatrix(1, 10, 1, 10): " << endl;
  cout <<  lt.SubMatrix(1, 10, 1, 10)    << endl;
  cout << "s.SubMatrix( 1, 10, 1, 10): " << endl;
  cout <<  s.SubMatrix( 1, 10, 1, 10)    << endl;
  return 0;
}
