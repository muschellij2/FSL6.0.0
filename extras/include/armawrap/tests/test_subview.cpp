#include "tests.hpp"

int main(int argc, char *argv[]) {

  Matrix         a(SIZE, SIZE); 
  DiagonalMatrix d(SIZE); 
  RowVector      r(SIZE);
  ColumnVector   c(SIZE);
  Matrix        ar(SIZE, SIZE*2);

  Matrix mi = IdentityMatrix(11);

  randu(a);
  randu(d);
  randu(r);
  randu(c);
  randu(ar);

  Matrix       smaa  = a.SubMatrix(20, 30, 30, 40);
  Matrix       smab  = a.SymSubMatrix(40, 60);
  RowVector    smac  = a.Row(50);
  ColumnVector smad  = a.Column(38);
  Matrix       smacd = SP(a.Column(38) * smac, a);
  RowVector    srar  = ar.Row(50);
  ColumnVector srac  = ar.Column(190);
  a.SubMatrix(50, 60, 60, 70) = mi;

  smacd.SubMatrix( 4,  8, 4,  8) = a.SubMatrix(4, 8, 4, 8);
  smacd.SubMatrix(10, 14, 8, 12) = d.SymSubMatrix(10, 14);
  
  printmat("a: ",     a);
  printmat("mi: ",    mi);
  printmat("r: ",     r);
  printmat("c: ",     c);
  printmat("smaa: ",  smaa);
  printmat("smab: ",  smab);
  printmat("smac: ",  smac);
  printmat("srar: ",  srar);
  printmat("srac: ",  srac);
  printmat("smad: ",  smad);
  printmat("smacd: ", smacd);

  cout << "a.SubMatrix(4, 8, 4, 8): "        << endl;
  cout <<  a.SubMatrix(4, 8, 4, 8)           << endl;
  cout << "smacd.SubMatrix(4, 8, 4, 8): "    << endl;
  cout <<  smacd.SubMatrix(4, 8, 4, 8)       << endl;
  cout << "d.SymSubMatrix(10, 14): "         << endl;
  cout <<  d.SymSubMatrix(10, 14)            << endl;
  cout << "smacd.SubMatrix(10, 14, 8, 12): " << endl;
  cout <<  smacd.SubMatrix(10, 14, 8, 12)    << endl;



  return 0;
}
