#include "tests.hpp"


void testAsScalar() {
  Matrix       asca(1, 1);
  RowVector    rsca(1);
  ColumnVector csca(1);

  randu(asca);
  randu(rsca);
  randu(csca);

  cout << "asca.AsScalar():  " << std::setprecision(4) << asca.AsScalar() << endl;
  cout << "rsca.AsScalar():  " << std::setprecision(4) << rsca.AsScalar() << endl;
  cout << "csca.AsScalar():  " << std::setprecision(4) << csca.AsScalar() << endl;
}

int main(int argc, char *argv[]) {

  #define SZ 16
  #define SZ_SQRT 4

  Matrix         a(SZ, SZ);
  DiagonalMatrix d(SZ);
  RowVector      r(SZ);
  ColumnVector   c(SZ);

  randu(a);
  randu(d);
  randu(r);
  randu(c);

  testAsScalar();

  printmat("a: ",  a);
  printmat("d: ",  d);
  printmat("r: ",  r);
  printmat("c: ",  c);

  RowVector      ar = a.AsRow();
  DiagonalMatrix ad = a.AsDiagonal();
  ColumnVector   ac = a.AsColumn();
  Matrix         am = a.AsMatrix(SZ_SQRT, SZ * SZ_SQRT);

  printmat("ar: ", ar);
  printmat("ad: ", ad);
  printmat("ac: ", ac);
  printmat("am: ", am);

  DiagonalMatrix rd = r.AsDiagonal();
  ColumnVector   rc = r.AsColumn();
  Matrix         rm = r.AsMatrix(2, SZ / 2);

  printmat("rd: ", rd);
  printmat("rc: ", rc);
  printmat("rm: ", rm);

  DiagonalMatrix cd = r.AsDiagonal();
  RowVector      cr = c.AsRow();
  Matrix         cm = c.AsMatrix(SZ_SQRT, SZ_SQRT);

  printmat("cd: ", cd);
  printmat("cr: ", cr);
  printmat("cm: ", cm);

  RowVector    dr = d.AsRow();
  ColumnVector dc = d.AsColumn();
  Matrix       dm = d.AsMatrix(SZ_SQRT, SZ_SQRT);

  printmat("dr: ", dr);
  printmat("dc: ", dc);
  printmat("dm: ", dm);
}
