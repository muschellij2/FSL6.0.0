#include "tests.hpp"

int main (int argc, char *argv[]) {
  Matrix       a(SIZE, SIZE);
  ColumnVector c(SIZE);
  RowVector    r(SIZE);

  randu(a);
  randu(c);
  randu(r);

  Matrix spa    = SP(a, a);
  Matrix spc    = SP(c, c);
  Matrix spr    = SP(r, r);
  Matrix suba   = a.SymSubMatrix(50, 100);
  Matrix kpa    = KP(a.SymSubMatrix(50, 100), a.SymSubMatrix(50, 100));
  Matrix kpsuba = KP(suba, suba);
  Matrix kpc    = KP(a, c);
  Matrix kpr    = KP(r, a);
  double aa     = DotProduct(a, a);
  double rr     = DotProduct(r, r);
  double cc     = DotProduct(c, c);
  double cr     = DotProduct(c, r);

  Matrix ai       = a.i();
  Matrix at       = a.t();
  RowVector ct    = c.t();
  ColumnVector rt = r.t();

  printmat("spa:      ", spa);
  printmat("spc:      ", spc);
  printmat("spr:      ", spr);
  printmat("kpa:      ", kpa);
  printmat("suba:     ", suba);
  printmat("kpsuba:   ", kpsuba);
  printmat("kpc:      ", kpc);
  printmat("kpr:      ", kpr);

  cout <<  "aa:       " << aa << endl;
  cout <<  "cc:       " << cc << endl;
  cout <<  "rr:       " << rr << endl;
  cout <<  "cr:       " << cr << endl;

  printmat("ai:       ", ai);
  printmat("at:       ", at);
  printmat("ct:       ", ct);
  printmat("rt:       ", rt);

  return 0;

}
