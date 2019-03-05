#include "tests.hpp"

int main(int argc, char *argv[]) {

  Matrix       m(43, 54);
  Matrix       z(43, 54);
  ColumnVector c(54);
  RowVector    r(43);

  randu(m);
  randu(z);
  randu(r);
  randu(c);

  Matrix         mrc = SP(m * c, r.AsColumn());
  IdentityMatrix im(48);

  im = 999;

  cout << "m:   " << m  .Sum() << " (" << m  .Nrows() << ", " << m  .Ncols() << ")" << endl;
  cout << "r:   " << r  .Sum() << " (" << r  .Nrows() << ", " << r  .Ncols() << ")" << endl;
  cout << "c:   " << c  .Sum() << " (" << c  .Nrows() << ", " << c  .Ncols() << ")" << endl;
  cout << "mrc: " << mrc.Sum() << " (" << mrc.Nrows() << ", " << mrc.Ncols() << ")" << endl;
  cout << "im:  " << im .Sum() << " (" << im .Nrows() << ", " << im .Ncols() << ")" << endl;

  return 0;
}
