#include "tests.hpp"


int main(int argc, char *argv[]) {

  Matrix       a(SIZE, SIZE);
  RowVector    r(SIZE);
  ColumnVector c(SIZE);

  randu(a);
  randu(r);
  randu(c);


  cout << "Before resize ... " << endl;
  cout << "a: " << Sum(a) << ", (" << a.Nrows() << ", " << a.Ncols() << ")" << endl;
  cout << "r: " << Sum(r) << ", (" << r.Nrows() << ", " << r.Ncols() << ")" << endl;
  cout << "c: " << Sum(c) << ", (" << c.Nrows() << ", " << c.Ncols() << ")" << endl;

  a.ReSize(SIZE / 2, SIZE / 2);
  r.ReSize(SIZE * 2);
  c.ReSize(SIZE / 3);

  randu(a);
  randu(r);
  randu(c);


  cout << "After first resize ... " << endl;
  cout << "a: " << Sum(a) << ", (" << a.Nrows() << ", " << a.Ncols() << ")" << endl;
  cout << "r: " << Sum(r) << ", (" << r.Nrows() << ", " << r.Ncols() << ")" << endl;
  cout << "c: " << Sum(c) << ", (" << c.Nrows() << ", " << c.Ncols() << ")" << endl;

  Matrix       aa(SIZE*3, 2);
  RowVector    rr(SIZE + 555);
  ColumnVector cc(SIZE / 4);

  a.ReSize(aa);
  r.ReSize(rr);
  c.ReSize(cc);

  randu(a);
  randu(r);
  randu(c);

  cout << "After second resize ... " << endl;
  cout << "a: " << Sum(a) << ", (" << a.Nrows() << ", " << a.Ncols() << ")" << endl;
  cout << "r: " << Sum(r) << ", (" << r.Nrows() << ", " << r.Ncols() << ")" << endl;
  cout << "c: " << Sum(c) << ", (" << c.Nrows() << ", " << c.Ncols() << ")" << endl;

  cout << "After CleanUp ... " << endl;
  a.CleanUp();
  cout << "a: " << Sum(a) << ", (" << a.Nrows() << ", " << a.Ncols() << ")" << endl;

  return 0;
}
