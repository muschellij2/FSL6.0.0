#include "tests.hpp"

using namespace NEWMAT;


int main(int argc, char *argv[]) {

  SymmetricMatrix s;
  DiagonalMatrix  d;
  RowVector       r(5);
  Matrix          m(5, 5);
  
  randu(r);
  randu(m);
  
  s << m;

  d << r;

  printmat("m: ", m);
  cout <<   m << endl;
  printmat("r: ", r);
  cout <<   r << endl;
  printmat("s: ", s);
  cout <<   s << endl;
  printmat("d: ", d);
  cout <<   d << endl;


  DiagonalMatrix dd(10);
  randu(dd);

  printmat("dd (before self insert): ", dd);
  cout <<   dd << endl; 

  dd << dd.SymSubMatrix(3, 7);

  printmat("dd (after self insert): ", dd);
  cout <<   dd << endl;
  
  return 0;
}
