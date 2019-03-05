#include "tests.hpp"


int main(int argc, char *argv[]) {

  Matrix       a(6, 6);
  ColumnVector p(6);
  ColumnVector q(6);
  RowVector    r(6);
  

  randu(a);
  randu(p);
  randu(q);
  randu(r);
 

  CroutMatrix c = a;

  cout << "c.i() *  p" << endl;          
  Matrix ap  = c.i() *  p;
  cout << "c.i() *  q" << endl;             
  Matrix aq  = c.i() *  q;
  cout << "c.i() * (p * r)" << endl;     
  Matrix apr = c.i() * (p * r);
  cout << "c.i()" << endl;               
  Matrix ai  = c.i();
  cout << "(c.i() * 1234.8) * p" << endl;
  Matrix ak  = (c.i() * 1234.8) * p;


  printmat("a: ",  a);
  cout << a << endl;

  printmat("ap: ", ap);
  cout << ap << endl;

  printmat("aq: ", aq);
  cout << aq << endl;

  printmat("apr: ", apr);
  cout << apr << endl; 
  
  printmat("ai: ", ai);
  cout << ai << endl;

  printmat("ak: ", ak);
  cout << ak << endl; 



  return 0;
}
