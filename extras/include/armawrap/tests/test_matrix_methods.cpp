#include <math.h>
#include "tests.hpp"

int main(int argc, char *argv[]) {

  Matrix       a = Matrix(SIZE, SIZE);
  RowVector    r = RowVector(SIZE);
  ColumnVector c = ColumnVector(SIZE);
  Matrix       z(SIZE, SIZE);

  z = 0;
  
  randu(a);
  randu(r);
  randu(c);

  cout << "a.SumSquare():        " << a.SumSquare()        << endl;
  cout << "SumSquare(a):         " << SumSquare(a)         << endl;
  cout << "r.SumSquare():        " << r.SumSquare()        << endl;
  cout << "SumSquare(r):         " << SumSquare(r)         << endl; 
  cout << "c.SumSquare():        " << c.SumSquare()        << endl;
  cout << "SumSquare(c):         " << SumSquare(c)         << endl; 
  cout << "a.SumAbsoluteValue(): " << a.SumAbsoluteValue() << endl;
  cout << "SumAbsoluteValue(a):  " << SumAbsoluteValue(a)  << endl;
  cout << "r.SumAbsoluteValue(): " << r.SumAbsoluteValue() << endl;
  cout << "SumAbsoluteValue(r):  " << SumAbsoluteValue(r)  << endl; 
  cout << "c.SumAbsoluteValue(): " << c.SumAbsoluteValue() << endl;
  cout << "SumAbsoluteValue(c):  " << SumAbsoluteValue(c)  << endl; 
  cout << "a.Sum():              " << a.Sum()              << endl;
  cout << "Sum(a):               " << Sum(a)               << endl;
  cout << "Sum(SP(a, a)):        " << Sum(SP(a, a))        << endl;
  cout << "r.Sum():              " << r.Sum()              << endl;
  cout << "Sum(r):               " << Sum(r)               << endl; 
  cout << "c.Sum():              " << c.Sum()              << endl;
  cout << "Sum(c):               " << Sum(c)               << endl; 
  cout << "a.Trace():            " << a.Trace()            << endl;
  cout << "Trace(a):             " << Trace(a)             << endl;
  cout << "a.Norm1():            " << a.Norm1()            << endl;
  cout << "Norm1(a):             " << Norm1(a)             << endl;
  cout << "r.Norm1():            " << r.Norm1()            << endl;
  cout << "Norm1(r):             " << Norm1(r)             << endl; 
  cout << "c.Norm1():            " << c.Norm1()            << endl;
  cout << "Norm1(c):             " << Norm1(c)             << endl; 
  cout << "a.NormInfinity():     " << a.NormInfinity()     << endl;
  cout << "NormInfinity(a):      " << NormInfinity(a)      << endl;
  cout << "r.NormInfinity():     " << r.NormInfinity()     << endl;
  cout << "NormInfinity(r):      " << NormInfinity(r)      << endl; 
  cout << "c.NormInfinity():     " << c.NormInfinity()     << endl;
  cout << "NormInfinity(c):      " << NormInfinity(c)      << endl; 
  cout << "a.NormFrobenius():    " << a.NormFrobenius()    << endl;
  cout << "NormFrobenius(a):     " << NormFrobenius(a)     << endl;
  cout << "r.NormFrobenius():    " << r.NormFrobenius()    << endl;
  cout << "NormFrobenius(r):     " << NormFrobenius(r)     << endl; 
  cout << "c.NormFrobenius():    " << c.NormFrobenius()    << endl;
  cout << "NormFrobenius(c):     " << NormFrobenius(c)     << endl; 
  cout << "a.Determinant():      " << std::setprecision(4) << a.Determinant()      << endl;
  cout << "Determinant(c):       " << std::setprecision(4) << Determinant(a)       << endl; 
  LogAndSign res1 = a.LogDeterminant();
  LogAndSign res2 = LogDeterminant(a);
  cout << "a.LogDeterminant():   " << res1.Sign() << "," << std::setprecision(3) <<  res1.LogValue() << endl;
  cout << "LogDeterminant(a):    " << res2.Sign() << "," << std::setprecision(3) <<  res2.LogValue() << endl;
  cout << "a.IsZero():           " << a.IsZero()                             << endl; 
  cout << "IsZero(a):            " << IsZero(a)                              << endl; 
  cout << "z.IsZero():           " << z.IsZero()                             << endl; 
  cout << "IsZero(z):            " << IsZero(z)                              << endl; 


  cout << "Before reverse..." << endl;
  printmat("a:   ", a);
  printmat("r:   ", r);
  printmat("c:   ", c);
  cout << "a first 3: " << a(1,    1)      << "," << a(1,    2)      << "," << a(1,    3)    << endl;
  cout << "a last 3:  " << a(SIZE, SIZE-2) << "," << a(SIZE, SIZE-1) << "," << a(SIZE, SIZE) << endl;
  cout << "r first 3: " << r(1)            << "," << r(2)            << "," << r(3)          << endl;
  cout << "r last 3:  " << r(SIZE-2)       << "," << r(SIZE-1)       << "," << r(SIZE)       << endl;
  cout << "c first 3: " << c(1)            << "," << c(2)            << "," << c(3)          << endl;
  cout << "c last 3:  " << c(SIZE-2)       << "," << c(SIZE-1)       << "," << c(SIZE)       << endl; 

  a = a.Reverse();
  r = r.Reverse();
  c = c.Reverse();

  cout << "After reverse..." << endl;
  printmat("a:   ", a);
  printmat("r:   ", r);
  printmat("c:   ", c);

  cout << "a first 3: " << a(1,    1)      << "," << a(1,    2)      << "," << a(1,    3)    << endl;
  cout << "a last 3:  " << a(SIZE, SIZE-2) << "," << a(SIZE, SIZE-1) << "," << a(SIZE, SIZE) << endl;
  cout << "r first 3: " << r(1)            << "," << r(2)            << "," << r(3)          << endl;
  cout << "r last 3:  " << r(SIZE-2)       << "," << r(SIZE-1)       << "," << r(SIZE)       << endl;
  cout << "c first 3: " << c(1)            << "," << c(2)            << "," << c(3)          << endl;
  cout << "c last 3:  " << c(SIZE-2)       << "," << c(SIZE-1)       << "," << c(SIZE)       << endl; 
}
