#include "tests.hpp"



using namespace NEWMAT;

template<class S, class T>
inline T Max(const S &a, const T &b) { if (a>b) return (T) a; else return b; }

ReturnMatrix pinv(const Matrix& mat2)
{
  // calculates the psuedo-inverse using SVD
  // note that the right-pinv(x') = pinv(x).t()
  Matrix mat(mat2);
  if ( mat2.Ncols() > mat2.Nrows() )
    mat=mat.t();
  Tracer tr("pinv");
  DiagonalMatrix D;
  Matrix U, V;
  cout << "SVD call .. " << endl;
  SVD(mat,D,U,V);
  cout << "Finished" << endl;
  printmat("mat: ", mat);
  cout <<   mat  << endl;
  printmat("D: ", D);
  cout <<   D  << endl;
  abs(U);
  abs(V);
  printmat("U: ", U);
  cout <<   U  << endl;
  printmat("V: ", V);
  cout <<   V  << endl;
  float tol;
  tol = MaximumAbsoluteValue(D) * Max(mat.Nrows(),mat.Ncols()) * 1e-16;
  for (int n=1; n<=D.Nrows(); n++) {
    if (fabs(D(n,n))>tol) D(n,n) = 1.0/D(n,n);
    else D(n,n) = 0.0; // reduce the number of columns because too close to singular
  }
  Matrix pinv = V * D * U.t();
  if ( mat2.Ncols() > mat2.Nrows() )
    pinv=pinv.t();
  pinv.Release();
  return pinv;
}

int main(int argc, char *argv[]) {

  Matrix       PM(12, 4);
  ColumnVector  P(12);
  ColumnVector RP(12);

  randu(PM);

  Matrix NP = pinv(PM) * (P - RP);

  for (int ri = 1; ri <= NP.Nrows(); ri++) {
    for (int ci = 1; ci <= NP.Ncols(); ci++) {
      if (fabs(NP(ri, ci)) < 0.00001) {
        NP(ri, ci) = 0;
      }
    }
  }

  cout << NP << endl;

  return 0;
}
