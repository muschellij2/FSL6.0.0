#include "tests.hpp"


int main (int argc, char *argv[] ){

  // FFTI(    F, G, X, Y)
  // RealFFT( X, F, G)
  // RealFFTI(F, G, X)


  // FFT(     X, Y, F, G)

  ColumnVector X(24);
  ColumnVector Y(24);

  ColumnVector X2, Y2, X3;
  ColumnVector F,  G;

  randu(X);
  randu(Y);

  printmat("X: ", X);
  cout <<  X << endl;
  printmat("Y: ", Y);
  cout <<  Y << endl;

  FFT(X, Y, F, G);
  cout << "FFT(X, Y, F, G)" << endl;
  printmat("F: ", F);
  cout <<  F << endl;
  printmat("G: ", G);
  cout <<  G << endl;  

  FFTI(F, G, X2, Y2);
  cout << "FFTI(F, G, X2, Y2)" << endl;

  printmat("X2: ", X2);
  cout <<  X2 << endl;
  printmat("Y2: ", Y2);
  cout <<  Y2 << endl;


  RealFFT(X, F, G);
  cout << "RealFFT(X, F, G)" << endl;
  printmat("F: ", F);
  cout << F << endl;
  printmat("G: ", G);
  cout << G << endl;

  RealFFTI(F, G, X3);
  cout << "RealFFTI(F, G, X3)" << endl;
  printmat("X3: ", X3);
  cout <<  X3 << endl;
  


  return 0;
}
