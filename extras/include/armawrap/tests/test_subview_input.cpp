#include "tests.hpp"


int main(int argc, char *argv[]) {
  
  Matrix                am( 4, 4);
  RowVector             ar( 8);
  ColumnVector          ac( 8);
  UpperTriangularMatrix aut(4);
  LowerTriangularMatrix alt(4);

  Real av1[] = {1, 2, 3, 4};
  //  Real av2[] = {9, 8, 7};
  Real av3[] = {21, 22};

  am  = 0;
  ar  = 0;
  ac  = 0;
  aut = 0;
  alt = 0;

  am.SubMatrix(2, 3, 2, 3) << av1;
  
  ar.Columns(3, 6) << av1;
  ac.Rows(   3, 6) << av1;

  ar.SubMatrix(1, 1, 1, 2) << av3;
  ac.SubMatrix(7, 8, 1, 1) << av3;

  am.Row(   1) << av1;
  am.Column(4) << av1;

  //alt.SubMatrix(2, 3, 2, 3) << av1;
  //aut.Row(   1) << av1;
  //aut.Row(   2) << av2;
  //aut.Column(3) << av3;

  printmat("am:  ", am);
  printmat("ar:  ", ar);
  printmat("ac:  ", ac);
  printmat("aut: ", aut);
  printmat("alt: ", alt);
  
  cout << "am:  " << endl;
  cout << am  << endl;
  cout << "ar:  " << endl;
  cout << ar  << endl;
  cout << "ac:  " << endl;
  cout << ac  << endl;
  cout << "aut: " << endl;
  cout << aut << endl;
  cout << "alt: " << endl;
  cout << alt << endl;


  return 0;
}
