#include "tests.hpp"


ReturnMatrix unwrapMatrix(const Matrix & m)
{
	ColumnVector munwrap(m.Nrows()*m.Ncols());
	 int count=0;
	for (int i =0; i<m.Nrows() ; i++)
		for (int j =0; j<m.Ncols() ; j++,count++)
			munwrap.element(count)=m.element(i,j);
	return munwrap;
}


ReturnMatrix matrixToDiagonalMatrix(Matrix & mat)
{
  ColumnVector vec = unwrapMatrix(mat);

	DiagonalMatrix out(vec.Storage());

	for (int row = 1; row <= vec.Storage(); row++)
		out(row) = vec(row);

	out.Release();
	return out;
}


int main(int argc, char *argv[]) {

  Matrix a(3, 4);

  randu(a);

  DiagonalMatrix m = matrixToDiagonalMatrix(a);

  printmat("a: ", a);
  printmat("m: ", m);

  cout << "a: " << endl;
  cout <<  a    << endl;
  cout << "m: " << endl;
  cout <<  m    << endl;

  return 0;
}
