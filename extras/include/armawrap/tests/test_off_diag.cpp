#include "tests.hpp"

Matrix data(5, 5);
Matrix tmp1(5, 5);
Matrix tmp2;

// This test is useless because I can't force the
// stack/heap to be filled with random values


void testDiagonalMatrix1() {

  // Assignment from a scalar
  DiagonalMatrix d(5);

  tmp2 = d * 1;

  d = 10;

  tmp1 << 10 << 0  << 0  << 0  << 0
       << 0  << 10 << 0  << 0  << 0
       << 0  << 0  << 10 << 0  << 0
       << 0  << 0  << 0  << 10 << 0
       << 0  << 0  << 0  << 0  << 10;

  tmp2 = d;

  if (tmp2 != tmp1) {
    cout << tmp2 << endl;
    throw std::string("testDiagonalMatrix1");
  }
}

void testDiagonalMatrix2() {

  // Assignment from another matrix
  DiagonalMatrix d;
  d = data;

  tmp1 << 1 << 0 << 0  << 0  << 0
       << 0 << 7 << 0  << 0  << 0
       << 0 << 0 << 13 << 0  << 0
       << 0 << 0 << 0  << 19 << 0
       << 0 << 0 << 0  << 0  << 25;

  tmp2 = d;
  if (tmp2 != tmp1) {
    cout << tmp2 << endl;
    throw std::string("testDiagonalMatrix2");
  }
}

void testDiagonalMatrix3() {

  DiagonalMatrix d;
  d = data.Row(2);

  tmp1 << 6 << 0 << 0  << 0  << 0
       << 0 << 7 << 0  << 0  << 0
       << 0 << 0 << 8  << 0  << 0
       << 0 << 0 << 0  << 9  << 0
       << 0 << 0 << 0  << 0  << 10;

  tmp2 = d;
  if (tmp2 != tmp1) {
    cout << tmp2 << endl;
    throw std::string("testDiagonalMatrix3");
  }
}

void testDiagonalMatrix4() {

  // Assignment from an expression
  DiagonalMatrix d;
  d = data * 10;

  tmp1 << 10 << 0  << 0   << 0   << 0
       << 0  << 70 << 0   << 0   << 0
       << 0  << 0  << 130 << 0   << 0
       << 0  << 0  << 0   << 190 << 0
       << 0  << 0  << 0   << 0   << 250;

  tmp2 = d;
  if (tmp2 != tmp1) {
    cout << tmp2 << endl;
    throw std::string("testDiagonalMatrix4");
  }
}

void testDiagonalMatrix5() {

  // Direct assignment;
  DiagonalMatrix d(5);
  for (int i = 1; i <= 5; i++) {
    d(i) = i;
  }

  tmp1 << 1 << 0 << 0 << 0 << 0
       << 0 << 2 << 0 << 0 << 0
       << 0 << 0 << 3 << 0 << 0
       << 0 << 0 << 0 << 4 << 0
       << 0 << 0 << 0 << 0 << 5;


  tmp2 = d;

  if (tmp2 != tmp1) {
    cout << tmp2 << endl;
    throw std::string("testDiagonalMatrix5");
  }
}

void testDiagonalMatrix6() {

  // Insertion
  DiagonalMatrix d(5);
  d << 1 << 2 << 3 << 4 << 5;

  tmp1 << 1 << 0 << 0 << 0 << 0
       << 0 << 2 << 0 << 0 << 0
       << 0 << 0 << 3 << 0 << 0
       << 0 << 0 << 0 << 4 << 0
       << 0 << 0 << 0 << 0 << 5;

  tmp2 = d;
  if (tmp2 != tmp1) {
    cout << tmp2 << endl;
    throw std::string("testDiagonalMatrix6");
  }
}


int main(int argc, char *argv[]) {

  data << 1  << 2  << 3  << 4  << 5
       << 6  << 7  << 8  << 9  << 10
       << 11 << 12 << 13 << 14 << 15
       << 16 << 17 << 18 << 19 << 20
       << 21 << 22 << 23 << 24 << 25;

  testDiagonalMatrix1();

  #if LIB == armawrap
  testDiagonalMatrix2();
  testDiagonalMatrix3();
  testDiagonalMatrix4();
  #endif

  testDiagonalMatrix5();
  testDiagonalMatrix6();
}
