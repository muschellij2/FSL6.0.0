#include <cassert>

#include "tests.hpp"


#if LIB == ARMAWRAP_API
int main(int argc, char *argv[]) {

  Matrix m(4, 4);

  assert(m.Row(1)    .Storage() == 4);
  assert(m.Row(2)    .Storage() == 4);
  assert(m.Row(3)    .Storage() == 4);
  assert(m.Row(4)    .Storage() == 4);
  assert(m.Column(1) .Storage() == 4);
  assert(m.Column(2) .Storage() == 4);
  assert(m.Column(3) .Storage() == 4);
  assert(m.Column(4) .Storage() == 4);
  assert(m.Rows(1, 2).Storage() == 8);
  assert(m.Rows(2, 4).Storage() == 12);

  assert(m.Columns(1, 2).Storage() == 8);
  assert(m.Columns(2, 4).Storage() == 12);

  assert(m.SubMatrix(1, 2, 1, 2).Storage() == 4);
  assert(m.SubMatrix(1, 3, 1, 2).Storage() == 6);
  assert(m.SubMatrix(1, 3, 1, 3).Storage() == 9);

  ColumnVector c(6);
  assert(c.Row(3).Storage()     == 1);
  assert(c.Rows(2, 4).Storage() == 3);

  RowVector r(6);
  assert(r.Column(3).Storage()     == 1);
  assert(r.Columns(2, 4).Storage() == 3);

  SymmetricMatrix s(5);
  assert(s.Row(1)   .Storage() == 1);
  assert(s.Row(2)   .Storage() == 2);
  assert(s.Row(3)   .Storage() == 3);
  assert(s.Row(4)   .Storage() == 4);
  assert(s.Row(5)   .Storage() == 5);
  assert(s.Column(1).Storage() == 5);
  assert(s.Column(2).Storage() == 4);
  assert(s.Column(3).Storage() == 3);
  assert(s.Column(4).Storage() == 2);
  assert(s.Column(5).Storage() == 1);

  assert(s.Rows(   1, 2).Storage() == 3);
  assert(s.Rows(   1, 3).Storage() == 6);
  assert(s.Rows(   1, 4).Storage() == 10);
  assert(s.Columns(1, 2).Storage() == 9);
  assert(s.Columns(1, 3).Storage() == 12);
  assert(s.Columns(1, 4).Storage() == 14);

  assert(s.SubMatrix(1, 2, 1, 2).Storage() == 3);
  assert(s.SubMatrix(1, 2, 1, 3).Storage() == 3);
  assert(s.SubMatrix(1, 3, 1, 2).Storage() == 5);

  assert(s.SubMatrix(2, 3, 2, 3).Storage() == 3);
  assert(s.SubMatrix(2, 3, 2, 4).Storage() == 3);
  assert(s.SubMatrix(2, 4, 2, 3).Storage() == 5);

  LowerTriangularMatrix lt(5);
  assert(lt.Row(1)   .Storage() == 1);
  assert(lt.Row(2)   .Storage() == 2);
  assert(lt.Row(3)   .Storage() == 3);
  assert(lt.Row(4)   .Storage() == 4);
  assert(lt.Row(5)   .Storage() == 5);
  assert(lt.Column(1).Storage() == 5);
  assert(lt.Column(2).Storage() == 4);
  assert(lt.Column(3).Storage() == 3);
  assert(lt.Column(4).Storage() == 2);
  assert(lt.Column(5).Storage() == 1);

  assert(lt.Rows(   1, 2).Storage() == 3);
  assert(lt.Rows(   1, 3).Storage() == 6);
  assert(lt.Rows(   1, 4).Storage() == 10);
  assert(lt.Columns(1, 2).Storage() == 9);
  assert(lt.Columns(1, 3).Storage() == 12);
  assert(lt.Columns(1, 4).Storage() == 14);

  assert(lt.SubMatrix(1, 2, 1, 2).Storage() == 3);
  assert(lt.SubMatrix(1, 2, 1, 3).Storage() == 3);
  assert(lt.SubMatrix(1, 3, 1, 2).Storage() == 5);

  assert(lt.SubMatrix(2, 3, 2, 3).Storage() == 3);
  assert(lt.SubMatrix(2, 3, 2, 4).Storage() == 3);
  assert(lt.SubMatrix(2, 4, 2, 3).Storage() == 5);

  UpperTriangularMatrix ut(5);
  assert(ut.Row(1)   .Storage() == 5);
  assert(ut.Row(2)   .Storage() == 4);
  assert(ut.Row(3)   .Storage() == 3);
  assert(ut.Row(4)   .Storage() == 2);
  assert(ut.Row(5)   .Storage() == 1);
  assert(ut.Column(1).Storage() == 1);
  assert(ut.Column(2).Storage() == 2);
  assert(ut.Column(3).Storage() == 3);
  assert(ut.Column(4).Storage() == 4);
  assert(ut.Column(5).Storage() == 5);

  assert(ut.Rows(   1, 2).Storage() == 9);
  assert(ut.Rows(   1, 3).Storage() == 12);
  assert(ut.Rows(   1, 4).Storage() == 14);
  assert(ut.Columns(1, 2).Storage() == 3);
  assert(ut.Columns(1, 3).Storage() == 6);
  assert(ut.Columns(1, 4).Storage() == 10);

  assert(ut.SubMatrix(1, 2, 1, 2).Storage() == 3);
  assert(ut.SubMatrix(1, 2, 1, 3).Storage() == 5);
  assert(ut.SubMatrix(1, 2, 1, 4).Storage() == 7);
  assert(ut.SubMatrix(1, 2, 1, 5).Storage() == 9);
  assert(ut.SubMatrix(1, 3, 1, 2).Storage() == 3);
  assert(ut.SubMatrix(1, 3, 1, 3).Storage() == 6);
  assert(ut.SubMatrix(1, 3, 1, 4).Storage() == 9);
  assert(ut.SubMatrix(1, 3, 1, 5).Storage() == 12);

  assert(ut.SubMatrix(2, 3, 2, 3).Storage() == 3);
  assert(ut.SubMatrix(2, 3, 2, 4).Storage() == 5);
  assert(ut.SubMatrix(2, 3, 2, 5).Storage() == 7);
  assert(ut.SubMatrix(2, 4, 2, 3).Storage() == 3);
  assert(ut.SubMatrix(2, 4, 2, 4).Storage() == 6);
  assert(ut.SubMatrix(2, 4, 2, 5).Storage() == 9);
  assert(ut.SubMatrix(2, 5, 2, 3).Storage() == 3);
  assert(ut.SubMatrix(2, 5, 2, 4).Storage() == 6);
  assert(ut.SubMatrix(2, 5, 2, 5).Storage() == 10);

}

#else
int main(int argc, char *argv[]) { }
#endif
