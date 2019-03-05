#include <cassert>

#include "tests.hpp"

#if LIB == ARMAWRAP_API
using namespace armawrap;

int main(int argc, char *argv[]) {

  assert(_tri(1) == 1);
  assert(_tri(2) == 3);
  assert(_tri(3) == 6);
  assert(_tri(4) == 10);
  assert(_tri(5) == 15);
  assert(_tri(6) == 21);
  assert(_tri(7) == 28);

  assert(_ut_diag(0, 0) == 0);
  assert(_ut_diag(1, 0) == 0);
  assert(_ut_diag(2, 0) == 0);
  assert(_ut_diag(2, 1) == 2);
  assert(_ut_diag(3, 0) == 0);
  assert(_ut_diag(3, 1) == 3);
  assert(_ut_diag(3, 2) == 5);
  assert(_ut_diag(4, 0) == 0);
  assert(_ut_diag(4, 1) == 4);
  assert(_ut_diag(4, 2) == 7);
  assert(_ut_diag(4, 3) == 9);
  assert(_ut_diag(5, 0) == 0);
  assert(_ut_diag(5, 1) == 5);
  assert(_ut_diag(5, 2) == 9);
  assert(_ut_diag(5, 3) == 12);
  assert(_ut_diag(5, 4) == 14);
  assert(_ut_diag(6, 0) == 0);
  assert(_ut_diag(6, 1) == 6);
  assert(_ut_diag(6, 2) == 11);
  assert(_ut_diag(6, 3) == 15);
  assert(_ut_diag(6, 4) == 18);
  assert(_ut_diag(6, 5) == 20);
  assert(_ut_diag(7, 0) == 0);
  assert(_ut_diag(7, 1) == 7);
  assert(_ut_diag(7, 2) == 13);
  assert(_ut_diag(7, 3) == 18);
  assert(_ut_diag(7, 4) == 22);
  assert(_ut_diag(7, 5) == 25);
  assert(_ut_diag(7, 6) == 27);

  assert(floor(_ut_diag_root(0, 0)) == 0);
  assert(floor(_ut_diag_root(1, 0)) == 0);
  assert(floor(_ut_diag_root(2, 0)) == 0);
  assert(floor(_ut_diag_root(2, 1)) == 0);
  assert(floor(_ut_diag_root(2, 2)) == 1);
  assert(floor(_ut_diag_root(3, 0)) == 0);
  assert(floor(_ut_diag_root(3, 1)) == 0);
  assert(floor(_ut_diag_root(3, 2)) == 0);
  assert(floor(_ut_diag_root(3, 3)) == 1);
  assert(floor(_ut_diag_root(3, 4)) == 1);
  assert(floor(_ut_diag_root(3, 5)) == 2);
  assert(floor(_ut_diag_root(4, 0)) == 0);
  assert(floor(_ut_diag_root(4, 1)) == 0);
  assert(floor(_ut_diag_root(4, 2)) == 0);
  assert(floor(_ut_diag_root(4, 3)) == 0);
  assert(floor(_ut_diag_root(4, 4)) == 1);
  assert(floor(_ut_diag_root(4, 5)) == 1);
  assert(floor(_ut_diag_root(4, 6)) == 1);
  assert(floor(_ut_diag_root(4, 7)) == 2);
  assert(floor(_ut_diag_root(4, 8)) == 2);
  assert(floor(_ut_diag_root(4, 9)) == 3);

  assert(floor(_ut_diag_root(5, 0))  == 0);
  assert(floor(_ut_diag_root(5, 1))  == 0);
  assert(floor(_ut_diag_root(5, 2))  == 0);
  assert(floor(_ut_diag_root(5, 3))  == 0);
  assert(floor(_ut_diag_root(5, 4))  == 0);
  assert(floor(_ut_diag_root(5, 5))  == 1);
  assert(floor(_ut_diag_root(5, 6))  == 1);
  assert(floor(_ut_diag_root(5, 7))  == 1);
  assert(floor(_ut_diag_root(5, 8))  == 1);
  assert(floor(_ut_diag_root(5, 9))  == 2);
  assert(floor(_ut_diag_root(5, 10)) == 2);
  assert(floor(_ut_diag_root(5, 11)) == 2);
  assert(floor(_ut_diag_root(5, 12)) == 3);
  assert(floor(_ut_diag_root(5, 13)) == 3);
  assert(floor(_ut_diag_root(5, 14)) == 4);

  assert(floor(_ut_diag_root(6, 0))  == 0);
  assert(floor(_ut_diag_root(6, 1))  == 0);
  assert(floor(_ut_diag_root(6, 2))  == 0);
  assert(floor(_ut_diag_root(6, 3))  == 0);
  assert(floor(_ut_diag_root(6, 4))  == 0);
  assert(floor(_ut_diag_root(6, 5))  == 0);
  assert(floor(_ut_diag_root(6, 6))  == 1);
  assert(floor(_ut_diag_root(6, 7))  == 1);
  assert(floor(_ut_diag_root(6, 8))  == 1);
  assert(floor(_ut_diag_root(6, 9))  == 1);
  assert(floor(_ut_diag_root(6, 10)) == 1);
  assert(floor(_ut_diag_root(6, 11)) == 2);
  assert(floor(_ut_diag_root(6, 12)) == 2);
  assert(floor(_ut_diag_root(6, 13)) == 2);
  assert(floor(_ut_diag_root(6, 14)) == 2);
  assert(floor(_ut_diag_root(6, 15)) == 3);
  assert(floor(_ut_diag_root(6, 16)) == 3);
  assert(floor(_ut_diag_root(6, 17)) == 3);
  assert(floor(_ut_diag_root(6, 18)) == 4);
  assert(floor(_ut_diag_root(6, 19)) == 4);
  assert(floor(_ut_diag_root(6, 20)) == 5);

  assert(_lt_diag(0)  == 0);
  assert(_lt_diag(1)  == 2);
  assert(_lt_diag(2)  == 5);
  assert(_lt_diag(3)  == 9);
  assert(_lt_diag(4)  == 14);
  assert(_lt_diag(5)  == 20);
  assert(_lt_diag(6)  == 27);
  assert(_lt_diag(7)  == 35);
  assert(_lt_diag(8)  == 44);
  assert(_lt_diag(9)  == 54);
  assert(_lt_diag(10) == 65);
  assert(_lt_diag(11) == 77);
  assert(_lt_diag(12) == 90);
  assert(_lt_diag(13) == 104);
  assert(_lt_diag(14) == 119);
  assert(_lt_diag(15) == 135);
  assert(_lt_diag(16) == 152);
  assert(_lt_diag(17) == 170);
  assert(_lt_diag(18) == 189);
  assert(_lt_diag(19) == 209);
  assert(_lt_diag(20) == 230);

  assert(ceil(_lt_diag_root(0))  == 0);
  assert(ceil(_lt_diag_root(1))  == 1);
  assert(ceil(_lt_diag_root(2))  == 1);
  assert(ceil(_lt_diag_root(3))  == 2);
  assert(ceil(_lt_diag_root(4))  == 2);
  assert(ceil(_lt_diag_root(5))  == 2);
  assert(ceil(_lt_diag_root(6))  == 3);
  assert(ceil(_lt_diag_root(7))  == 3);
  assert(ceil(_lt_diag_root(8))  == 3);
  assert(ceil(_lt_diag_root(9))  == 3);
  assert(ceil(_lt_diag_root(10)) == 4);
  assert(ceil(_lt_diag_root(11)) == 4);
  assert(ceil(_lt_diag_root(12)) == 4);
  assert(ceil(_lt_diag_root(13)) == 4);
  assert(ceil(_lt_diag_root(14)) == 4);
  assert(ceil(_lt_diag_root(15)) == 5);
  assert(ceil(_lt_diag_root(16)) == 5);
  assert(ceil(_lt_diag_root(17)) == 5);
  assert(ceil(_lt_diag_root(18)) == 5);
  assert(ceil(_lt_diag_root(19)) == 5);
  assert(ceil(_lt_diag_root(20)) == 5);
  assert(ceil(_lt_diag_root(21)) == 6);
}
#else
int main(int argc, char *argv[]) { }
#endif
