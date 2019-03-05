#include <string>
#include <vector>
#include "newimage/newimageall.h"
#include "newmat.h"
#include "topup_matrices.h"

namespace TOPUP {

// Global functions

NEWMAT::Matrix MovePar2Matrix(const NEWMAT::ColumnVector&     mp,
                              const NEWIMAGE::volume<float>&  vol) 
{
  if (mp.Nrows() != 6) throw TopupMatrixException("MovePar2Matrix: mp must have 6 elements");

  NEWMAT::ColumnVector tmp(6);
  tmp(1) = mp(4); tmp(2) = mp(5); tmp(3) = mp(6);
  tmp(4) = mp(1); tmp(5) = mp(2); tmp(6) = mp(3);
  
  NEWMAT::ColumnVector cntr(3);
  cntr(1) = ((vol.xsize()-1)*vol.xdim())/2.0;
  cntr(2) = ((vol.ysize()-1)*vol.ydim())/2.0;
  cntr(3) = ((vol.zsize()-1)*vol.zdim())/2.0;

  NEWMAT::Matrix mat(4,4);
  MISCMATHS::construct_rotmat_euler(tmp,6,mat,cntr);

  return(mat);
}

NEWMAT::ColumnVector Matrix2MovePar(const NEWMAT::Matrix&           M,
				    const NEWIMAGE::volume<float>&  vol)
{
  NEWMAT::ColumnVector mp(6); mp = 0.0;
  NEWMAT::ColumnVector rot(3);

  MISCMATHS::rotmat2euler(rot,M);
  mp.Rows(4,6) = rot;

  NEWMAT::Matrix MM = MovePar2Matrix(mp,vol);
  mp.Rows(1,3) = M.SubMatrix(1,3,4,4) - MM.SubMatrix(1,3,4,4);

  return(mp);  
}

} // End namespace TOPUP
