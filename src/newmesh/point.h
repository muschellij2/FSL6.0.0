/*  Copyright (C) 1999-2004 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Oxford
    University Innovation ("OUI"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    Innovation@innovation.ox.ac.uk quoting reference DE/9564. */

#ifndef _point
#define _point

#include <iostream>
#include <cmath>
#include "newmat.h"

using namespace std;
using namespace NEWMAT;

namespace NEWMESH{

  class Pt {

  public:

  double X;
  double Y; 
  double Z;
  
  Pt() : X(0), Y(0), Z(0){};
  Pt(double x, double y, double z) : X(x), Y(y), Z(z){};
  Pt (const Pt& p) : X(p.X), Y(p.Y), Z(p.Z){};


  inline Pt operator =(const Pt& p)
    {
	X = p.X;
	Y = p.Y;
	Z = p.Z;
	return *this;
    }
    
    
  
    inline void operator+=(const Pt p)      
    {
      X=X+p.X;
      Y=Y+p.Y;
      Z=Z+p.Z;
    }

    inline void operator*=(const double d)    
    {
      X*=d;
      Y*=d;
      Z*=d;
    }

  inline void operator/=(const double d) 
    {
      if (d!=0)
	{
	  X/=d;
	  Y/=d;
	  Z/=d;
	}
      else cerr << "division by zero" << endl;
    }

   const inline double norm() const /// transferred from Vec by emma
   {
      return (sqrt(X*X + Y*Y + Z*Z));
   }

  void normalize(){      /// transferred from Vec by emma
    double n = norm();
    if (n!=0){
    X/=n;
    Y/=n;
    Z/=n;}
  }

  inline bool operator==(const Pt &p) const
    {
      return((fabs(X-p.X)<1e-8) && (fabs(Y-p.Y)<1e-8) && (fabs(Z-p.Z)<1e-8));
    }

};

// Pt operators added by emma

  const double operator|(const Pt &p1, const Pt &p2);    // dot product
  const Pt operator*(const Pt &p1, const Pt &p2);        // cross product
  const Pt operator*(const Pt &v, const double &d);
  const Pt operator/(const Pt &v, const double &d);   //NB all these were const in original meshclass, but I don't think they need/should be
  const Pt operator -(const Pt&p1, const Pt &p2);
  const Pt operator +(const Pt&p1, const Pt &p2);
  const Pt operator*(const Matrix &M,const Pt &p1); 
  const Pt operator*(const Pt &p1, const Matrix &M); 
}

#endif
