/*  resampler.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  

    Some sections of code inspired by A. Petrovic.
*/
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
/*  Adaptive Barycentric code used with the permission of Tim Coulson under the below licence:
 *  Copyright (C) 2014  Washington University School of Medicine
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be included
 *  in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 *  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 *  CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 *  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 *  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.*/

/*  RESAMPLING CLASS, PERFORMS DATA INTERPOLATION ON MESH SURFACES */
/*  ASSUMED THAT MESHES ARE SPHERES AS NEIGHBOURHOODS ARE DEFINED USING RELATIONS CLASS */
/* IF RELATIONS WAS MADE MORE GENERAL (use aprox nearest neighbours class and geodesic distances) & DISTANCES WERE CALCULATED IN RESAMPLER AS GEODESICS */
/* THEN THIS CLASS COULD RESAMPLE TO ANY MESH */
/* ALSO CONTAINS RBF CLASS USED FOR TPS MESH INTERPOLATION, & FUNCTIONS FOR WARP INTERPOLATIONS */

#ifndef resampler_h
#define resampler_h

#include <fstream>


#include "utils/log.h"
#include "Relations.h"

#include "miscmaths/bfmatrix.h"
#include <boost/shared_ptr.hpp>
//#include <time.h>


// RESAMPLING AND SMOOTHING FUNCTIONS /////////////////

using namespace MISCMATHS;
using namespace std;
using namespace NEWMESH;


#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
/* should define sruct and mem just once somewhere  - in a class? */

namespace NEWMESH {
 
  //RBF interpolation methods for Thin Plate Spline interpolation see Boostein Principal warps: Thin plate splines ...
  //-------------------------
  class RBF{
    
  private:
   
   Matrix        cent;
   ColumnVector  f;
   ColumnVector coeff ;

   ReturnMatrix fitit(Matrix& cent, ColumnVector& f);
   ReturnMatrix direct(Matrix& cent);
   double eval_direct(Matrix& cent, ColumnVector coeff, double uX, double uY, double uZ);
   
 public:
   
   RBF(){};
   ~RBF(){};
   
   void set_point(double X, double Y, double Z, double fn);
   double interpolate(double X, double Y, double Z,bool estimate_coeff=true);
   inline void clean() {cent.CleanUp(); f.CleanUp();}
  };
  
 
  class resampler{  
   
    string _method;

  public:
   
   
    resampler(){_method="NN";} // options for data resampling: Nearest Neighbour (NN), linear, gaussian, sum or mean (of points within kernel range), barycentric and adaptive barycentric 
    ~resampler(){}
   
   RELATIONS  initialize_relations(const double &,NEWMESH::newmesh & ,const  NEWMESH::newmesh &); 
   
   void set_method(const string &);
   string get_method(){return _method;};  


   double calc_weight(const double &,const double &); /// weights use euclidean distances at this time.
   double guess_angular_spacing(const int &);
   //// UTILITIES FOR GENERIC BFMATRIX DATA ////
   void resampledata(NEWMESH::newmesh& ,const NEWMESH::newmesh&, boost::shared_ptr<NEWMESH::newmesh> &, boost::shared_ptr<BFMatrix> &, const double &,boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >());
   void resampledata(NEWMESH::newmesh& ,const NEWMESH::newmesh&,  boost::shared_ptr<BFMatrix> &, const double &, boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >()); 
   
   void nearest_neighbour_interpolation(NEWMESH::newmesh &,const  NEWMESH::newmesh &, boost::shared_ptr<BFMatrix> &, boost::shared_ptr<NEWMESH::newmesh> &, boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >());
   void barycentric_interpolation(NEWMESH::newmesh &,const  NEWMESH::newmesh  &, boost::shared_ptr<BFMatrix> &,  boost::shared_ptr<NEWMESH::newmesh> &, boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >()); // also adaptive barycentric

   vector<std::map<int,double> > get_all_barycentric_weights(NEWMESH::newmesh &MESHA, NEWMESH::newmesh &MESHB, boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >());
   void reverse_barycentric_weights(vector<std::map<int,double> > &, const int &, const int &);

   vector<std::map<int,double> > get_all_adap_barycentric_weights(NEWMESH::newmesh &IN, NEWMESH::newmesh &ref, boost::shared_ptr<NEWMESH::newmesh> &, boost::shared_ptr<RELATIONS > _rel=boost::shared_ptr<RELATIONS >());
   map<int,double> get_barycentric_weight_for_ind(const int &,const double &,const Pt &,const NEWMESH::newmesh &,const RELATIONS &);
   map<int,double> get_barycentric_weight_for_ind(int,double , NEWMESH::newmesh &, NEWMESH::newmesh &,RELATIONS &);

   //upsample and downsample differently coded for speed.
   void upsample_w_interpolation(NEWMESH::newmesh &,const NEWMESH::newmesh &, const double &,  boost::shared_ptr<BFMatrix> &, boost::shared_ptr<NEWMESH::newmesh> &, boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >());
   void downsample_w_interpolation(NEWMESH::newmesh &,const NEWMESH::newmesh &, const double &,  boost::shared_ptr<BFMatrix> &, boost::shared_ptr<NEWMESH::newmesh> &, boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >());
   void return_all_triangles(vector<vector<int> > &, vector<vector<int> > &, NEWMESH::newmesh &,const NEWMESH::newmesh &);
   void smooth_data(const double & ,  boost::shared_ptr<BFMatrix> &, NEWMESH::newmesh &);   

   //// UTILITIES FOR FULL MATRIX DATA ////
   void resampledata(NEWMESH::newmesh& ,const NEWMESH::newmesh&, boost::shared_ptr<NEWMESH::newmesh>&, Matrix &, const double &, boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >());
   void resampledata(NEWMESH::newmesh& ,const NEWMESH::newmesh&, Matrix &, const double &,boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >()); 
  
   //// UTILITIES FOR SPARSE MATRIX DATA ////
   void resampledata(NEWMESH::newmesh&, const NEWMESH::newmesh&, boost::shared_ptr<NEWMESH::newmesh>&, SpMat<double> &, const double &, boost::shared_ptr<RELATIONS > _rel =boost::shared_ptr<RELATIONS >());
   
   // scalar data 
   
   void resample_scalar(NEWMESH::newmesh &,const NEWMESH::newmesh &,const double &, boost::shared_ptr<NEWMESH::newmesh> EXCL= boost::shared_ptr<NEWMESH::newmesh>());  
    
 };
 

  //// RESAMPLE ANATOMICAL MESH ACCORDING TO SPHERICAL RESAMPLING
  newmesh mesh_resample(const NEWMESH::newmesh &,NEWMESH::newmesh,NEWMESH::newmesh&,bool _adapbary=false, boost::shared_ptr<RELATIONS > _rel=boost::shared_ptr<RELATIONS >());

 ///////////// PROJECTING WARP TO NEW MESHES /////////////////////////
  void barycentric_mesh_interpolation(NEWMESH::newmesh &,NEWMESH::newmesh &,NEWMESH::newmesh &,boost::shared_ptr<RELATIONS > _rel=boost::shared_ptr<RELATIONS >());
  void upsample_transform_RBF(NEWMESH::newmesh &,NEWMESH::newmesh &,NEWMESH::newmesh &, const double&); // upsample registration output to high res sphere once again
  void upsample_transform_RBF(NEWMESH::newmesh &,NEWMESH::newmesh &,NEWMESH::newmesh &, const RELATIONS &); // upsample registration output to high res sphere once again
  
  void project(NEWMESH::newmesh &, NEWMESH::newmesh&,const NEWMESH::newmesh&,const double&);  // project transformation using a locally affine transformation
  Matrix affine_transform(NEWMESH::newmesh &, NEWMESH::newmesh , NEWMESH::newmesh );
  void affine_transform(NEWMESH::newmesh &,Matrix);

  Pt project(const vector<int> &,const NEWMESH::Pt&,const NEWMESH::newmesh &,const NEWMESH::newmesh &);
  Matrix return_affine_transform(const vector<int>&,const NEWMESH::newmesh &,const NEWMESH::newmesh &);
  Matrix return_affine_transform(const vector<Pt>&,const vector<Pt>&);


  ///////////////// MISCELLANEOUS //////////////
  double Calculate_MVD(const NEWMESH::newmesh &); /// mean vertex displacement function - ONLY FOR REGULAR MESH - this cheats and uses the first vertex spacings only
  double Calculate_MaxVD(const NEWMESH::newmesh &); // calculates max vertex displacement, searches whole mesh, uses euclidean spacings
  double Calculate_MaxVD(const NEWMESH::newmesh &, ColumnVector &); /// as above
  
    
  double computevertexArea(const int &,const NEWMESH::newmesh &); // averages adjoining face areas for each vertex
  
  double computeArea(const NEWMESH::Pt &v0, const NEWMESH::Pt &v1, const NEWMESH::Pt &v2); // calculates area of triangular mesh face
  
  ////  these functions output barycentric weights in different formats and need combining or reducing
  double barycentric(const NEWMESH::Pt&, const NEWMESH::Pt&,const  NEWMESH::Pt&,const  NEWMESH::Pt&,const double &,const double &,const double &); /// inspired by spherical demons  
  Pt barycentric(const NEWMESH::Pt&, const NEWMESH::Pt&,const  NEWMESH::Pt&,const  NEWMESH::Pt&,const Pt &,const Pt &,const Pt &); /// inspired by spherical demons  

  SpMat<double> barycentric(const NEWMESH::Pt&, const NEWMESH::Pt&, const NEWMESH::Pt&, const NEWMESH::Pt&, const int &,const int &,const int &,const SpMat<double> &);   
  void barycentric2(const NEWMESH::Pt&, const NEWMESH::Pt&, const NEWMESH::Pt&, const NEWMESH::Pt&, const int &,const int &,const int &,const BFMatrix &, vector<double> &);
  std::map<int,double> get_barycentric_weights(const NEWMESH::Pt&, const NEWMESH::Pt&, const NEWMESH::Pt&, const NEWMESH::Pt&, const int &,const int &,const int &);
} 

#endif
