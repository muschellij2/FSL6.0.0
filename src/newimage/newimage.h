/*  Templated image storage class

    Mark Jenkinson, FMRIB Image Analysis Group

    Copyright (C) 2000 University of Oxford  */

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

#if !defined(__newimage_h)
#define __newimage_h

#define volume4D volume

#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <stdexcept>
#include "newmatap.h"
#include "positerators.h"
#include "NewNifti/NewNifti.h"
#include "miscmaths/miscmaths.h"
#include "miscmaths/kernel.h"
#include "miscmaths/splinterpolator.h"
#include <boost/type_traits/is_integral.hpp>


using namespace NEWMAT;
using namespace MISCMATHS;
using namespace std;

namespace NEWIMAGE {


  #define SUBSET -1

  const bool USEASSERT=false;  // used inside imthrow to determine type of action

  void imthrow(const string& msg, int nierrnum, const bool quiet=false);

  enum extrapolation { zeropad, constpad, extraslice, mirror, periodic,
                       boundsassert, boundsexception, userextrapolation };

  enum interpolation { nearestneighbour, trilinear, sinc, userkernel,
		       userinterpolation, spline };

  enum threshtype { inclusive , exclusive };

#define FSL_TYPE_ANALYZE         10
#define FSL_TYPE_NIFTI            1
#define FSL_TYPE_NIFTI2           2
#define FSL_TYPE_NIFTI_PAIR      11
#define FSL_TYPE_NIFTI2_PAIR     12
#define FSL_TYPE_ANALYZE_GZ     100
#define FSL_TYPE_NIFTI_GZ       101
#define FSL_TYPE_NIFTI2_GZ      102
#define FSL_TYPE_NIFTI_PAIR_GZ  111
#define FSL_TYPE_NIFTI2_PAIR_GZ 112

#define FSL_RADIOLOGICAL        -1
#define FSL_NEUROLOGICAL         1
#define FSL_INCONSISTENT         0
#define FSL_ZERODET           -101

template<class T> class ShadowVolume;
template<class T> class volume;

template <class T>
int read_volumeROI(volume<T>& target, const string& filename,
		   short& dtype,
		   int64_t x0, int64_t y0, int64_t z0, int64_t t0, int64_t d50, int64_t d60, int64_t d70,
		   int64_t x1, int64_t y1, int64_t z1, int64_t t1, int64_t d51, int64_t d61, int64_t d71,
		   const bool swap2radiological=true, const bool readAs4D=false);

#pragma interface
  template <class T>
  class volume {
  private:
    T* Data;
    T* DataEnd;
    virtual int isShadow() { return 0; }
    mutable bool data_owner;
    mutable double maskDelimiter;
    int64_t nElements;

    int64_t ColumnsX;
    int64_t RowsY;
    int64_t SlicesZ;
    int64_t dim4;
    int64_t dim5;
    int64_t dim6;
    int64_t dim7;
    vector<int64_t> originalSizes;

    float Xdim;
    float Ydim;
    float Zdim;
    float p_TR;
    float pxdim5;
    float pxdim6;
    float pxdim7;

    // the following matrices map voxel coords to mm coords (nifti conventions)
    Matrix StandardSpaceCoordMat;
    Matrix RigidBodyCoordMat;
    int StandardSpaceTypeCode;
    int RigidBodyTypeCode;

    mutable int IntentCode;
    mutable float IntentParam1;
    mutable float IntentParam2;
    mutable float IntentParam3;

    mutable int SliceOrderingCode;

   int64_t no_voxels;

    mutable SPLINTERPOLATOR::Splinterpolator<T> splint;

    mutable kernel interpkernel;
    mutable extrapolation p_extrapmethod;

    mutable interpolation p_interpmethod;
    mutable int splineorder;
    mutable T (*p_userextrap)(const volume<T>& , int64_t, int64_t, int64_t);
    mutable float (*p_userinterp)(const volume<T>& , float, float, float);
    mutable T padvalue; // was p_padval in volume4D

    mutable T extrapval;                 // the reference target for all extrapolations
    mutable std::vector<bool> ep_valid;  // Indicates if extrapolation can produce "valid" values.

    mutable float displayMaximum;
    mutable float displayMinimum;
    mutable char auxFile[24];

    // Internal functions
    inline T* basicptr(int64_t x, int64_t y, int64_t z) {
        return (nsfbegin() + (z*RowsY + y)*ColumnsX + x);
    }
    inline T* basicptr(int64_t x, int64_t y, int64_t z) const
      { return (Data + (z*RowsY + y)*ColumnsX + x); }
    int initialize(int64_t xsize, int64_t ysize, int64_t zsize, T *d, bool d_owner); //3D
    int initialize(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize, T *d, bool d_owner); //4D
    int initialize(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize, int64_t d5, int64_t d6, int64_t d7, T *d, bool d_owner); //Master 7D
    void setdefaultproperties();
    void enforcelimits(std::vector<int>& lims) const;
    const T& extrapolate(int64_t x, int64_t y, int64_t z) const;
    float kernelinterpolation(const float x, const float y,
			      const float z) const;
    float splineinterpolate(float x, float y, float z) const;
    float spline_interp1partial(float x, float y, float z, int dir, float *deriv) const;
    float spline_interp3partial(float x, float y, float z, float *dfdx, float *dfdy, float *dfdz) const;

    template <class S, class D> friend
    void copybasicproperties(const volume<S>& source, volume<D>& dest);

    template <class S> friend
    int read_volumeROI(volume<S>& target, const string& filename,
		     short& dtype,
		     int64_t x0, int64_t y0, int64_t z0, int64_t t0, int64_t d50, int64_t d60, int64_t d70,
		     int64_t x1, int64_t y1, int64_t z1, int64_t t1, int64_t d51, int64_t d61, int64_t d71,
		     const bool swap2radiological, const bool readAs4D);

    template <class S> friend
      int read_volume_hdr_only(volume<S>&, const std::string&);

    void basic_swapdimensions(int dim1, int dim2, int dim3, bool keepLRorder);

#ifdef EXPOSE_TREACHEROUS
  public:
#endif
    // sampling_mat should now be avoided - use newimagevox2mm_mat instead
    std::vector<NiftiExtension> extensions;
    bool RadiologicalFile;
    Matrix sampling_mat() const;
    void set_sform(int sform_code, const Matrix& snewmat);
    void set_qform(int qform_code, const Matrix& qnewmat);
    int  left_right_order() const; // see also newimagevox2mm_mat()
    void swapLRorder();
    void setLRorder(int LRorder);
    void makeradiological();
    void makeneurological();
    Matrix sform_mat() const { return StandardSpaceCoordMat; }
    int sform_code() const { return StandardSpaceTypeCode; }
    Matrix qform_mat() const { return RigidBodyCoordMat; }
    int qform_code() const { return RigidBodyTypeCode; }
    void set_data_owner(bool d_owner) const { data_owner=d_owner; }
    public:
    typedef T* nonsafe_fast_iterator;
    inline nonsafe_fast_iterator nsfbegin()
      {       this->negateSpline(); return Data; }
    inline nonsafe_fast_iterator nsfend()
    { return DataEnd; }
   public:
    // CONSTRUCTORS AND DESTRUCTORS (including copy, = and reinitialize)
    volume();
    volume(const volume<T>& source, const bool copyingData=true);
    volume(int64_t xsize, int64_t ysize, int64_t zsize);
    volume(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize);
    volume(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize, int64_t d5, int64_t d6, int64_t d7);
    ~volume();
    void destroy();
    virtual const volume<T>& operator=(const volume<T>& source);
    virtual const volume<T>& operator=(const ShadowVolume<T>& source);

    void reinitialize(const volume<T>& source, const bool copyingData=true);
    int reinitialize(int64_t xsize, int64_t ysize, int64_t zsize);
    int reinitialize(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize,T *d=NULL, bool d_owner=true);
    int reinitialize(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize, int64_t d5, int64_t d6, int64_t d7);
    void copyproperties(const volume<T>& source);
    int copydata(const volume<T>& source, const bool isOwner=true);

    // BASIC PROPERTIES
    inline int64_t xsize() const { return ColumnsX; }
    inline int64_t ysize() const { return RowsY; }
    inline int64_t zsize() const { return SlicesZ; }
    inline int64_t tsize() const { return dim4; }
    inline int64_t size1() const { return xsize(); }
    inline int64_t size2() const { return ysize(); }
    inline int64_t size3() const { return zsize(); }
    inline int64_t size4() const { return tsize(); }
    inline int64_t size5() const { return dim5; }
    inline int64_t size6() const { return dim6; }
    inline int64_t size7() const { return dim7; }
    int64_t size(int n) const;

    int dimensionality() const;
    void throwsIfNot3D() const { if (this->dimensionality()>3) imthrow("3D only method called by higher-dimensional volume.",75); }

    inline int64_t maxx() const { return xsize()-1; }
    inline int64_t maxy() const { return ysize()-1; }
    inline int64_t maxz() const { return zsize()-1; }
    inline int64_t maxt() const { return tsize()-1; }
    inline int64_t minx() const { return 0; }
    inline int64_t miny() const { return 0; }
    inline int64_t minz() const { return 0; }
    inline int64_t mint() const { return 0; }

    inline float xdim() const { return Xdim; }
    inline float ydim() const { return Ydim; }
    inline float zdim() const { return Zdim; }
    inline float tdim() const { return p_TR; }
    inline float TR() const { return p_TR; }
    inline float pixdim1() const { return Xdim; }
    inline float pixdim2() const { return Ydim; }
    inline float pixdim3() const { return Zdim; }
    inline float pixdim4() const { return p_TR; }
    inline float pixdim5() const { return pxdim5; }
    inline float pixdim6() const { return pxdim6; }
    inline float pixdim7() const { return pxdim7; }

    inline float getDisplayMaximum() const { return displayMaximum; }
    inline float getDisplayMinimum() const { return displayMinimum; }
    inline string getAuxFile() const { return string(auxFile); }

    void setxdim(float x) { Xdim = fabs(x); }
    void setydim(float y) { Ydim = fabs(y); }
    void setzdim(float z) { Zdim = fabs(z); }
    void settdim(float tr) { p_TR = fabs(tr); }
    void setd5dim(float d5) { pxdim5 = d5; }
    void setd6dim(float d6) { pxdim6 = d6; }
    void setd7dim(float d7) { pxdim7 = d7; }
    void setTR(float tr)   { settdim(tr); }
    void setdims(float x, float y, float z)
      { setxdim(x); setydim(y); setzdim(z); }
    void setdims(float x, float y, float z, float tr)
      { setxdim(x); setydim(y); setzdim(z); settdim(tr); }
    void setdims(float x, float y, float z, float tr, float d5, float d6, float d7)
      { setxdim(x); setydim(y); setzdim(z); settdim(tr); setd5dim(d5); setd6dim(d6); setd7dim(d7); }
    void setDisplayMaximumMinimum(const float maximum, const float minimum) const {  displayMaximum=maximum; displayMinimum=minimum; }
    void setDisplayMaximum(const float maximum) const { setDisplayMaximumMinimum(maximum,displayMinimum); }
    void setDisplayMinimum(const float minimum) const { setDisplayMaximumMinimum(displayMaximum,minimum); }
    void setAuxFile(const string fileName) { strncpy(auxFile,fileName.c_str(),24); }
    int64_t nvoxels() const { return no_voxels; }
    int64_t ntimepoints() const { return tsize(); }
    int64_t totalElements() const { return nElements; }

    vector<int64_t> ptrToCoord(size_t offset) const;

    // MATRIX <-> VOLUME CONVERSIONS
    ReturnMatrix vec(const volume<T>& mask) const;
    ReturnMatrix vec() const;
    void insert_vec(const ColumnVector& pvec, const volume<T>& pmask);
    void insert_vec(const ColumnVector& pvec);
    vector<int64_t> labelToCoord(const int64_t label) const;
    ReturnMatrix matrix(const volume<T>& mask, bool using_mask=true) const;
    ReturnMatrix matrix(const volume<T>& mask, vector<int64_t>& voxelLabels) const;
    ReturnMatrix matrix() const;
    void setmatrix(const Matrix& newmatrix, const volume<T>& mask,
		   const T pad=0, bool using_mask=true);
    void setmatrix(const Matrix& newmatrix);
    volume<int> vol2matrixkey(const volume<T>& mask); //returns a volume with numbers in relating to matrix colnumbers
    ReturnMatrix matrix2volkey(volume<T>& mask);

    // SECONDARY PROPERTIES
    // maps *NEWIMAGE* voxel coordinates to mm (consistent with FSLView mm)
    // NB: do not try to determine left-right order from this matrix
    // sampling_mat should now be avoided - use newimagevox2mm_mat instead
    Matrix newimagevox2mm_mat() const;
    Matrix niftivox2newimagevox_mat() const;

    int intent_code() const { return IntentCode; }
    float intent_param(int n) const;
    void set_intent(int intent_code, float p1, float p2, float p3) const;

    inline double maskThreshold() const { return maskDelimiter; }
    T min() const;
    T max() const;
    double sum() const;
    double sumsquares() const;
    double mean() const { return sum()/(Max(1.0,(double)totalElements()));}
    double variance() const { double n(totalElements());
		return (n/(n-1))*(sumsquares()/n - mean()*mean()); }
    double stddev() const { return sqrt(variance()); }
    T robustmin() const;
    T robustmax() const;
    vector<T> robustlimits(const volume<T>& mask) const;
    ColumnVector principleaxis(int n) const;
    Matrix principleaxes_mat() const;
    T percentile(float pvalue) const;  // argument in range [0.0 , 1.0]
    ColumnVector histogram(int nbins) const;
    ColumnVector histogram(int nbins, T minval, T maxval) const;
    ColumnVector cog(const string& coordtype="voxel") const;

    // SECONDARY PROPERTIES (using mask)
    T min(const volume<T>& mask) const;
    T max(const volume<T>& mask) const;
    T min(const volume<T>& mask, vector<int64_t>& coords) const;
    T max(const volume<T>& mask, vector<int64_t>& coords) const;
    double sum(const volume<T>& mask) const;
    double sumsquares(const volume<T>& mask) const;
    double mean(const volume<T>& mask) const;
    double variance(const volume<T>& mask) const;
    double stddev(const volume<T>& mask) const { return sqrt(variance(mask)); }
    T robustmin(const volume<T>& mask) const;
    T robustmax(const volume<T>& mask) const;
    T percentile(float pvalue, const volume<T>& mask) const;  // arg in [0,1]
    ColumnVector histogram(int nbins, const volume<T>& mask) const;
    ColumnVector histogram(int nbins, T minval, T maxval, const volume<T>& mask)
      const;

    // DATA ACCESS FUNCTIONS (iterators)

    typedef const T* fast_const_iterator;

    inline fast_const_iterator fbegin() const { return Data; }
    inline fast_const_iterator fend() const { return DataEnd; }
    inline fast_const_iterator fbegin(const int64_t offset) const { return Data + offset*nvoxels(); }
    inline fast_const_iterator fend(const int64_t offset) const { return Data + (offset+1)*nvoxels(); }

    // BASIC DATA ACCESS FUNCTIONS
    template<class A, class B, class C>
    inline bool in_bounds( A x, B y, C z) const {
      if ( boost::is_integral<A>::value && boost::is_integral<B>::value && boost::is_integral<C>::value )
	return ( (x>=0) && (y>=0) && (z>=0) && (x<ColumnsX) && (y<RowsY) && (z<SlicesZ) );
      return ( (floor(x)>=0) && (floor(y)>=0) && (floor(z)>=0) && (ceil(x)<ColumnsX) && (ceil(y)<RowsY) && (ceil(z)<SlicesZ) );
    }

    inline bool in_bounds(int64_t t) const
      { return ( (t>=0) && (t<this->tsize()) ); }

    template<class A, class B, class C, class D>
    inline bool in_bounds(A x, B y, C z, D t) const
      { return ( this->in_bounds(t) && this->in_bounds(x,y,z) ); }

    inline bool in_bounds(int64_t x, int64_t y, int64_t z, int64_t t, int64_t d5, int64_t d6=0, int64_t d7=0) const {
      return (this->in_bounds(t) && this->in_bounds(x,y,z) && (d5>=0) && (d6>=0) && (d7>=0) && (d5<dim5) && (d6<dim6) && (d7<dim7));}

    bool in_extraslice_bounds(float x, float y, float z) const
    {
      int64_t ix=((int64_t) floor(x));
      int64_t iy=((int64_t) floor(y));
      int64_t iz=((int64_t) floor(z));
      return((ix>=-1) && (iy>=-1) && (iz>=-1) && (ix<ColumnsX) && (iy<RowsY) && (iz<SlicesZ));
    }
    inline bool valid(int64_t x, int64_t y, int64_t z) const
    {
      return((ep_valid[0] || (x>=0 && x<ColumnsX)) && (ep_valid[1] || (y>=0 && y<RowsY)) && (ep_valid[2] || (z>=0 && z<SlicesZ)));
    }
    bool valid(float x, float y, float z, double tol=1e-8) const
    {
      // int ix=((int) floor(x));
      // int iy=((int) floor(y));
      // int iz=((int) floor(z));
      return((ep_valid[0] || (x+tol >= 0.0 && x <= ColumnsX-1+tol)) &&
             (ep_valid[1] || (y+tol >= 0.0 && y <= RowsY-1+tol)) &&
             (ep_valid[2] || (z+tol >= 0.0 && z <= SlicesZ-1+tol)));
      // return((ep_valid[0] || (ix>=0 && (ix+1)<ColumnsX)) && (ep_valid[1] || (iy>=0 && (iy+1)<RowsY)) && (ep_valid[2] || (iz>=0 && (iz+1)<SlicesZ)));
    }
    inline T& operator()(int64_t x, int64_t y, int64_t z) {
	    if (in_bounds(x,y,z)) return *(basicptr(x,y,z));
	    else                  return const_cast<T& > (extrapolate(x,y,z));
    }
    inline const T& operator()(int64_t x, int64_t y, int64_t z) const {
      if (in_bounds(x,y,z)) return *(basicptr(x,y,z));
	        else                  return extrapolate(x,y,z);
    }

   inline T& operator()(int64_t x, int64_t y, int64_t z, int64_t t)
   { if (!in_bounds(t)) imthrow("Out of Bounds (time index)",5);
     else if (!in_bounds(x,y,z)) return const_cast<T& > (operator[](t).extrapolate(x,y,z));
	return *(Data + ((zsize()*t + z)*ysize() + y)*xsize() + x); }
    inline const T& operator()(int64_t x, int64_t y, int64_t z, int64_t t) const
    { if (!in_bounds(t)) imthrow("Out of Bounds (time index)",5);
      else if (!in_bounds(x,y,z)) return const_cast<T& > (operator[](t).extrapolate(x,y,z));
        return *(Data + ((zsize()*t + z)*ysize() + y)*xsize() + x); }

    inline T& operator()(int64_t x, int64_t y, int64_t z, int64_t t, int64_t d5, int64_t d6=0, int64_t d7=0)
    { if (!in_bounds(x,y,z,t,d5,d6,d7)) imthrow("Out of Bounds (7D index)",5);
	return *(Data + ((((((size6()*d7+d6)*size5()+d5)*tsize()+t)*zsize()*t+z)*ysize()+y)*xsize()+ x)); }
    inline const T& operator()(int64_t x, int64_t y, int64_t z, int64_t t, int64_t d5, int64_t d6=0, int64_t d7=0) const
    { if (!in_bounds(x,y,z,t,d5,d6,d7)) imthrow("Out of Bounds (7D index)",5);
        return *(Data + ((((((size6()*d7+d6)*size5()+d5)*tsize()+t)*zsize()*t+z)*ysize()+y)*xsize()+ x)); }



    float interpolate(float x, float y, float z) const;
    float interpolate(float x, float y, float z, bool *ep) const;
    float interp1partial(// Input
                         float x, float y, float z,  // Co-ordinates to get value for
                         int     dir,                   // Direction for partial, 0->x, 1->y, 2->z
                         // Output
                         float  *pderiv                // Derivative returned here
                         ) const;
    float interp3partial(// Input
                         float x, float y, float z,             // Co-ordinate to get value for
                         // Output
                         float *dfdx, float *dfdy, float *dfdz  // Partials
                         ) const;

    inline T& value(int64_t x, int64_t y, int64_t z)
      { return *(basicptr(x,y,z)); }
    inline const T& value(int64_t x, int64_t y, int64_t z) const
      { return *(basicptr(x,y,z)); }

    inline T& value(int64_t x, int64_t y, int64_t z, int64_t t)
    { return *(Data + ((zsize()*t + z)*ysize() + y)*xsize() + x); }
    inline const T& value(int64_t x, int64_t y, int64_t z, int64_t t) const
    { return *(Data + ((zsize()*t + z)*ysize() + y)*xsize() + x); }
    inline T& value(int64_t x, int64_t y, int64_t z, int64_t t, int64_t d5, int64_t d6=0, int64_t d7=0)
    { return *(Data + ((((((size6()*d7+d6)*size5()+d5)*tsize()+t)*zsize()+z)*ysize()+y)*xsize()+ x)); }
    inline const T& value(int64_t x, int64_t y, int64_t z, int64_t t, int64_t d5, int64_t d6=0, int64_t d7=0) const
    { return *(Data + ((((((size6()*d7+d6)*size5()+d5)*tsize()+t)*zsize()+z)*ysize()+y)*xsize()+ x)); }

    float interpolatevalue(float x, float y, float z) const;
    ColumnVector ExtractRow(int j, int k) const;
    ColumnVector ExtractColumn(int i, int k) const;
    void SetRow(int j, int k, const ColumnVector& row);
    void SetColumn(int j, int k, const ColumnVector& col);

    // ROI and higher dimensional accesses
    void addvolume(const volume<T>& source);
    void deletevolume(int64_t t);
    void clear();  // deletes all volumes
    ReturnMatrix voxelts(int64_t x, int64_t y, int64_t z) const;
    void setvoxelts(const ColumnVector& ts, int64_t x, int64_t y, int64_t z);
    void copySubVolume(const int64_t t, volume<T>& newVolume);
    ShadowVolume<T> operator[](const int64_t t);
    const ShadowVolume<T> operator[](const int64_t t) const;
    const volume<T> constSubVolume(const int64_t t) const;
    void replaceSubVolume(const int64_t t, const volume<T>& newVolume);
    void copyROI(const volume<T>& source,const int64_t minx, const int64_t miny, const int64_t minz, const int64_t mint, const int64_t maxx, const int64_t maxy, const int64_t maxz, const int64_t maxt);
    volume<T> ROI(const int64_t minx, const int64_t miny, const int64_t minz, const int64_t mint, const int64_t maxx, const int64_t maxy, const int64_t maxz, const int64_t maxt) const;

    // SECONDARY FUNCTIONS
    void setextrapolationmethod(extrapolation extrapmethod) const { p_extrapmethod = extrapmethod; }
    extrapolation getextrapolationmethod() const { return(p_extrapmethod); }

    void setpadvalue(T padval) const { padvalue = padval; }
    T getpadvalue() const { return padvalue; }
    T backgroundval() const;

    void defineuserextrapolation(T (*extrap)(
             const volume<T>& , int64_t, int64_t, int64_t)) const;
    void setinterpolationmethod(interpolation interpmethod) const;
    interpolation getinterpolationmethod() const { return p_interpmethod; }
    void setsplineorder(int order) const;
    void negateSpline() const {splineorder=-abs(splineorder);}
    int getsplineorder() const { return(splineorder); }
    void forcesplinecoefcalculation() const;
    void setextrapolationvalidity(bool xv, bool yv, bool zv) const { ep_valid[0]=xv; ep_valid[1]=yv; ep_valid[2]=zv; }
    std::vector<bool> getextrapolationvalidity() const { return(ep_valid); }
    void defineuserinterpolation(float (*interp)(
             const volume<T>& , float, float, float)) const;
    void definekernelinterpolation(const ColumnVector& kx,
				   const ColumnVector& ky,
				   const ColumnVector& kz,
				   int wx, int wy, int wz) const;  // full-width
    void definekernelinterpolation(const volume<T>& vol) const;
    void definesincinterpolation(const string& sincwindowtype,
				 int w, int nstore=1201) const;  // full-width
    void definesincinterpolation(const string& sincwindowtype,
				 int wx, int wy, int wz, int nstore=1201) const;
                                  // full-width

    inline void getneighbours(int64_t x, int64_t y, int64_t z,
			      T &v000, T &v001, T &v010,
			      T &v011, T &v100, T &v101,
			      T &v110, T &v111) const;
    inline void getneighbours(int64_t x, int64_t y, int64_t z,
			      T &v000, T &v010,
			      T &v100, T &v110) const;



    // ARITHMETIC FUNCTIONS
    T operator=(T val);
    const volume<T>& operator+=(T val);
    const volume<T>& operator-=(T val);
    const volume<T>& operator*=(T val);
    const volume<T>& operator/=(T val);
    const volume<T>& operator+=(const volume<T>& source);
    const volume<T>& operator-=(const volume<T>& source);
    const volume<T>& operator*=(const volume<T>& source);
    const volume<T>& operator/=(const volume<T>& source);

    volume<T> operator+(T num) const;
    volume<T> operator-(T num) const;
    volume<T> operator*(T num) const;
    volume<T> operator/(T num) const;
    volume<T> operator+(const volume<T>& vol2) const;
    volume<T> operator-(const volume<T>& vol2) const;
    volume<T> operator*(const volume<T>& vol2) const;
    volume<T> operator/(const volume<T>& vol2) const;

    template <class S>
    friend volume<S> operator+(S num, const volume<S>& vol);
    template <class S>
    friend volume<S> operator-(S num, const volume<S>& vol);
    template <class S>
    friend volume<S> operator*(S num, const volume<S>& vol);
    template <class S>
    friend volume<S> operator/(S num, const volume<S>& vol);
    template <class S>
    friend volume<S> operator-(const volume<S>& vol);

    // Comparisons. These are used for "spatial" purposes
    // so that if data is identical and all the "spatial
    // fields of the header are identical then the volumes
    // are considered identical.
    template <class S>
    friend bool operator==(const volume<S>& v1, const volume<S>& v2);
    template <class S>
    friend bool operator!=(const volume<S>& v1, const volume<S>& v2); // { return(!(v1==v2)); }

    // GENERAL MANIPULATION

    void binarise(T lowerth, T upperth, threshtype tt=inclusive, const bool invert=false);
    void binarise(T thresh) { this->binarise(thresh,this->max(),inclusive); }
    void threshold(T lowerth, T upperth, threshtype tt=inclusive);
    void threshold(T thresh) { this->threshold(thresh,this->max(),inclusive); }
    // valid entries for dims are +/- 1, 2, 3 (and for newx, etc they are x, -x, y, -y, z, -z)
    void swapdimensions(int dim1, int dim2, int dim3, bool keepLRorder=false);
    void swapdimensions(const string& newx, const string& newy, const string& newz, const bool keepLRorder=false);
    Matrix swapmat(int dim1, int dim2, int dim3) const;
    Matrix swapmat(const string& newx, const string& newy, const string& newz) const;


    // CONVERSION FUNCTIONS
    template <class S, class D> friend
      void copyconvert(const volume<S>& source, volume<D>& dest, const bool copyData);
  };

template<class T>
class ShadowVolume : public volume<T> {
 private:
  void setinterpolationmethod(interpolation interp) const { imthrow("Called private shadow method",101); }
  void setextrapolationmethod(extrapolation extrapmethod) const { imthrow("Called private shadow method",101); }

  virtual int isShadow() { return assigned ? assigned : assigned++; }
  int assigned;
 ShadowVolume() : assigned(0) {};
 public:
  const T operator=(const T val);
  const ShadowVolume<T>& operator=(const volume<T>& source);
  const ShadowVolume<T>& operator=(const ShadowVolume<T>& source);
 ShadowVolume(const volume<T>& source) : assigned(0) { this->reinitialize(source.xsize(),source.ysize(),source.zsize(),source.tsize(),const_cast<T*>(source.fbegin()),false); this->copyproperties(source); };
};

template<class T>
class baseIterator {
protected:
  T* ptr;
  T* const ptrBegin;
  T* const ptrEnd;
  const volume<T>* const source;
  baseIterator() : ptr(NULL), ptrBegin(NULL),  ptrEnd(NULL), source(NULL) { };
public:
  inline bool isValid() { return ptr != ptrEnd; }
 baseIterator(const volume<T>& input) : ptr((T*)input.fbegin()), ptrBegin(ptr), ptrEnd((T*)input.fend()), source(&input) { }
  baseIterator& operator++() { ++ptr; return *this; }
  baseIterator operator++(int) {baseIterator temp(*this); operator++(); return temp;}
  inline T& operator*() { return *ptr;}
};

template<class T>
  class cyclicIterator : public baseIterator<T> {
 private:
  cyclicIterator();
 public:
 cyclicIterator(const volume<T>& input) : baseIterator<T>(input) {};
  cyclicIterator& operator++() { if ( ++baseIterator<T>::ptr ==  baseIterator<T>::ptrEnd ) baseIterator<T>::ptr=baseIterator<T>::ptrBegin;  return *this; }
  cyclicIterator operator++(int) {cyclicIterator temp(*this); operator++(); return temp;}
 };

template<class T,class U>
class maskedIterator : public baseIterator<T> {
private:
  const U* mPtr;
  const U* const mPtrEnd;
  const volume<U>* const mask;
  bool usingMask;
  maskedIterator() : mPtr(NULL) {};
public:
  size_t offset() { return baseIterator<T>::ptr - baseIterator<T>::ptrBegin; }
  maskedIterator(const volume<T>& input,const volume<U>& inputMask,const string& exceptionSource="");
  maskedIterator& operator++();
  maskedIterator operator++(int) {maskedIterator temp(*this); operator++(); return temp;}
 };

template<class T, class U>
  maskedIterator<T,U>::maskedIterator(const volume<T>& input,const volume<U>& inputMask, const string& exceptionSource) : baseIterator<T>(input), mPtrEnd(inputMask.fend()), mask(&inputMask)
{
  mPtr=mask->fbegin(); //Should be safe even for NULL mask
  usingMask=( mask->totalElements() != 0 );

  if ( usingMask && !samesize(input,inputMask,SUBSET) )
    throw runtime_error(exceptionSource+"maskedIterator: mask and volume must be the same size");
  if ( usingMask && *mPtr <= mask->maskThreshold() )  //forward data pointer to first valid mask voxel
    operator++();
}

template<class T, class U>
maskedIterator<T,U>& maskedIterator<T,U>::operator++() {
  do {
  ++baseIterator<T>::ptr;
  if ( usingMask && ++mPtr == mPtrEnd ) { //always increment mPtr if used, reset mask for dimensional subset wrapping
     mPtr=mask->fbegin();
     if ( baseIterator<T>::ptr == baseIterator<T>::ptrEnd ) //Sanity check to prevent infinite loops for all-zero masks!
       return *this;
  }
  } while ( usingMask && *mPtr <= mask->maskThreshold() );
return *this;
}

template<class T>
  class safeLazyIterator : public baseIterator<T> { //Safe: cannot dereference if invalid, Lazy: invalidate source lazy flag on dereference
public:
  bool isValid() { return  baseIterator<T>::ptr != NULL &&  baseIterator<T>::ptr != baseIterator<T>::source->fend(); }
  safeLazyIterator(const volume<T>& input) : baseIterator<T>(input) {}
  safeLazyIterator& operator++() { ++baseIterator<T>::ptr ; return *this; }
  safeLazyIterator operator++(int) {safeLazyIterator temp(*this); operator++(); return temp;}
  T& operator*() { if(isValid()) return *baseIterator<T>::ptr; else throw runtime_error("Attempted to dereference an invalid safeLazyIterator");} //Need to invalidate source lazy flag
 };

template<class T>
  class lazyIterator : public baseIterator<T> {
public:
  lazyIterator(const volume<T>& input) : baseIterator<T>(input) {}
  lazyIterator& operator++() { ++baseIterator<T>::ptr ; return *this; }
  lazyIterator operator++(int) {lazyIterator temp(*this); operator++(); return temp;}
  T& operator*() { return *baseIterator<T>::ptr; } //Need to invalidate source
};



  // HELPER FUNCTIONS
  template <class S, class D>
  void convertbuffer(const S* source, D* dest, size_t len);

  template <class S, class D>
  void convertbuffer(const S* source, D* dest, size_t len, float slope, float intercept);

  template <class S1, class S2>
  bool samesize(const volume<S1>& vol1, const volume<S2>& vol2);
  template <class S1, class S2>
  bool samesize(const volume<S1>& vol1, const volume<S2>& vol2, bool checkdim);
  template <class S1, class S2>
  bool samesize(const volume<S1>& vol1, const volume<S2>& vol2, int ndims, bool checkdim=false);

  template <class S1, class S2>
  bool samedim(const volume<S1>& vol1, const volume<S2>& vol2, int ndims);

//////////////////////////////////////////////////////////////////////
///////////////////////// INLINE DEFINITIONS /////////////////////////
//////////////////////////////////////////////////////////////////////

  template <class T>
  inline void volume<T>::getneighbours(int64_t x, int64_t y, int64_t z,
					    T &v000, T &v001, T &v010,
					    T &v011, T &v100, T &v101,
					    T &v110, T &v111) const {
    T *ptr = basicptr(x,y,z);
    v000 = *ptr;
    ptr++;
    v100 = *ptr;
    ptr+= ColumnsX;
    v110 = *ptr;
    ptr--;
    v010 = *ptr;
    ptr += RowsY * ColumnsX;
    v011 = *ptr;
    ptr++;
    v111 = *ptr;
    ptr-= ColumnsX;
    v101 = *ptr;
    ptr--;
    v001 = *ptr;
  }


  template <class T>
  inline void volume<T>::getneighbours(int64_t x, int64_t y, int64_t z,
					    T &v000, T &v010,
					    T &v100, T &v110) const {
    T *ptr = basicptr(x,y,z);
    v000 = *ptr;
    ptr++;
    v100 = *ptr;
    ptr+= ColumnsX;
    v110 = *ptr;
    ptr--;
    v010 = *ptr;
  }





////////////////////////////////////////////////////////////////////////
/////////////////////////// HELPER FUNCTIONS ///////////////////////////
////////////////////////////////////////////////////////////////////////

  template <class T>
  int64_t no_mask_voxels(const volume<T>&  mask)
  {
     int64_t n(0);
     for(maskedIterator<T,T> it(mask,mask);it.isValid();++it)
       n++;
     return n;
  }


  template <class S, class D>
  void convertbuffer(const S* source, D* dest, size_t len)
  {
    D* dptr=dest;
    for (const S* sptr=source; sptr<(source+len); sptr++) {
      *dptr = (D) *sptr;
      dptr++;
    }
  }

  template <class S, class D>
  void convertbuffer(const S* source, D* dest, size_t len, float slope, float intercept)
  {
    D* dptr=dest;
    for (const S* sptr=source; sptr<(source+len); sptr++) {
      *dptr = (D) ((*sptr) * slope + intercept);
      dptr++;
    }
  }



  template <class S1, class S2>
  bool samesize(const volume<S1>& vol1, const volume<S2>& vol2, int ndims, bool checkdim)
  {
    int n=ndims;

    if (ndims==SUBSET) { n = Min(vol1.dimensionality(),vol2.dimensionality()); }
    bool same(vol1.xsize()==vol2.xsize());
    if (n>=2) same = same && ( vol1.ysize()==vol2.ysize());
    if (n>=3) same = same && ( vol1.zsize()==vol2.zsize());
    if (n>=4) same = same && ( vol1.tsize()==vol2.tsize());
    if (n>=5) same = same && ( vol1.size5()==vol2.size5());
    if (n>=6) same = same && ( vol1.size6()==vol2.size6());
    if (n>=7) same = same && ( vol1.size7()==vol2.size7());
    if (checkdim)
      same = same && samedim(vol1,vol2,n);
    return(same);
  }

  template <class S1, class S2>
  bool samesize(const volume<S1>& vol1, const volume<S2>& vol2, bool checkdim)
  {
    return samesize(vol1,vol2,7,checkdim);
  }

  template <class S1, class S2>
  bool samesize(const volume<S1>& vol1, const volume<S2>& vol2)
  {
    return samesize(vol1,vol2,7,false);
  }


  template <class S1, class S2>
  bool samedim(const volume<S1>& vol1, const volume<S2>& vol2, int ndims)
  {
    bool same=true;
    same = (std::fabs(vol1.xdim()-vol2.xdim())<1e-6);
    if (ndims>=2) same = same && (std::fabs(vol1.ydim()-vol2.ydim())<1e-6);
    if (ndims>=3) same = same && (std::fabs(vol1.zdim()-vol2.zdim())<1e-6);
    if (ndims>=4) same = same && (std::fabs(vol1.tdim()-vol2.tdim())<1e-6);
    if (ndims>=5) same = same && (std::fabs(vol1.pixdim5()-vol2.pixdim5())<1e-6);
    if (ndims>=6) same = same && (std::fabs(vol1.pixdim6()-vol2.pixdim6())<1e-6);
    if (ndims>=7) same = same && (std::fabs(vol1.pixdim7()-vol2.pixdim7())<1e-6);
    return same;
  }



////////////////////////////////////////////////////////////////////////
/////////////////////////// FRIEND FUNCTIONS ///////////////////////////
////////////////////////////////////////////////////////////////////////

  template <class S, class D>
  void copybasicproperties(const volume<S>& source, volume<D>& dest)
  {
    // set up properties
    dest.Xdim = source.Xdim;
    dest.Ydim = source.Ydim;
    dest.Zdim = source.Zdim;
    dest.p_TR = source.p_TR;
    dest.pxdim5 = source.pxdim5;
    dest.pxdim6 = source.pxdim6;
    dest.pxdim7 = source.pxdim7;

    dest.originalSizes = source.originalSizes;

    dest.StandardSpaceCoordMat = source.StandardSpaceCoordMat;
    dest.RigidBodyCoordMat = source.RigidBodyCoordMat;
    dest.StandardSpaceTypeCode = source.StandardSpaceTypeCode;
    dest.RigidBodyTypeCode = source.RigidBodyTypeCode;

    dest.RadiologicalFile = source.RadiologicalFile;

    dest.IntentCode = source.IntentCode;
    dest.IntentParam1 = source.IntentParam1;
    dest.IntentParam2 = source.IntentParam2;
    dest.IntentParam3 = source.IntentParam3;

    dest.SliceOrderingCode = source.SliceOrderingCode;

    dest.interpkernel = source.interpkernel;
    dest.p_interpmethod = source.p_interpmethod;
    dest.p_extrapmethod = source.p_extrapmethod;

    dest.padvalue = (D) source.padvalue;
    dest.splineorder = source.splineorder;
    dest.ep_valid = source.ep_valid;

    dest.displayMaximum=source.displayMaximum;
    dest.displayMinimum=source.displayMinimum;
    dest.setAuxFile(source.getAuxFile());
  }


  template <class S, class D>
  void copyconvert(const volume<S>& source, volume<D>& dest, const bool copyData=true)
  {
    // set up basic size and data storage
    if ( dest.totalElements() != source.totalElements() )
      dest.reinitialize(source.xsize(),source.ysize(),source.zsize(),source.tsize(),source.size5(),source.size6(),source.size7());
    // set up properties (except lazy ones)
    copybasicproperties(source,dest);
    // now copy across the data
    if ( copyData )
      convertbuffer(source.Data, dest.Data, source.totalElements() );
  }

  template <class S>
  volume<S> operator+(S num, const volume<S>& vol)
    { return (vol + num); }

  template <class S>
  volume<S> operator-(S num, const volume<S>& vol)
  {
    volume<S> tmp = vol;
    tmp=num;
    tmp-=vol;
    return tmp;
  }

  template <class S>
  volume<S> operator*(S num, const volume<S>& vol)
    { return (vol * num); }

  template <class S>
  volume<S> operator/(S num, const volume<S>& vol)
  {
    volume<S> tmp = vol;
    tmp=num;
    tmp/=vol;
    return tmp;
  }

  template <class S>
  volume<S> operator-(const volume<S>& vol)
  {
    return(vol * (static_cast<S>(-1)));
  }

  template <class S>
  bool operator==(const volume<S>& v1,
		  const volume<S>& v2)
  {
    // Check relevant parts of header
    if (!samesize(v1,v2,true)) return(false);
    if (v1.sform_code() != v2.sform_code()) return(false);
    if (v1.sform_mat() != v2.sform_mat()) return(false);
    if (v1.qform_code() != v2.qform_code()) return(false);
    if (v1.qform_mat() != v2.qform_mat()) return(false);
    // Check data
    for (typename volume<S>::fast_const_iterator it1=v1.fbegin(), it_end=v1.fend(), it2=v2.fbegin(); it1 != it_end; ++it1, ++it2) {
      if ((*it1) != (*it2)) return(false);
    }

    return(true);
  }
  template <class S>
  bool operator!=(const volume<S>& v1, const volume<S>& v2) {return(!(v1==v2));}

  template <class T>
  const T ShadowVolume<T>::operator=(const T val)
  {
    fill(this->nsfbegin(),this->nsfend(),val);  // use the STL
    return val;
  }


  template < class S, class V >
  vector<double> calculateSums(const S& inputVolume, const V& mask );
  template < class T, class V>
  vector<T> calculateExtrema(const volume<T>& inputVolume, vector<int64_t>& coordinates, const volume<V>& mask );

}  // end namespace

#endif
