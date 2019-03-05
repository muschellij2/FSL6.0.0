/*  newimage.cc

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

#include <complex>
#include <cassert>
#include <sstream>
#include <iostream>
#include <boost/utility/enable_if.hpp>
#include <boost/static_assert.hpp>
#include "newmatio.h"
#include "newimage.h"


using namespace NEWMAT;
using namespace MISCMATHS;

namespace NEWIMAGE {
  void imthrow(const string& msg, int nierrnum, const bool quiet) {
    if(!quiet)
      cerr << "Image Exception : #" << nierrnum << " :: " << msg << endl;
    if (USEASSERT) { bool no_error=false; assert(no_error); }
    else { throw Exception(msg.data()); }
  }

  template <class T>
  T calc_backgroundval(const volume<T>& vol);

  template <class T>
  ColumnVector calc_cog(const volume<T>& vol);

  template <class T>
  std::vector<T> calc_robustlimits(const volume<T>& vol);

  template <class T>
  Matrix calc_principleaxes(const volume<T>& vol);

  template <class T>
  std::vector<T> calc_percentiles(const volume<T>& vol);

  template < class T, class V>
  std::vector<T> calc_percentiles(const volume<T>& vol, const volume<T>& mask,
				  const std::vector<float>& percentilepvals);

  template <class T, class V>
  int64_t imageHistogram(const volume<T>& vol, const int bins,
			 const T min, const T max, ColumnVector& hist, const volume<V>& mask);

  template <class T, class V>
  int64_t imageHistogram(const volume<T>& vol, const int bins,
			 const T min, const T max, ColumnVector& hist);

  // Declaration of helper functions.

  SPLINTERPOLATOR::ExtrapolationType translate_extrapolation_type(extrapolation ep);


  template<class T>
  vector<int64_t> volume<T>::ptrToCoord(size_t offset) const
  {
  vector<int64_t> output;
    output.push_back(offset%xsize());
    offset/=xsize();
    output.push_back(offset%ysize());
    offset/=ysize();
    output.push_back(offset%zsize());
    offset/=zsize();
    output.push_back(offset%tsize());
    offset/=tsize();
    output.push_back(offset%size5());
    offset/=size5();
    output.push_back(offset%size6());
    offset/=size6();
    output.push_back(offset);
    return(output);
  }

  template <class T>
  int64_t volume<T>::size(int n) const
  {
    if (n==1) return size1();
    if (n==2) return size2();
    if (n==3) return size3();
    if (n==4) return size4();
    if (n==5) return size5();
    if (n==6) return size6();
    if (n==7) return size7();
    return 0;
  }


  // CONSTRUCTORS (not including copy constructor - see under copying)
 template <class T>
  int volume<T>::initialize(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize, int64_t d5, int64_t d6, int64_t d7, T *d, bool d_owner)
  {
    if ( this->isShadow() )
      imthrow("Attempted to reinitialise shadow volume",99);
    this->destroy();
    SlicesZ = zsize;
    RowsY = ysize;
    ColumnsX = xsize;
    dim4=tsize;
    dim5=d5;
    dim6=d6;
    dim7=d7;
    nElements=xsize*ysize*zsize*tsize*d5*d6*d7;
    no_voxels = xsize*ysize*zsize;
    maskDelimiter=0.5;
    // decide whether to allocate new memory or not, depending on validity of pointer
    if (nElements > 0) {
      if (d != NULL) {
	Data = d;
	DataEnd = d+nElements;
	data_owner = d_owner;
      } else {
	try {
	  Data = new T[nElements];
	  DataEnd = Data+nElements;
	} catch(...) { Data=NULL; DataEnd=NULL;}
	if (Data==NULL) { imthrow("Out of memory",99); }
	data_owner = true;
      }
    } else {
      Data = NULL;
      DataEnd = NULL;
      data_owner = false;
    }
    setdefaultproperties();
    return 0;
  }

  template <class T>
  int volume<T>::initialize(int64_t xsize, int64_t ysize, int64_t zsize, T *d, bool d_owner)
  {
      return initialize(xsize,ysize,zsize,1,1,1,1,d,d_owner);
  }

  template <class T>
  int volume<T>::initialize(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize, T *d, bool d_owner)
  {
      return initialize(xsize,ysize,zsize,tsize,1,1,1,d,d_owner);
  }

  template <class T>
  void volume<T>::setdefaultproperties()
    {
      Xdim = 1.0;
      Ydim = 1.0;
      Zdim = 1.0;
      p_TR = 1.0;
      pxdim5=1.0;
      pxdim6=1.0;
      pxdim7=1.0;

      originalSizes.resize(8,-1); //if originalSizes[0] is _not_ -1 then there has been alteration

      StandardSpaceCoordMat = IdentityMatrix(4);
      RigidBodyCoordMat = IdentityMatrix(4);
      StandardSpaceTypeCode = NIFTI_XFORM_UNKNOWN;
      RigidBodyTypeCode = NIFTI_XFORM_UNKNOWN;
      RadiologicalFile = true;

      IntentCode = NIFTI_INTENT_NONE;
      IntentParam1 = 0.0;
      IntentParam2 = 0.0;
      IntentParam3 = 0.0;

      SliceOrderingCode = NIFTI_SLICE_UNKNOWN;

      p_interpmethod = trilinear;
      p_extrapmethod = zeropad;
      splineorder = -3; //Negative splineorder indicates splines no longer valid

      padvalue = (T) 0;
      extrapval = padvalue;
      p_userinterp = 0;
      p_userextrap = 0;
      ep_valid.resize(3);
      ep_valid[0] = false; ep_valid[1] = false; ep_valid[2] = false;

      displayMaximum=0;
      displayMinimum=0;
      strncpy(auxFile,string("").c_str(),24);

    }

  template <class T>
  volume<T>::volume() : Data(0), DataEnd(0), data_owner(false)
    {
      this->initialize(0,0,0,0,0,0,0,0,false);
    }

  template <class T>
  volume<T>::volume(int64_t xsize, int64_t ysize, int64_t zsize) : Data(0), DataEnd(0), data_owner(false)
    {
      this->initialize(xsize,ysize,zsize,1,1,1,1,0,true);
    }

  template <class T>
  volume<T>::volume(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize) : Data(0), DataEnd(0), data_owner(false)
    {
      this->initialize(xsize,ysize,zsize,tsize,1,1,1,0,true);
    }


  template <class T>
  volume<T>::volume(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize, int64_t d5, int64_t d6, int64_t d7) : Data(0), DataEnd(0), data_owner(false)
    {
      this->initialize(xsize,ysize,zsize,tsize,d5,d6,d7,0,true);
    }

  template <class T>
  int volume<T>::reinitialize(int64_t xsize, int64_t ysize, int64_t zsize)
    {
      return this->initialize(xsize,ysize,zsize,0,true);
    }

  template <class T>
  int volume<T>::reinitialize(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize, T *d, bool d_owner)
    {
      return this->initialize(xsize,ysize,zsize,tsize,d,d_owner);
    }

  template <class T>
  int volume<T>::reinitialize(int64_t xsize, int64_t ysize, int64_t zsize, int64_t tsize, int64_t d5, int64_t d6, int64_t d7)
    {
      return this->initialize(xsize,ysize,zsize,tsize,d5,d6,d7,0,true);
    }



  template <class T>
  void volume<T>::destroy()
  {
    if ( data_owner && Data != NULL ) delete [] Data;
    Data = NULL;
    DataEnd = NULL;
    data_owner=false;
    // make the volume of zero size now (to prevent access to the null data pointer)
    nElements=0;
    ColumnsX=0;
    RowsY=0;
    SlicesZ=0;
    dim4=0;
    dim5=0;
    dim6=0;
    dim7=0;
  }

  template <class T>
  volume<T>::~volume()
    {  this->destroy(); }

  // COPYING AND CONVERSION FUNCTIONS

  template <class T>
  void volume<T>::reinitialize(const volume<T>& source, const bool copyingData)
    {
      this->initialize(source.xsize(),source.ysize(),source.zsize(),source.tsize(),source.size5(),source.size6(),source.size7(),0,false);
      if ( copyingData )
	      this->copydata(source);
      this->copyproperties(source);
    }

  template <class T>
  volume<T>::volume(const volume<T>& source, const bool copyingData) : Data(0), DataEnd(0), data_owner(false)
    {
      this->reinitialize(source, copyingData);
    }

  template <class T>
  ShadowVolume<T> volume<T>::operator[](const int64_t t) {
    if ( !in_bounds(t) )
       imthrow("Invalid t index in [] operator",61);
    ShadowVolume<T> newShadowVolume(constSubVolume(t));
    newShadowVolume.copyproperties(*this);
    return newShadowVolume;
  }

  template <class T>
  const ShadowVolume<T> volume<T>::operator[](const int64_t t) const {
    if ( !in_bounds(t) )
       imthrow("Invalid t index in [] operator",61);
    ShadowVolume<T> newShadowVolume(constSubVolume(t));
    newShadowVolume.copyproperties(*this);
    return newShadowVolume;
  }


  template <class T>
  const volume<T>& volume<T>::operator=(const volume<T>& source)
  {
    this->reinitialize(source);
    return *this;
  }

  template <class T>
  const volume<T>& volume<T>::operator=(const ShadowVolume<T>& source)
  {
    this->reinitialize(source);
    return *this;
  }

  template <class T>
  const ShadowVolume<T>& ShadowVolume<T>::operator=(const ShadowVolume<T>& source)
  {
    this->copydata(source,false);
    return *this;
  }

  template <class T>
  const ShadowVolume<T>& ShadowVolume<T>::operator=(const volume<T>& source)
  {
    this->copydata(source,false);
    return *this;
  }

  template <class T>
  int volume<T>::copydata(const volume<T>& source, const bool isOwner) {
    if (nElements != source.nElements) {
      imthrow("Attempted to copydata with non-matching sizes",2);
    }
    copy(source.Data, source.Data + nElements, nsfbegin());  // use the STL
    data_owner = isOwner;
    return 0;
  }

  template <class T>
  void volume<T>::copyproperties(const volume<T>& source)
  { //the user extraps require matching template in source and dest
    p_userextrap = source.p_userextrap;
    p_userinterp = source.p_userinterp;
    copybasicproperties(source,*this);
  }


  // VOLUME ACCESS FUNCTIONS (INSERT AND DELETE)

  template <class T>
  int volume<T>::dimensionality() const
  {
    int dim=0;
    if (xsize()>=1) dim=1;
    if (ysize()>1) dim=2;
    if (zsize()>1) dim=3;
    if (tsize()>1) dim=4;
    if (size5()>1) dim=5;
    if (size6()>1) dim=6;
    if (size7()>1) dim=7;
    return dim;
  }

  template <class T>
  void volume<T>::copySubVolume(const int64_t t, volume<T>& newSubVolume) {
    newSubVolume.reinitialize(xsize(),ysize(),zsize());
    newSubVolume.copyproperties(*this);
    copy(this->Data + t*nvoxels(), this->Data + (t+1)*nvoxels(), newSubVolume.Data);  // use the STL
  }

  template <class T>
  const volume<T> volume<T>::constSubVolume(const int64_t t) const {
    volume<T> newSubVolume;
    newSubVolume.initialize(xsize(),ysize(),zsize(),Data+t*nvoxels(),false);
    newSubVolume.copyproperties(*this);
    return newSubVolume;
  }

  template <class T>
  void volume<T>::replaceSubVolume(const int64_t t, const volume<T>& newVolume) {
    if (newVolume.dimensionality()>3) { imthrow("Attempted to replaceSubVolume with non-3D input",2); }
    if (!samesize(*this,newVolume,SUBSET)) { imthrow("Attempted to replaceSubVolume with non-matching sizes",2); }
    copy(newVolume.Data, newVolume.Data + nvoxels(), this->Data + t*nvoxels());  // use the STL
  }

  template <class T>
  void volume<T>::addvolume(const volume<T>& source)
  {
    if ( Data == NULL ) { //Blank volume
      *this=source;
      return;
    }

    if (!samesize(*this,source,3)) { imthrow("Attempted to addvolume when 3D sizes were not the same",2); }
    int concatdim=Max(this->dimensionality(),source.dimensionality());  // default to concat in highest non-trivial dim
    // check dimensionalities: cannot have mismatching dims for any dim less than the highest non-trivial dim
    bool errflag=false;
    if (this->size4()!=source.size4()) { if (concatdim>4) errflag=true; }
    if (this->size5()!=source.size5()) { if (concatdim>5) errflag=true; }
    if (this->size6()!=source.size6()) { if (concatdim>6) errflag=true; }
    if (errflag) { imthrow("Attempted to addvolume with inconsistent dimensions",60); }
    if (concatdim<4) concatdim=4;  // if it is only two 3D vols then still concat in the 4th dim

    int64_t out4=this->size4(), out5=this->size5(), out6=this->size6(), out7=this->size7();
    if (concatdim==4) { out4+=source.size4(); }
    if (concatdim==5) { out5+=source.size5(); }
    if (concatdim==6) { out6+=source.size6(); }
    if (concatdim==7) { out7+=source.size7(); }
    volume output(this->xsize(),this->ysize(),this->zsize(),out4,out5,out6,out7);
    nonsafe_fast_iterator dit=output.nsfbegin();
    for (fast_const_iterator sit=this->fbegin(), sitend=this->fend(); sit!=sitend; ++sit, ++dit) { *dit = *sit; }
    for (fast_const_iterator sit=source.fbegin(), sitend=source.fend(); sit!=sitend; ++sit, ++dit) { *dit = *sit; }
    output.copyproperties(*this);
    this->reinitialize(output);
  }

  template <class T>
  void volume<T>::deletevolume(int64_t t)
  {
    if ((t>=this->tsize()) || (t<0)) { imthrow("Invalid t index in deletevolume",61); }
    // check that there is only one non-trivial dimension in dims4,5,6,7  (the dimension we shall delete over)
    int64_t proddims=size4()*size5()*size6()*size7();
    int deldim=this->dimensionality();
    if (deldim<4) { this->destroy(); return; }  // if you delete the last volume then everything's gone!
    if (proddims!=this->size(deldim)) { imthrow("Attempted to deletevolume with inconsistent dimensions",60); }
    // determine new output size
    int64_t out4=this->size4(), out5=this->size5(), out6=this->size6(), out7=this->size7();
    if (deldim==4) { out4--; }
    if (deldim==5) { out5--; }
    if (deldim==6) { out6--; }
    if (deldim==7) { out7--; }
    volume output(this->xsize(),this->ysize(),this->zsize(),out4,out5,out6,out7);
    nonsafe_fast_iterator dit=output.nsfbegin(), sit=this->nsfbegin(), sitend1, sitend2=this->nsfend();
    if (t>0) { // initial t-1 vols to copy
      for (sitend1=this->nsfbegin() + t*nvoxels(); sit!=sitend1; ++sit, ++dit) { *dit = *sit; }
    }
    if (t<this->tsize()-1) {  // remaining vols to copy
      for (sit+=nvoxels() ; sit!=sitend2; ++sit, ++dit) { *dit = *sit; }
    }
    output.copyproperties(*this);
    this->reinitialize(output);
  }


  template <class T>
  void volume<T>::clear()
  {
    this->destroy();
  }


  // Volume->ColumnVector

  template<class T>
  ReturnMatrix volume<T>::vec(const volume<T>& mask) const
  {
    if (this->dimensionality()>3) { imthrow("volume<T>::vec: input vol must be 3D",7); }
    if (!samesize(mask,*this,3)) {imthrow("volume<T>::vec: Mask and volume of different size",3);}
    ColumnVector ovec(xsize()*ysize()*zsize());
    for (int vindx=0,k=0; k<zsize(); k++) {
      for (int j=0; j<ysize(); j++) {
	for (int i=0; i<xsize(); i++) {
          ovec.element(vindx) = (mask(i,j,k)>0) ? (*this)(i,j,k) : 0.0;
          vindx++;
	}
      }
    }
    ovec.Release();
    return ovec;
  }

  // Code multiplication to avoid allocating mask volume

  template <class T>
  ReturnMatrix volume<T>::vec() const
  {
    if (this->dimensionality()>3) { imthrow("volume<T>::vec: input vol must be 3D",7); }
    ColumnVector ovec(xsize()*ysize()*zsize());
    for (int vindx=0, k=0; k<zsize(); k++) {
      for (int j=0; j<ysize(); j++) {
	for (int i=0; i<xsize(); i++) {
          ovec.element(vindx) = (*this)(i,j,k);
          vindx++;
	}
      }
    }
    ovec.Release();
    return ovec;
  }

  // Insert ColumnVector into volume

  template <class T>
  void volume<T>::insert_vec(const ColumnVector&  pvec,
                             const volume<T>&     mask)
  {
    if (pvec.Nrows() != no_mask_voxels(mask)) {
      imthrow("volume<T>::insert_vec: Size mismatch between ColumnVector and mask volume",3);
    }
    if (!samesize(mask,*this,3)) {
      this->reinitialize(mask.xsize(),mask.ysize(),mask.zsize());
      this->copyproperties(mask);
      this->operator=((T)0);
    }
    for (int vindx=0, k=0; k<zsize(); k++) {
      for (int j=0; j<ysize(); j++) {
	for (int i=0; i<xsize(); i++) {
	  (*this)(i,j,k) = (mask(i,j,k) > 0) ? ((T) pvec.element(vindx)) : ((T) 0);
          vindx++;
	}
      }
    }
  }

  template <class T>
  void volume<T>::insert_vec(const ColumnVector&  pvec)
  {
    if (pvec.Nrows() != xsize()*ysize()*zsize()) {
      imthrow("volume<T>::insert_vec: Size mismatch between ColumnVector and image volume",3);
    }
    for (int vindx=0, k=0; k<zsize(); k++) {
      for (int j=0; j<ysize(); j++) {
	for (int i=0; i<xsize(); i++) {
	  (*this)(i,j,k) = ((T) pvec.element(vindx));
          vindx++;
	}
      }
    }
  }

 template <class T>
  vector<int64_t> volume<T>::labelToCoord(const int64_t label) const
  {
    vector<int64_t> coordinates;
    coordinates.push_back(label%this->xsize());
    coordinates.push_back( (floor) ( ( label%( this->xsize()*this->ysize() ) ) / this->xsize() ));
    coordinates.push_back( (floor) ( label / ( this->xsize()*this->ysize() ) ) );
    return coordinates;
  }
  // EXTRAPOLATION AND INTERPOLATION

  template <class T>
  void volume<T>::setinterpolationmethod(interpolation interpmethod) const
    {
      p_interpmethod = interpmethod;
      // define a default sinc kernel if no kernel has previously been defined
      if ( (interpmethod == sinc) && (interpkernel.kernelvals()==0) ) {
	string sincwindowtype = "blackman";
	this->definesincinterpolation(sincwindowtype,7);
      }
    }

  template<class T>
  void volume<T>::setsplineorder(int order) const
  {
    if (order > 7) imthrow("setsplineorder: Only splines of order up to 7 allowed",10);
    splineorder = order;
  }

  template <class T>
  bool in_neigh_bounds(const volume<T>& vol, int64_t x, int64_t y, int64_t z)
    {  return ( (x>=0) && (y>=0) && (z>=0) &&
		(x<(vol.xsize()-1)) && (y<(vol.ysize()-1)) &&
		(z<(vol.zsize()-1)) ); }

  inline float q_tri_interpolation(float v000, float v001, float v010,
				   float v011, float v100, float v101,
				   float v110, float v111,
				   float dx, float dy, float dz)

    {
      float temp1, temp2, temp3, temp4, temp5, temp6;
      temp1 = (v100 - v000)*dx + v000;
      temp2 = (v101 - v001)*dx + v001;
      temp3 = (v110 - v010)*dx + v010;
      temp4 = (v111 - v011)*dx + v011;
      // second order terms
      temp5 = (temp3 - temp1)*dy + temp1;
      temp6 = (temp4 - temp2)*dy + temp2;
      // final third order term
      return (temp6 - temp5)*dz + temp5;
    }

  //////// Kernel Interpolation Call /////////

  template <class T>
  float volume<T>::kernelinterpolation(const float x, const float y,
				       const float z) const
  {
    const kernelstorage* storedkernel = interpkernel.kernelvals();
    // sanity check on kernel
    if (storedkernel==0) {
      cerr << "ERROR: Must set kernel parameters before using interpolation!"
	   << endl;
      return (float) extrapolate(0,0,0);
    }

    // kernel half-width  (i.e. range is +/- w)
    int wx=storedkernel->widthx();
    int wy=storedkernel->widthy();
    int wz=storedkernel->widthz();
    ColumnVector kernelx = storedkernel->kernelx();
    ColumnVector kernely = storedkernel->kernely();
    ColumnVector kernelz = storedkernel->kernelz();
    float *storex = storedkernel->storex;
    float *storey = storedkernel->storey;
    float *storez = storedkernel->storez;

    int64_t ix0, iy0, iz0;
    ix0 = (int64_t) floor(x);
    iy0 = (int64_t) floor(y);
    iz0 = (int64_t) floor(z);

    float convsum=0.0, interpval=0.0, kersum=0.0;

    for (int d=-wz; d<=wz; d++) {
      storez[d+wz] = kernelval((z-iz0+d),wz,kernelz);
    }
    for (int d=-wy; d<=wy; d++) {
      storey[d+wy] = kernelval((y-iy0+d),wy,kernely);
    }
    for (int d=-wx; d<=wx; d++) {
      storex[d+wx] = kernelval((x-ix0+d),wx,kernelx);
    }

    int64_t xj, yj, zj;
    for (int64_t z1=iz0-wz; z1<=iz0+wz; z1++) {
      zj=iz0-z1+wz;
      for (int64_t y1=iy0-wy; y1<=iy0+wy; y1++) {
	yj=iy0-y1+wy;
	for (int64_t x1=ix0-wx; x1<=ix0+wx; x1++) {
	  if (in_bounds(x1,y1,z1)) {
	    xj=ix0-x1+wx;
	    float kerfac = storex[xj] * storey[yj] * storez[zj];
	    convsum += this->operator()(x1,y1,z1) * kerfac;
	    kersum += kerfac;
	  }
	}
      }
    }

    if ( (fabs(kersum)>1e-9) ) {
      interpval = convsum / kersum;
    } else {
      interpval = (float) extrapolate(ix0,iy0,iz0);
    }
    return interpval;

  }

  // The following routines are used to obtain an interpolated intensity value and
  // either a selected partial derivative (dx, dy or dz) or all partial derivatives
  // at the same location. The routine returning all derivatives is useful for
  // non-linear registration of one subject to another (or to an atlas) and the
  // routines returning a single derivative are useful e.g. for distortion correction.
  // Puss J

  template <class T>
  float volume<T>::interp1partial(// Input
                                  float x, float y, float z,    // Co-ordinates to get value for
                                  int     dir,                  // Direction for partial, 0->x, 1->y, 2->z
                                  // Output
                                  float  *pderiv)               // Derivative returned here
  const
  {
    if (getinterpolationmethod() != trilinear && getinterpolationmethod() != spline) {
      imthrow("Derivatives only implemented for tri-linear and spline interpolation",10);
    }
    if (dir < 0 || dir > 2) {
      imthrow("Ivalid derivative direction",11);
    }
    if (getinterpolationmethod() == trilinear) {
      int64_t ix = ((int64_t) floor(x));
      int64_t iy = ((int64_t) floor(y));
      int64_t iz = ((int64_t) floor(z));
      float dx = x - ((float) ix);
      float dy = y - ((float) iy);
      float dz = z - ((float) iz);
      float v000, v001, v010, v011, v100, v101, v110, v111;
      if (!in_neigh_bounds(*this,ix,iy,iz)) {   // We'll have to do some extrapolation
	v000 = (float) this->operator()(ix,iy,iz);
	v001 = (float) this->operator()(ix,iy,iz+1);
	v010 = (float) this->operator()(ix,iy+1,iz);
	v011 = (float) this->operator()(ix,iy+1,iz+1);
	v100 = (float) this->operator()(ix+1,iy,iz);
	v101 = (float) this->operator()(ix+1,iy,iz+1);
	v110 = (float) this->operator()(ix+1,iy+1,iz);
	v111 = (float) this->operator()(ix+1,iy+1,iz+1);
      }
      else {
	T t000, t001, t010, t011, t100, t101, t110, t111;
	this->getneighbours(ix,iy,iz,t000,t001,t010,t011,t100,t101,t110,t111);
	v000 = ((float) t000); v001 = ((float) t001); v010 = ((float) t010);
	v011 = ((float) t011); v100 = ((float) t100); v101 = ((float) t101);
	v110 = ((float) t110); v111 = ((float) t111);
      }
      // The (seemingly silly) code multiplication below is to
      // ensure that in no case does calculating one of the partials
      // neccessitate any calculation over and above just calculating
      // the interpolated value.
      float tmp11, tmp12, tmp13, tmp14;
      float tmp21, tmp22;
      if (dir == 0) {            // df/dx
	float onemdz = 1.0-dz;
	tmp11 = onemdz*v000 + dz*v001;
	tmp12 = onemdz*v010 + dz*v011;
	tmp13 = onemdz*v100 + dz*v101;
	tmp14 = onemdz*v110 + dz*v111;
	tmp21 = (1.0-dy)*tmp11 + dy*tmp12;
	tmp22 = (1.0-dy)*tmp13 + dy*tmp14;
	*pderiv = tmp22 - tmp21;
	return((1.0-dx)*tmp21 + dx*tmp22);
      }
      else if (dir == 1) {       // df/dy
	float onemdz = 1.0-dz;
	tmp11 = onemdz*v000 + dz*v001;
	tmp12 = onemdz*v010 + dz*v011;
	tmp13 = onemdz*v100 + dz*v101;
	tmp14 = onemdz*v110 + dz*v111;
	tmp21 = (1.0-dx)*tmp11 + dx*tmp13;
	tmp22 = (1.0-dx)*tmp12 + dx*tmp14;
	*pderiv = tmp22 - tmp21;
	return((1.0-dy)*tmp21 + dy*tmp22);
      }
      else if (dir == 2) {       // df/dz
	float onemdy = 1.0-dy;
	tmp11 = onemdy*v000 + dy*v010;
	tmp12 = onemdy*v001 + dy*v011;
	tmp13 = onemdy*v100 + dy*v110;
	tmp14 = onemdy*v101 + dy*v111;
	tmp21 = (1.0-dx)*tmp11 + dx*tmp13;
	tmp22 = (1.0-dx)*tmp12 + dx*tmp14;
	*pderiv = tmp22 - tmp21;
	return((1.0-dz)*tmp21 + dz*tmp22);
      }
    }
    else if (getinterpolationmethod() == spline) {
      return(spline_interp1partial(x,y,z,dir,pderiv));
    }
    return(-1.0); // Should not be reached. Just to stop compiler from complaining.
  }

  template <class T>
  float volume<T>::interp3partial(// Input
                                  float x, float y, float z,              // Co-ordinates to get value for
                                   // Output
                                  float *dfdx, float *dfdy, float *dfdz)  // Partials
  const
  {
    if (getinterpolationmethod() != trilinear && getinterpolationmethod() != spline) {
      imthrow("interp3partial: Derivatives only implemented for tri-linear and spline interpolation",10);
    }
    if (getinterpolationmethod() == trilinear) {
      int64_t ix = ((int64_t) floor(x));
      int64_t iy = ((int64_t) floor(y));
      int64_t iz = ((int64_t) floor(z));
      float dx = x - ((float) ix);
      float dy = y - ((float) iy);
      float dz = z - ((float) iz);
      float v000, v001, v010, v011, v100, v101, v110, v111;
      if (!in_neigh_bounds(*this,ix,iy,iz)) {   // We'll have to do some extrapolation
	v000 = (float) this->operator()(ix,iy,iz);
	v001 = (float) this->operator()(ix,iy,iz+1);
	v010 = (float) this->operator()(ix,iy+1,iz);
	v011 = (float) this->operator()(ix,iy+1,iz+1);
	v100 = (float) this->operator()(ix+1,iy,iz);
	v101 = (float) this->operator()(ix+1,iy,iz+1);
	v110 = (float) this->operator()(ix+1,iy+1,iz);
	v111 = (float) this->operator()(ix+1,iy+1,iz+1);
      }
      else {
	T t000, t001, t010, t011, t100, t101, t110, t111;
	this->getneighbours(ix,iy,iz,t000,t001,t010,t011,t100,t101,t110,t111);
	v000 = ((float) t000); v001 = ((float) t001); v010 = ((float) t010);
	v011 = ((float) t011); v100 = ((float) t100); v101 = ((float) t101);
	v110 = ((float) t110); v111 = ((float) t111);
      }
      //
      // And do linear interpolation with calculation of all partials
      //
      float onemdz = 1.0-dz;
      float onemdy = 1.0-dy;
      float tmp11 = onemdz*v000 + dz*v001;
      float tmp12 = onemdz*v010 + dz*v011;
      float tmp13 = onemdz*v100 + dz*v101;
      float tmp14 = onemdz*v110 + dz*v111;
      *dfdx = onemdy*(tmp13-tmp11) + dy*(tmp14-tmp12);
      *dfdy = (1.0-dx)*(tmp12-tmp11) + dx*(tmp14-tmp13);
      tmp11 = onemdy*v000 + dy*v010;
      tmp12 = onemdy*v001 + dy*v011;
      tmp13 = onemdy*v100 + dy*v110;
      tmp14 = onemdy*v101 + dy*v111;
      float tmp21 = (1.0-dx)*tmp11 + dx*tmp13;
      float tmp22 = (1.0-dx)*tmp12 + dx*tmp14;
      *dfdz = tmp22 - tmp21;
      return(onemdz*tmp21 + dz*tmp22);
    }
    else if (getinterpolationmethod() == spline) {
      return(spline_interp3partial(x,y,z,dfdx,dfdy,dfdz));
    }
    return(0.0);  // To silence compiler.
  }

  template <class T>
  float volume<T>::spline_interp1partial(// Input
                                         float x, float y, float z,    // Co-ordinates to get value for
                                         int     dir,                  // Direction for partial, 0->x, 1->y, 2->z
                                         // Output
                                         float  *deriv)               // Derivative returned here
  const
  {
    if (!in_bounds(x,y,z)) {
      extrapolation ep = getextrapolationmethod();
      if (ep == boundsassert) { *deriv=0.0; assert(false); extrapval = padvalue; return(extrapval); }
      else if (ep == boundsexception) imthrow("splineinterpolate: Out of bounds",1);
      else if (ep == zeropad) { *deriv=0.0; extrapval = static_cast<T>(0.0); return(extrapval); }
      else if (ep == constpad) { *deriv=0.0; extrapval = padvalue; return(extrapval); }
    }

    T         partial = static_cast<T>(0.0);
    float     rval = 0.0;
    if (!splint.Valid() || getsplineorder() != splint.Order() || translate_extrapolation_type(getextrapolationmethod()) != splint.Extrapolation(0)) {
      forcesplinecoefcalculation();
    }
    rval = static_cast<float>(splint(x,y,z,dir,&partial));
    *deriv = static_cast<float>(partial);
    return(rval);
  }

  template<class T>
  float volume<T>::spline_interp3partial(// Input
                                         float x, float y, float z,              // Co-ordinates to get value for
                                         // Output
                                         float *dfdx, float *dfdy, float *dfdz)  // Partials
  const
  {
    if (!in_bounds(x,y,z)) {
      extrapolation ep = getextrapolationmethod();
      if (ep == boundsassert) { *dfdx=0.0; *dfdy=0.0; *dfdz=0.0; assert(false); extrapval = padvalue; return(extrapval); }
      else if (ep == boundsexception) imthrow("splineinterpolate: Out of bounds",1);
      else if (ep == zeropad) { *dfdx=0.0; *dfdy=0.0; *dfdz=0.0; extrapval = static_cast<T>(0.0); return(extrapval); }
      else if (ep == constpad) { *dfdx=0.0; *dfdy=0.0; *dfdz=0.0; extrapval = padvalue; return(extrapval); }
    }

    static std::vector<T>   partials(3,0);
    float                   rval = 0.0;
    if (!splint.Valid() || getsplineorder() != splint.Order() || translate_extrapolation_type(getextrapolationmethod()) != splint.Extrapolation(0)) {
      forcesplinecoefcalculation();
    }
    rval = static_cast<float>(splint.ValAndDerivs(x,y,z,partials));
    *dfdx = static_cast<float>(partials[0]);
    *dfdy = static_cast<float>(partials[1]);
    *dfdz = static_cast<float>(partials[2]);
    return(rval);
  }

  // TODO!!!!!!!! NEED TO PUT BACK A CACHE VALIDITY FLAG (LAZY-LITE) FOR SPLINES

  template<class T>
  float volume<T>::splineinterpolate(float x, float y, float z) const
  {
    this->throwsIfNot3D();
    extrapolation ep = getextrapolationmethod();
    if (!in_bounds(x,y,z)) {
      if (ep == boundsassert) { assert(false); extrapval = padvalue; return(extrapval); }
      else if (ep == boundsexception) imthrow("splineinterpolate: Out of bounds",1);
      else if (ep == zeropad) { extrapval = static_cast<T>(0.0); return(extrapval); }
      else if (ep == constpad) { extrapval = padvalue; return(extrapval); }
    }
    if (ep == extraslice) if (!in_extraslice_bounds(x,y,z)) { extrapval = padvalue; return(extrapval); }
    if (!splint.Valid() || getsplineorder() != splint.Order() || translate_extrapolation_type(ep) != splint.Extrapolation(0)) {
      forcesplinecoefcalculation();
    }
    return(static_cast<float>(splint(x,y,z)));
  }

  template <class T>
  float volume<T>::interpolate(float x, float y, float z) const
    {
      this->throwsIfNot3D();
      int64_t ix, iy, iz;
      switch (p_interpmethod) {
      case userinterpolation:
	if (p_userinterp == 0) {
	  imthrow("No user interpolation method set",7);
	} else {
	  return (*p_userinterp)(*this,x,y,z);
	}
      case nearestneighbour:
	ix=MISCMATHS::round(x); iy=MISCMATHS::round(y); iz=MISCMATHS::round(z);
	return this->operator()(ix,iy,iz);
      case trilinear:
	{
	  ix=(int64_t) floor(x); iy=(int64_t) floor(y); iz=(int64_t) floor(z);
	  if (in_neigh_bounds(*this,ix,iy,iz)) return interpolatevalue(x,y,z);
	  float dx=x-ix, dy=y-iy, dz=z-iz;
	  float v000=0, v001=0, v010=0, v011=0, v100=0, v101=0, v110=0, v111=0;
	  v000 = (float) this->operator()(ix,iy,iz);
	  v001 = (float) this->operator()(ix,iy,iz+1);
	  v010 = (float) this->operator()(ix,iy+1,iz);
	  v011 = (float) this->operator()(ix,iy+1,iz+1);
	  v100 = (float) this->operator()(ix+1,iy,iz);
	  v101 = (float) this->operator()(ix+1,iy,iz+1);
	  v110 = (float) this->operator()(ix+1,iy+1,iz);
	  v111 = (float) this->operator()(ix+1,iy+1,iz+1);
	  return q_tri_interpolation(v000,v001,v010,v011,v100,v101,v110,v111,
				     dx,dy,dz);
	}
      case sinc:
      case userkernel:
	{
	  return kernelinterpolation(x,y,z);
	}
      case spline:
        {
          return(splineinterpolate(x,y,z));
	}
      default:
	imthrow("Invalid interpolation method",6);
      }
      return 0.0;  // Should never get to here
    }


  template <class T>
  float volume<T>::interpolatevalue(float x, float y, float z) const
    {
      this->throwsIfNot3D();
      int64_t ix, iy, iz;
      switch (p_interpmethod) {
      case userinterpolation:
	if (p_userinterp == 0) {
	  imthrow("No user interpolation method set",7);
	} else {
	  return (*p_userinterp)(*this,x,y,z);
	}
      case nearestneighbour:
	ix=MISCMATHS::round(x); iy=MISCMATHS::round(y); iz=MISCMATHS::round(z);
	return value(ix,iy,iz);
      case trilinear:
	{
	  ix=(int64_t) floor(x); iy=(int64_t) floor(y); iz=(int64_t) floor(z);
	  float dx=x-ix, dy=y-iy, dz=z-iz;
	  T t000=0, t001=0, t010=0, t011=0, t100=0, t101=0, t110=0, t111=0;
	  float v000, v001, v010, v011, v100, v101, v110, v111;
	  this->getneighbours(ix,iy,iz,t000,t001,t010,t011,t100,t101,t110,t111);
	  v000=(float) t000; v001=(float) t001; v010=(float) t010;
	  v011=(float) t011; v100=(float) t100; v101=(float) t101;
	  v110=(float) t110; v111=(float) t111;
	  return q_tri_interpolation(v000,v001,v010,v011,v100,v101,v110,v111,
				     dx,dy,dz);
	}
      case sinc:
      case userkernel:
	{
	  return kernelinterpolation(x,y,z);
	}
      case spline:
        {
          return(splineinterpolate(x,y,z));
	}
      default:
	imthrow("Invalid interpolation method",6);
      }
      return 0.0;  // Should never get to here
    }

  template <class T>
  ColumnVector volume<T>::ExtractRow(int j, int k) const
  {
    if (j<0 || j>ysize()-1 || k<0 || k>zsize()-1) imthrow("ExtractRow: index out of range",3);
    ColumnVector rval(xsize());
    for (int i=0; i<xsize(); i++) rval(i+1) = (*this)(i,j,k);
    return(rval);
  }

  template <class T>
  ColumnVector volume<T>::ExtractColumn(int i, int k) const
  {
    if (i<0 || i>xsize()-1 || k<0 || k>zsize()-1) imthrow("ExtractColumn: index out of range",3);
    ColumnVector rval(ysize());
    for (int j=0; j<ysize(); j++) rval(j+1) = (*this)(i,j,k);
    return(rval);
  }

  template <class T>
  void volume<T>::SetRow(int j, int k, const ColumnVector& row)
  {
    if (j<0 || j>ysize()-1 || k<0 || k>zsize()-1) imthrow("SetRow: index out of range",3);
    if (row.Nrows() != xsize()) imthrow("SetRow: mismatched row vector",3);
    for (int i=0; i<xsize(); i++) (*this)(i,j,k) = row(i+1);
  }

  template <class T>
  void volume<T>::SetColumn(int i, int k, const ColumnVector& col)
  {
    if (i<0 || i>xsize()-1 || k<0 || k>zsize()-1) imthrow("SetColumn: index out of range",3);
    if (col.Nrows() != ysize()) imthrow("SetRow: mismatched row vector",3);
    for (int j=0; j<ysize(); j++) (*this)(i,j,k) = col(j+1);
  }


  int64_t mirrorclamp(int64_t x, int64_t x1, int64_t x2) {
    if (x2<x1) return mirrorclamp(x,x2,x1);
    if (x1==x2) return x1;
    int64_t x3 = 2*x2 - x1 + 1;
    int64_t nx = periodicclamp(x,x1,x3);
    if (nx > x2)
      nx = 2*x2 + 1 - nx;
    return nx;
  }

  template <class T>
  const T& volume<T>::extrapolate(int64_t x, int64_t y, int64_t z) const
  {
    switch (getextrapolationmethod()) {
    case userextrapolation:
      if (p_userextrap == 0) {
        imthrow("No user extrapolation method set",7);
      } else {
        extrapval = (*p_userextrap)(*this,x,y,z);
        return extrapval;
      }
    case zeropad:
      extrapval = (T) 0;
      return extrapval;
    case constpad:
      extrapval = padvalue;
      return extrapval;
    default:
      ; // do nothing
    }
    this->throwsIfNot3D();
    int64_t nx=x, ny=y, nz=z;
    switch (getextrapolationmethod()) {
    case periodic:
      nx = periodicclamp(x,0,maxx());
      ny = periodicclamp(y,0,maxy());
      nz = periodicclamp(z,0,maxz());
      return value(nx,ny,nz);
    case mirror:
      nx = mirrorclamp(x,0,maxx());
      ny = mirrorclamp(y,0,maxy());
      nz = mirrorclamp(z,0,maxz());
      return value(nx,ny,nz);
    case extraslice:
      if (nx==-1) { nx=0; }
      else { if (nx==xsize()) nx=maxx(); }
      if (ny==-1) { ny=0; }
      else { if (ny==ysize()) ny=maxy(); }
      if (nz==-1) { nz=0; }
      else { if (nz==zsize()) nz=maxz(); }
      if (in_bounds(nx,ny,nz)) { return value(nx,ny,nz); }
      else { extrapval = padvalue; return extrapval; }
    case boundsexception:
      if (!in_bounds(x,y,z)) {
        ostringstream msg;
        msg << "Out of Bounds at ("<<x<<","<<y<<","<<z<<")";
        imthrow(msg.str(),1);
      } else {
        return extrapval;
      }
    case boundsassert:
      assert(in_bounds(x,y,z));
      return extrapval;
    default:
      imthrow("Invalid extrapolation method",6);
    }

    return extrapval;
  }


  template <class T>
  void volume<T>::defineuserinterpolation(float (*interp)(
                         const volume<T>& , float, float, float)) const
  {
    p_userinterp = interp;
  }


  template <class T>
  void volume<T>::defineuserextrapolation(T (*extrap)(
                         const volume<T>& , int64_t, int64_t, int64_t)) const
  {
    p_userextrap = extrap;
  }


  template <class T>
  void volume<T>::definekernelinterpolation(const ColumnVector& kx,
					    const ColumnVector& ky,
					    const ColumnVector& kz,
					    int wx, int wy, int wz) const
  {
    // takes full-widths and converts all to half-widths
    int hwx = (wx-1)/2;
    int hwy = (wy-1)/2;
    int hwz = (wz-1)/2;
    interpkernel.setkernel(kx,ky,kz,hwx,hwy,hwz);
  }

  template <class T>
  void volume<T>::definekernelinterpolation(const volume<T>& vol) const
  {
    // copying like this is safe
    interpkernel = vol.interpkernel;
  }

  // Support Functions

  template <class T>
  void volume<T>::definesincinterpolation(const string& sincwindowtype,
					  int w, int nstore) const
  {
    // full width
    this->definesincinterpolation(sincwindowtype,w,w,w,nstore);
  }

  template <class T>
  void volume<T>::definesincinterpolation(const string& sincwindowtype,
					  int wx, int wy, int wz,
					  int nstore) const
  {
    // full widths
    if (nstore<1) nstore=1;
    ColumnVector kx, ky, kz;
    // calculate kernels
    kx = sinckernel1D(sincwindowtype,wx,nstore);
    ky = sinckernel1D(sincwindowtype,wy,nstore);
    kz = sinckernel1D(sincwindowtype,wz,nstore);

    this->definekernelinterpolation(kx,ky,kz,wx,wy,wz);
  }


  template<class T>
  void volume<T>::forcesplinecoefcalculation() const
  { //TODO eventually modify the spliterpolator class to use
    //cout << "forcesplinecoefcalculation: " << getsplineorder() << " " << splint.Order() << " " << splint.Valid() << endl;
    this->throwsIfNot3D();
    //if (splint.Valid() && -getsplineorder() == splint.Order() )
    //  imthrow("forcesplinecoefcalculation: requested spline for incoherent image",1);
    std::vector<unsigned int>                        dim(3,0);
    dim[0] = (unsigned int)xsize(); dim[1] = (unsigned int)ysize(); dim[2] =  (unsigned int)zsize();
    std::vector<SPLINTERPOLATOR::ExtrapolationType>  ep(3,SPLINTERPOLATOR::Mirror);
    for (unsigned int i=0; i<3; i++) ep[i] = translate_extrapolation_type(getextrapolationmethod());
    setsplineorder(abs(getsplineorder()));
    splint = SPLINTERPOLATOR::Splinterpolator<T> (fbegin(),dim,ep,getsplineorder(),false);
  }

  SPLINTERPOLATOR::ExtrapolationType translate_extrapolation_type(extrapolation ep)
  {
    switch (ep) {
    case zeropad:
      return(SPLINTERPOLATOR::Zeros);
      break;
    case extraslice:
      return(SPLINTERPOLATOR::Constant);  // It is constant for a given column, hence name.
      break;
    case mirror:
      return(SPLINTERPOLATOR::Mirror);
      break;
    case periodic:
      return(SPLINTERPOLATOR::Periodic);
      break;
    case boundsassert: case boundsexception: // We deal with this at the actual interpolation, and for now just return something
      return(SPLINTERPOLATOR::Zeros);
      break;
    case constpad: // Not implemented in splinterpolator, so I'll deal with this too at the actual interpolation.
      return(SPLINTERPOLATOR::Zeros);
      break;
    case userextrapolation:
      imthrow("translate_extrapolation_type: userextrapolation not implemented for spline interpolation",10);
      break;
    default:
      imthrow("translate_extrapolation_type: I am lost",10);
      break;
    }
    return(SPLINTERPOLATOR::Zeros);
  }





  // PROPERTIES


  template <class T>
  Matrix volume<T>::sampling_mat() const
  {
    Matrix samp=IdentityMatrix(4);
    samp(1,1) = xdim();
    samp(2,2) = ydim();
    samp(3,3) = zdim();
    // NOTE: no origin information is contained in this matrix!
    return samp;
  }

  template <class T>
  void volume<T>::set_sform(int sform_code, const Matrix& snewmat)
  {
    StandardSpaceTypeCode = sform_code;
    StandardSpaceCoordMat = snewmat;
  }


  template <class T>
  void volume<T>::set_qform(int qform_code, const Matrix& qnewmat)
  {
    RigidBodyTypeCode = qform_code;
    RigidBodyCoordMat = qnewmat;
  }


  template <class T>
  float volume<T>::intent_param(int n) const
  {
    float retval=0;
    if (n==1) { retval = IntentParam1; }
    if (n==2) { retval = IntentParam2; }
    if (n==3) { retval = IntentParam3; }
    return retval;
  }


  template <class T>
  void volume<T>::set_intent(int intent_code, float p1, float p2, float p3) const
  {
    IntentCode = intent_code;
    IntentParam1 = p1;
    IntentParam2 = p2;
    IntentParam3 = p3;
  }


  template <class T>
  ColumnVector volume<T>::principleaxis(int n) const
  {
    Matrix tmp = calc_principleaxes(*this);
    ColumnVector res = tmp.SubMatrix(1,3,n,n);
    return res;
  }

  template <class T>
  Matrix volume<T>::principleaxes_mat() const
  {
    return calc_principleaxes(*this);
  }


  int pval_index_end() { return -1; }

  template <class T>
  int get_pval_index(const std::vector<T>& pvals, float p)
  {
    int idx=0;
    while (idx < (int) pvals.size()) {
      // success if p is near pvals[idx] by a relative factor of 0.001 or less
      if ( fabs((p-pvals[idx])/Max(1e-5,Min(pvals[idx],1-pvals[idx]))) < 0.001 )
	return idx;
      else
	idx++;
    }
    return pval_index_end();
  }


  template <class T>
  T volume<T>::percentile(float pvalue, const volume<T>& mask) const
  {
    if ((pvalue>1.0) || (pvalue<0.0))
      { imthrow("Percentiles must be in the range [0.0,1.0]",4); }
    std::vector<float> pvaluevec;
    std::vector<T> retval;
    pvaluevec.push_back(pvalue);
    retval = calc_percentiles(*this,mask,pvaluevec);
    return retval[0];
  }


  template <class T>
  T volume<T>::percentile(float pvalue) const
  {
    if ((pvalue>1.0) || (pvalue<0.0))
      { imthrow("Percentiles must be in the range [0.0,1.0]",4); }
    std::vector<float> pvaluevec;
    pvaluevec.push_back(pvalue);
    return calc_percentiles(*this,pvaluevec)[0];
  }



  template <class T>
  std::vector<T> percentile_vec(std::vector<T>& hist,
				const std::vector<float>& percentilepvals)
  {
    unsigned int numbins = hist.size();
    if (numbins==0) {
      hist.push_back((T) 0);
      return hist;
    }

    sort(hist.begin(),hist.end());

    std::vector<T> outputvals(percentilepvals.size());
    for (unsigned int n=0; n<percentilepvals.size(); n++) {
      unsigned int percentile =
	(unsigned int) (((float) numbins) * percentilepvals[n]);
      if (percentile>=numbins)  percentile=numbins-1;
      outputvals[n] = hist[percentile];
    }
    return outputvals;
  }


  template < class T, class V>
  std::vector<T> calc_percentiles(const volume<T>& vol, const volume<V>& mask,
				  const std::vector<float>& percentilepvals)
  {
    bool usingMask( mask.totalElements() != 0 );
    if ( usingMask && !samesize(vol,mask,SUBSET) )
      imthrow("mask and vol have different sizes in calc_percentiles",3);
    std::vector<T> hist;
    if (!usingMask) hist.reserve(vol.totalElements());

    typename volume<V>::fast_const_iterator maskIt( mask.fbegin() ), maskEnd( mask.fend() );
    for ( typename volume<T>::fast_const_iterator it=vol.fbegin(), itEnd=vol.fend(); it < itEnd; it++ ) {
      if ( !usingMask || *(maskIt++) > 0.5 ) {
	hist.push_back(*it);
      }
      if ( usingMask && maskIt == maskEnd )
	maskIt=mask.fbegin(); //reset mask for 3D mask applied to a 4D volume
    }
    return percentile_vec(hist,percentilepvals);
  }


  template <class T>
  std::vector<T> calc_percentiles(const volume<T>& vol, const std::vector<float>& percentilepvals)
  {
    return calc_percentiles(vol,volume<char>(),percentilepvals);
  }


  template <class T>
  ColumnVector volume<T>::histogram(int nbins, T minval, T maxval) const
  {
    ColumnVector hist;
    imageHistogram(*this,nbins,minval,maxval,hist);
    return hist;
  }


  template <class T>
  ColumnVector volume<T>::histogram(int nbins) const
  {
    return histogram(nbins,robustmin(),robustmax());
  }

  template <class T>
  ColumnVector volume<T>::histogram(int nbins, T minval, T maxval,
				    const volume<T>& mask) const
  {
    ColumnVector hist;
    imageHistogram(*this,nbins,minval,maxval,hist,mask);
    return hist;
  }

  template <class T>
  ColumnVector volume<T>::histogram(int nbins, const volume<T>& mask) const
  {
    return histogram(nbins,robustmin(),robustmax(),mask);
  }

  template <class T>
  double volume<T>::mean(const volume<T>& mask) const
  {
    return sum(mask)/(Max((double) no_mask_voxels(mask)*tsize(),1.0));
  }


  template <class T>
  double volume<T>::variance(const volume<T>& mask) const
  {
    if (no_mask_voxels(mask)>0) {
      double n=(double) no_mask_voxels(mask)*tsize();
      return Max ( 0 , (n/Max(1.0,n-1))*(sumsquares(mask)/n - mean(mask)*mean(mask)) );
    } else {
      cerr << "ERROR:: Empty mask image" << endl;
      return 0;
    }
  }



  template <class T>
  double volume<T>::sum(const volume<T>& mask) const
  {
    return calculateSums(*this, mask)[0];
  }

  template <class T>
  double volume<T>::sumsquares(const volume<T>& mask) const
  {
    return calculateSums(*this, mask)[1];
  }

  template <class T>
  ColumnVector volume<T>::cog(const string& coordtype) const
  {
    // for coordtype="scaled_mm" return the old style, otherwise
    //  return newimage voxel coordinates
    ColumnVector retcog;
    retcog = calc_cog(*this);
    if (coordtype=="scaled_mm") {
      ColumnVector v(4);
      v << retcog(1) << retcog(2) << retcog(3) << 1.0;
      v = this->sampling_mat() * v;
      retcog(1) = v(1); retcog(2) = v(2); retcog(3) = v(3);
    }
    return retcog;
  }

  // the following calculates a robust background by taking the 10th percentile
  //  of the edge voxels only (from the first 3D volume *only*)
  template <class T>
  T calc_bval(const volume<T>& vol, int64_t edgewidth)
    {
      int64_t zb = vol.zsize(), yb = vol.ysize(), xb = vol.xsize();
      int64_t ewx, ewy, ewz, numbins;
      ewx = edgewidth;  ewy = edgewidth;  ewz = edgewidth;
      if (ewx >= xb)  ewx=xb-1;
      if (ewy >= yb)  ewy=yb-1;
      if (ewz >= zb)  ewz=zb-1;
      numbins = 2*(xb-2*ewx)*(yb-2*ewy)*ewz + 2*(xb-2*ewx)*zb*ewy + 2*yb*zb*ewx;
      std::vector<T> hist(numbins);
      // put the edge voxel values into the histogram
      int64_t hindx = 0;
      // put in the faces
      // xy faces (small lids of the box)
      for (int64_t e=0; e<ewz; e++) {
	for (int64_t x=ewx; x<xb-ewx; x++) {
	  for (int64_t y=ewy; y<yb-ewy; y++) {
	    hist[hindx++] = vol.value(x,y,e);
	    hist[hindx++] = vol.value(x,y,zb-1-e);
	  }
	}
      }
      // xz faces (smallish edge faces)
      for (int64_t e=0; e<ewy; e++) {
	for (int64_t x=ewx; x<xb-ewx; x++) {
	  for (int64_t z=0; z<zb; z++) {
	    hist[hindx++] = vol.value(x,e,z);
	    hist[hindx++] = vol.value(x,yb-1-e,z);
	  }
	}
      }
      // yz faces (large edge faces)
      for (int64_t e=0; e<ewx; e++) {
	for (int64_t y=0; y<yb; y++) {
	  for (int64_t z=0; z<zb; z++) {
	    hist[hindx++] = vol.value(e,y,z);
	    hist[hindx++] = vol.value(xb-1-e,y,z);
	  }
	}
      }
      sort(hist.begin(),hist.end());
      int64_t percentile10 = numbins / 10;
      T v10 = hist[percentile10];
      return v10;
    }

  template <class T>
  T calc_backgroundval(const volume<T>& vol)
  {
    return calc_bval(vol,2);
  }


  // Calculate on first 3D volume *only*
  template <class T>
  ColumnVector calc_cog(const volume<T>& vol)
    {
      ColumnVector v_cog(3);
      v_cog(1)=0.0;
      v_cog(2)=0.0;
      v_cog(3)=0.0;
      double val=0, total=0, vx=0, vy=0, vz=0, tot=0;
      T vmin=vol.min();
      int64_t n=0, nlim;
      nlim = (int64_t) sqrt((double) vol.nvoxels());
      if (nlim<1000) nlim=1000;
      for (int64_t z=0; z<vol.zsize(); z++) {
	for (int64_t y=0; y<vol.ysize(); y++) {
	  for (int64_t x=0; x<vol.xsize(); x++) {
	    val = (double) (vol(x,y,z) - vmin);
	    vx += val*x;
	    vy += val*y;
	    vz += val*z;
	    tot += val;
	    n++;
	    if (n>nlim) {
	      n=0; total+=tot; v_cog(1)+=vx; v_cog(2)+=vy; v_cog(3)+=vz;
	      tot=0; vx=0; vy=0; vz=0;
	    }
	  }
	}
      }
      total+=tot; v_cog(1)+=vx; v_cog(2)+=vy; v_cog(3)+=vz;
      if (fabs(total) < 1e-5) {
	cerr << "WARNING::in calculating COG, total = 0.0" << endl;
	total = 1.0;
      }
      v_cog(1) /= total;
      v_cog(2) /= total;
      v_cog(3) /= total;
      // Leave these values in (newimage) voxel coordinates
      return v_cog;
    }


  // Calculate on first 3D volume *only*
  template <class T>
  Matrix calc_principleaxes(const volume<T>& vol)
    {
      SymmetricMatrix m2(3);
      m2 = 0;
      double val=0, total=0, tot=0;
      double mxx=0, mxy=0, mxz=0, myy=0, myz=0, mzz=0, mx=0, my=0, mz=0;
      ColumnVector mean(3);
      mean = 0;
      T vmin=vol.min();

      int64_t n=0, nlim;
      nlim = (int64_t) sqrt((double) vol.nvoxels());
      if (nlim<1000) nlim=1000;
      for (int64_t z=0; z<vol.zsize(); z++) {
	for (int64_t y=0; y<vol.ysize(); y++) {
	  for (int64_t x=0; x<vol.xsize(); x++) {
	    val = (double) (vol(x,y,z) - vmin);
	    mxx += val*x*x;
	    mxy += val*x*y;
	    mxz += val*x*z;
	    myy += val*y*y;
	    myz += val*y*z;
	    mzz += val*z*z;
	    mx += val*x;
	    my += val*y;
	    mz += val*z;
	    tot += val;
	    n++;
	    if (n>nlim) {
	      n=0; total+=tot; m2(1,1)+=mxx; m2(1,2)+=mxy; m2(1,3)+=mxz;
	      m2(2,2)+=myy; m2(2,3)+=myz; m2(3,3)+=mzz;
	      mean(1)+=mx; mean(2)+=my; mean(3)+=mz;
	      tot=0; mxx=0; mxy=0; mxz=0; myy=0; myz=0; mzz=0; mx=0; my=0; mz=0;
	    }
	  }
	}
      }
      total+=tot; m2(1,1)+=mxx; m2(1,2)+=mxy; m2(1,3)+=mxz;
      m2(2,2)+=myy; m2(2,3)+=myz; m2(3,3)+=mzz;
      mean(1)+=mx; mean(2)+=my; mean(3)+=mz;

      if (fabs(total) < 1e-5) {
	cerr << "WARNING::in calculating Principle Axes, total = 0.0" << endl;
	total = 1.0;
      }
      m2 /= total;
      mean /= total;
      // Now adjust for voxel dimensions
      Matrix samp(3,3);
      samp = vol.sampling_mat().SubMatrix(1,3,1,3);
      m2 << samp * m2 * samp;
      mean = samp*mean;
      // Now make it central (taking off the cog)
      Matrix meanprod(3,3);
      for (int k1=1; k1<=3; k1++) {
	for (int k2=1; k2<=3; k2++) {
	  meanprod(k1,k2) = mean(k1)*mean(k2);
	}
      }
      m2 << m2 - meanprod;

      Matrix paxes;
      DiagonalMatrix evals;
      Jacobi(m2,evals,paxes);
      // Force the eigenvalues (and vectors) to be in descending order
      ColumnVector ptemp;
      float etemp;
      // brute force sort the eigen values and vectors
      // find inedx of least e-value
      int kmin=1;
      for (int k=2; k<=3; k++) {
	if (evals(k,k) < evals(kmin,kmin)) kmin = k;
      }
      // put the least in position 1
      etemp = evals(1,1);
      ptemp = paxes.SubMatrix(1,3,1,1);
      evals(1,1) = evals(kmin,kmin);
      paxes.SubMatrix(1,3,1,1) = paxes.SubMatrix(1,3,kmin,kmin);
      evals(kmin,kmin) = etemp;
      paxes.SubMatrix(1,3,kmin,kmin) = ptemp;
      // check if remaining ones require swapping
      if (evals(3,3) < evals(2,2)) {
	etemp = evals(2,2);
	ptemp = paxes.SubMatrix(1,3,2,2);
	evals(2,2) = evals(3,3);
	paxes.SubMatrix(1,3,2,2) = paxes.SubMatrix(1,3,3,3);
	evals(3,3) = etemp;
	paxes.SubMatrix(1,3,3,3) = ptemp;
      }
      return paxes;
    }


  template <class T, class V>
  int64_t imageHistogram(const volume<T>& vol, const int bins,
			 const T min, const T max, ColumnVector& hist, const volume<V>& mask)
  {
    maskedIterator<T,V> it(vol,mask,"imageHistogram: ");
    if ( !it.isValid() )
      cerr << "findHistogram: mask is empty" << endl;
    if (max<=min) return -1;
    int64_t validsize(0);
    if (hist.Nrows()!=bins) hist.ReSize(bins);
    hist=0;    // zero histogram
    // create histogram; the MIN is so that the maximum value falls in the last valid bin, not the (last+1) bin
    double fA = ((double)bins)/(max-min);
    double fB = ( ((double)bins) * ((double)(-min)) ) / (max-min);
    for (; it.isValid(); ++it ) {
      ++hist(Max(0, Min( (int)(fA*(*it) + fB), bins-1) ) + 1);
      ++validsize;
    }
    return validsize;
  }

  template <class T>
  int64_t imageHistogram(const volume<T>& vol, int bins,
			 T& min, T& max, ColumnVector& hist)
  {
    return imageHistogram(vol,bins,min,max,hist,volume<char>());
  }

  template <class T, class S, class R>
  void find_thresholds(const S& vol, T& minval, T& maxval, const R& mask, bool use_mask=true)
  {
    // STEVE SMITH'S CODE (in ancient times) - adapted for newimage by MARK JENKINSON & MATTHEW WEBSTER
  int HISTOGRAM_BINS=1000;
  ColumnVector hist(HISTOGRAM_BINS);
  int MAX_PASSES=10;
  int top_bin=0, bottom_bin=0, count, pass=1,
    lowest_bin=0, highest_bin=HISTOGRAM_BINS-1;
  int64_t validsize;
  //vector<int64_t> coordinates;
  //vector<T> extrema(calculateExtrema(vol,coordinates,mask));
  //T thresh98=0, thresh2=0, min(extrema[1]), max(extrema[0]);

   T thresh98=0, thresh2=0, min, max;
  if (use_mask) { min=vol.min(mask), max=vol.max(mask); }
  else { min=vol.min();  max=vol.max(); }

  if (hist.Nrows()!=HISTOGRAM_BINS) { hist.ReSize(HISTOGRAM_BINS); }

  while ( (pass==1) ||
	  ( (double) (thresh98 - thresh2) < (((double) (max - min)) / 10.0) ) ) // test for very long tails
    // find histogram and thresholds
    {
      if (pass>1) // redo histogram with new min and max
	{
	  // increase range slightly from the 2-98% range found
	  bottom_bin=Max(bottom_bin-1,0);
	  top_bin=Min(top_bin+1,HISTOGRAM_BINS-1);

	  // now set new min and max on the basis of this new range
	  T tmpmin = (T)( min + ((double)bottom_bin/(double)(HISTOGRAM_BINS))*(max-min) );
	  max = (T)( min + ((double)(top_bin+1)/(double)(HISTOGRAM_BINS))*(max-min) );
	  min=tmpmin;
	}

      if (pass==MAX_PASSES || min==max)  // give up and revert to full range ...
	{
	  if (use_mask) { min=vol.min(mask); max=vol.max(mask); }
 	  else { min=vol.min();  max=vol.max(); }
	}

      if (use_mask) validsize = imageHistogram(vol,HISTOGRAM_BINS,min,max,hist,mask);
      else validsize = imageHistogram(vol,HISTOGRAM_BINS,min,max,hist);

      if (validsize<1)
	{
          minval=thresh2=min;
	  maxval=thresh98=max;
	  return;
	}

      if (pass==MAX_PASSES)  /* ... _but_ ignore end bins */
	{
	  validsize-= MISCMATHS::round(hist(lowest_bin+1)) +
	    MISCMATHS::round(hist(highest_bin+1));
	  lowest_bin++;
	  highest_bin--;
	}

      if (validsize<0) /* ie zero range */
	{
	  thresh2=thresh98=min;
	  break;
	}

      double fA = (max-min)/(double)(HISTOGRAM_BINS);

      for(count=0, bottom_bin=lowest_bin; count<validsize/50; bottom_bin++)
	count+=MISCMATHS::round(hist(bottom_bin + 1));
      bottom_bin--;
      thresh2 =  min + (T)((double)bottom_bin*fA);

      for(count=0, top_bin=highest_bin; count<validsize/50; top_bin--)
	count+=MISCMATHS::round(hist(top_bin + 1));
      top_bin++;
      thresh98 = min + (T)((double)(top_bin+1)*fA);

      if (pass==MAX_PASSES) break;
      pass++;
    }
    minval=thresh2;
    maxval=thresh98;

  }


  template <class T, class S>
  void find_thresholds(const S& vol, T& minval, T& maxval) {
    return find_thresholds(vol,minval,maxval,vol,false);
  }

  // Return to non-SS land




  template <class T>
  std::vector<T> calc_robustlimits(const volume<T>& vol)
  {
    std::vector<T> rlimits(2);
    T minval=0, maxval=0;
    find_thresholds(vol,minval,maxval);
    //find_robust_limits(vol,1000,hist,minval,maxval);  // MJ version
    rlimits[0] = minval;
    rlimits[1] = maxval;
    return rlimits;
  }

  template <class T>
  std::vector<T> calc_robustlimits(const volume<T>& vol, const volume<T>& mask)
  {
    std::vector<T> rlimits(2);
    if ( mask.nvoxels()>0 && no_mask_voxels(mask)==0) {
      cerr << "WARNING:: Empty mask image in calc_robustlimits" << endl;
      rlimits[0]=0;
      rlimits[1]=0;
      return rlimits;
    }
    T minval=0, maxval=0;
    find_thresholds(vol,minval,maxval,mask);
    rlimits[0] = minval;
    rlimits[1] = maxval;
    return rlimits;
  }

  template <class T>
  vector<T> volume<T>::robustlimits(const volume<T>& mask) const {
    return calc_robustlimits(*this,mask);
  }

  template <class T>
  T volume<T>::min() const {
    vector<int64_t> coords;
    return calculateExtrema(*this,coords,volume<char>())[1];
  }

  template <class T>
  T volume<T>::min(const volume<T>& mask) const {
    vector<int64_t> coords;
    return calculateExtrema(*this,coords,mask)[1];
  }

  template <class T>
  T volume<T>::min(const volume<T>& mask, vector<int64_t>& coords) const {
    return calculateExtrema(*this, coords, mask)[1];
  }

  template <class T>
  T volume<T>::max() const {
    vector<int64_t> coords;
    return calculateExtrema(*this,coords,volume<char>())[0];
  }

  template <class T>
  T volume<T>::max(const volume<T>& mask) const {
    vector<int64_t> coords;
    return calculateExtrema(*this,coords,mask)[0];
  }

  template <class T>
  T volume<T>::max(const volume<T>& mask, vector<int64_t>& coords) const {
    return calculateExtrema(*this, coords, mask)[0];
  }

  template <class T>
  double volume<T>::sum() const {
    return calculateSums(*this,volume<char>())[0];
  }
  template <class T>
  double volume<T>::sumsquares() const {
    return calculateSums(*this, volume<char>())[1];
  }

  template <class T>
  T volume<T>::robustmin() const { return calc_robustlimits(*this)[0]; }
  template <class T>
  T volume<T>::robustmax() const { return calc_robustlimits(*this)[1]; }
  template <class T>
  T volume<T>::backgroundval() const { return calc_backgroundval(*this); }


  template <class T>
  T volume<T>::robustmin(const volume<T>& mask) const
  {
    return calc_robustlimits(*this,mask)[0];
  }

  template <class T>
  T volume<T>::robustmax(const volume<T>& mask) const
  {
    return calc_robustlimits(*this,mask)[1];
  }



  // GENERAL MANIPULATION

  template <class T>
  void volume<T>::threshold(T lowerth, T upperth, threshtype tt)
  {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	if ( ! ( ((tt==inclusive) && ((*it)>= lowerth) && ((*it)<= upperth))
		 || ((tt==exclusive) && ((*it)> lowerth) && ((*it)< upperth)) ) )
	{
	  *it = 0;
	}
      }
  }



  template <class T>
  void volume<T>::binarise(T lowerth, T upperth, threshtype tt, bool invert)
  {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	if ( ((tt==inclusive) && ((*it)>= lowerth) && ((*it)<= upperth))
	  || ((tt==exclusive) && ((*it)> lowerth) && ((*it)< upperth)) )
	{
	  *it = !invert;
	} else {
	  *it = invert;
	}
      }
  }



  template <class S>
  inline S swapval(S xval, S yval, S zval, int dim)
  {
    switch (dim) {
    case 1:
      return xval;
    case -1:
      return -xval;
    case 2:
      return yval;
    case -2:
      return -yval;
    case 3:
      return zval;
    case -3:
      return -zval;
    }
    return (S) 0;  // should never get here
  }


  template <class T>
  inline int64_t coordval(const volume<T>& vol, int64_t x, int64_t y, int64_t z, int dim)
  {
    switch (dim) {
    case 1:
      return x;
    case -1:
      return vol.xsize() - x - 1;
    case 2:
      return y;
    case -2:
      return vol.ysize() - y - 1;
    case 3:
      return z;
    case -3:
      return vol.zsize() - z - 1;
    }
    return 0;  // should never get here
  }


  int dimarg(const string& val)
  {
    if (val=="x") {
      return 1;
    } else if (val=="x-" || val=="-x") {
      return -1;
    } else if (val=="y") {
      return 2;
    } else if (val=="y-" || val=="-y") {
      return -2;
    } else if (val=="z") {
      return 3;
    } else if (val=="z-" || val=="-z") {
      return -3;
    } else {
      return 0;
    }
  }


  template <class T>
  void setrow(Matrix& affmat, int rownum, int dimnum, const volume<T>& invol)
  {
    if (dimnum==1 || dimnum==-1) {
      affmat(rownum,1)=1*sign(dimnum); affmat(rownum,2)=0; affmat(rownum,3)=0;
    }
    if (dimnum==2 || dimnum==-2) {
      affmat(rownum,1)=0; affmat(rownum,2)=1*sign(dimnum); affmat(rownum,3)=0;
    }
    if (dimnum==3 || dimnum==-3) {
      affmat(rownum,1)=0; affmat(rownum,2)=0; affmat(rownum,3)=1*sign(dimnum);
    }
    if (dimnum>0) return;
    float fov=0.0;
    if (dimnum==-1) {
      fov = (invol.xsize() -1) * invol.xdim();
    }
    if (dimnum==-2) {
      fov = (invol.ysize() -1) * invol.ydim();
    }
    if (dimnum==-3) {
      fov = (invol.zsize() -1) * invol.zdim();
    }
    affmat(rownum,4)=fov;
  }


  template <class T>
  Matrix volume<T>::swapmat(int dim1, int dim2, int dim3) const
  {
    Matrix affmat(4,4);
    affmat = 0.0;
    affmat(4,4)=1.0;
    setrow(affmat,1,dim1,*this);
    setrow(affmat,2,dim2,*this);
    setrow(affmat,3,dim3,*this);
    return affmat;
  }


  template <class T>
  Matrix volume<T>::swapmat(const string& newx, const string& newy, const string& newz) const
  {
    return this->swapmat(dimarg(newx),dimarg(newy),dimarg(newz));
  }


  template <class T>
  void volume<T>::swapdimensions(const string& newx, const string& newy, const string& newz, const bool keepLRorder)
  {
    this->swapdimensions(dimarg(newx),dimarg(newy),dimarg(newz), keepLRorder);
  }


  template <class T>
  void volume<T>::swapdimensions(int dim1, int dim2, int dim3, bool keepLRorder)
  {
    basic_swapdimensions(dim1,dim2,dim3, keepLRorder);
  }


  template <class T>
  void volume<T>::basic_swapdimensions(int dim1, int dim2, int dim3, bool keepLRorder)
  {
    // valid entries for dims are +/- 1, 2, 3 (corresponding to +/- x,y,z)
    if ( (dim1>3) || (dim1<-3) || (dim1==0) ||
	 (dim2<-3) || (dim2>3) || (dim2==0) ||
	 (dim3<-3) || (dim3>3) || (dim3==0) )
      {
	imthrow("Invalid dimension numbers entered to swapdimensions",8);
      }

    if ( (std::abs(dim1)==std::abs(dim2)) || (std::abs(dim1)==std::abs(dim3))
	 || (std::abs(dim2)==std::abs(dim3)) )
      {
	imthrow("Dimension numbers were not a permutation in swapdimensions",8);
      }
    int64_t sx = std::abs(swapval(this->xsize(),this->ysize(),this->zsize(),dim1));
    int64_t sy = std::abs(swapval(this->xsize(),this->ysize(),this->zsize(),dim2));
    int64_t sz = std::abs(swapval(this->xsize(),this->ysize(),this->zsize(),dim3));
    volume<T> swapvol(sx,sy,sz);

    for(int64_t d7=0;d7<this->size7() && swapvol.totalElements() > 1;d7++)
      for(int64_t d6=0;d6<this->size6();d6++)
	for(int64_t d5=0;d5<this->size5();d5++)
	  for(int64_t t=0;t<this->tsize();t++) { //For each 3D subvolume copy across the individual voxels to their new location and then copy entire data block back
	    for (int64_t z=0; z<this->zsize(); z++)
	      for (int64_t y=0; y<this->ysize(); y++)
		for (int64_t x=0; x<this->xsize(); x++) {
		  int64_t nx = coordval(*this,x,y,z,dim1);
		  int64_t ny = coordval(*this,x,y,z,dim2);
		  int64_t nz = coordval(*this,x,y,z,dim3);
		  swapvol(nx,ny,nz)=this->value(x,y,z,t,d5,d6,d7);
		}
	    copy(swapvol.fbegin(),swapvol.fend(),&this->value(0,0,0,t,d5,d6,d7));
	  }




    // now fix up the spatial properties
    // if a LR flip has happened then for all properties retain the unflipped values
    // therefore the data really has flipped, as otherwise it views identically
    if (keepLRorder && (this->swapmat(dim1,dim2,dim3).Determinant() < 0)) {
      // arbitrarily choose x to flip (if necessary)
      dim1*=-1;
    }

    float dx = swapval(this->xdim(), this->ydim(), this->zdim(), dim1);
    float dy = swapval(this->xdim(), this->ydim(), this->zdim(), dim2);
    float dz = swapval(this->xdim(), this->ydim(), this->zdim(), dim3);

    Matrix swappedSamplingMatrix(IdentityMatrix(4));
    swappedSamplingMatrix(1,1)=fabs(dx);
    swappedSamplingMatrix(2,2)=fabs(dy);
    swappedSamplingMatrix(3,3)=fabs(dz);

    // fix sform and qform matrices
    //   NB: sform and qform are voxel->mm but swapmat is mm->mm (flirt mm),
    //   hence sampling mats
    Matrix newS = this->sform_mat() * this->sampling_mat().i() * this->swapmat(dim1,dim2,dim3).i() * swappedSamplingMatrix;
    Matrix newQ = this->qform_mat() * this->sampling_mat().i() * this->swapmat(dim1,dim2,dim3).i() * swappedSamplingMatrix;
    this->ColumnsX=sx;
    this->RowsY=sy;
    this->SlicesZ=sz;
    this->setdims(dx,dy,dz);
    this->set_sform(this->sform_code(), newS);
    this->set_qform(this->qform_code(), newQ);
  }

  template <class T>
  int volume<T>::left_right_order() const
  {
     /* Determines if the image is stored in neurological or radiological convention */
    int order=FSL_RADIOLOGICAL;
    float dets=-1.0, detq=-1.0, det=-1.0;
    if (qform_code()!=NIFTI_XFORM_UNKNOWN) {
      detq = qform_mat().Determinant();
      det = detq;
    }
    if (sform_code()!=NIFTI_XFORM_UNKNOWN) {
      dets = sform_mat().Determinant();
      det = dets;
    }

    if (det<0.0) order=FSL_RADIOLOGICAL;
    else order=FSL_NEUROLOGICAL;
    /* check for inconsistency if both are set */
    if ( (sform_code()!=NIFTI_XFORM_UNKNOWN) &&
	 (qform_code()!=NIFTI_XFORM_UNKNOWN) ) {
      if (dets * detq < 0.0) order=FSL_INCONSISTENT;
      if (fabs(dets * detq)<1e-12)  order=FSL_ZERODET;
    }
    if (fabs(det)<1e-12) order=FSL_ZERODET;
   return order;
  }


  template <class T>
  Matrix volume<T>::newimagevox2mm_mat() const
  {
    Matrix vox2mm_mat;
    short code=NIFTI_XFORM_UNKNOWN;  // not used currently, but it's nice!
    if (sform_code()!=NIFTI_XFORM_UNKNOWN) {
      vox2mm_mat = sform_mat();
      code = sform_code();
    } else if (qform_code()!=NIFTI_XFORM_UNKNOWN) {
      vox2mm_mat = qform_mat();
      code = qform_code();
    } else {
      /* default case - for FSLView is positive voxel to mm scalings */
      vox2mm_mat=IdentityMatrix(4);
      vox2mm_mat(1,1)=pixdim1();
      vox2mm_mat(2,2)=pixdim2();
      vox2mm_mat(3,3)=pixdim3();
      code=NIFTI_XFORM_UNKNOWN;
    }
    return vox2mm_mat;
  }


  template <class T>
  Matrix volume<T>::niftivox2newimagevox_mat() const
  {
    Matrix vox2vox=IdentityMatrix(4);
    if ((!RadiologicalFile) && (this->left_right_order()==FSL_RADIOLOGICAL)) {
      vox2vox = (this->sampling_mat()).i() * this->swapmat(-1,2,3) *
	(this->sampling_mat());
    }
    return vox2vox;
  }



  template <class T>
  void volume<T>::swapLRorder()
  {
    basic_swapdimensions(-1,2,3,false);
  }


  template <class T>
  void volume<T>::setLRorder(int LRorder)
  {
    if (LRorder != this->left_right_order()) { this->swapLRorder(); }
  }


  template <class T>
  void volume<T>::makeradiological()
  {
    // use existing matrices to determine the order and if necessary swap
    //  all data and matrices
    if (this->left_right_order()==FSL_NEUROLOGICAL) { this->swapLRorder(); }
  }


  template <class T>
  void volume<T>::makeneurological()
  {
    // use existing matrices to determine the order and if necessary swap
    //  all data and matrices
    if (this->left_right_order()==FSL_RADIOLOGICAL) { this->swapLRorder(); }
  }




  // ARITHMETIC OPERATIONS


  template <class T>
  T volume<T>::operator=(T val)
  {
    fill(nsfbegin(),nsfend(),val);  // use the STL
    return val;
  }

  template <class T>
  const volume<T>& volume<T>::operator+=(T val)
  {
    for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	 it!=itend; ++it) {
      *it += val;
    }
    return *this;
  }

  template <class T>
  const volume<T>& volume<T>::operator-=(T val)
  {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	*it -= val;
      }
    return *this;
  }

  template <class T>
  const volume<T>& volume<T>::operator*=(T val)
  {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	*it *= val;
      }
    return *this;
  }

  template <class T>
  const volume<T>& volume<T>::operator/=(T val)
  {
      for (nonsafe_fast_iterator it=nsfbegin(), itend=nsfend();
	   it!=itend; ++it) {
	*it /= val;
      }
    return *this;
  }


  template <class T>
  const volume<T>& volume<T>::operator+=(const volume<T>& source)
  {
    if ( ( this->dimensionality() < source.dimensionality() ) || !samesize(*this,source,SUBSET) )
      imthrow("Attempted to add images of different sizes",3);
    cyclicIterator<T> srcIt(source);
    for ( typename volume<T>::nonsafe_fast_iterator it=this->nsfbegin(), itEnd=this->nsfend(); it != itEnd; ++it, ++srcIt )
      *it += *srcIt;
    return *this;
  }


  template <class T>
  const volume<T>& volume<T>::operator-=(const volume<T>& source)
  {
    if ( ( this->dimensionality() < source.dimensionality() ) || !samesize(*this,source,SUBSET) )
      imthrow("Attempted to subtract images of different sizes",3);
    cyclicIterator<T> srcIt(source);
    for ( typename volume<T>::nonsafe_fast_iterator it=this->nsfbegin(), itEnd=this->nsfend(); it != itEnd; ++it, ++srcIt )
      *it -= *srcIt;
    return *this;
  }


  template <class T>
  const volume<T>& volume<T>::operator*=(const volume<T>& source)
  {
    if ( ( this->dimensionality() < source.dimensionality() ) || !samesize(*this,source,SUBSET) )
      imthrow("Attempted to multiply images of different sizes",3);
    cyclicIterator<T> srcIt(source);
    for ( typename volume<T>::nonsafe_fast_iterator it=this->nsfbegin(), itEnd=this->nsfend(); it != itEnd; ++it, ++srcIt )
      *it *= *srcIt;
    return *this;
  }


  template <class T>
  const volume<T>& volume<T>::operator/=(const volume<T>& source)
  {
    if ( ( this->dimensionality() < source.dimensionality() ) || !samesize(*this,source,SUBSET) )
      imthrow("Attempted to divide images of different sizes",3);
    cyclicIterator<T> srcIt(source);
    for ( typename volume<T>::nonsafe_fast_iterator it=this->nsfbegin(), itEnd=this->nsfend(); it != itEnd; ++it, ++srcIt )
      *it /= *srcIt;
    return *this;
  }


  template <class T>
  volume<T> volume<T>::operator+(T num) const
  {
    volume<T> tmp = *this;
    tmp+=num;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator-(T num) const
  {
    volume<T> tmp = *this;
    tmp-=num;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator*(T num) const
  {
    volume<T> tmp;
    tmp = *this;
    tmp*=num;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator/(T num) const
  {
    volume<T> tmp = *this;
    tmp/=num;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator+(const volume<T>& vol2) const
  {
    volume<T> tmp(this->dimensionality()>=vol2.dimensionality() ? *this: vol2);
    tmp+=this->dimensionality()<vol2.dimensionality() ? *this: vol2;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator-(const volume<T>& vol2) const
  {
    volume<T> tmp(this->dimensionality()>=vol2.dimensionality() ? *this: vol2);
    tmp-=this->dimensionality()<vol2.dimensionality() ? *this: vol2;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator*(const volume<T>& vol2) const
  {
    volume<T> tmp(this->dimensionality()>=vol2.dimensionality() ? *this: vol2);
    tmp*=this->dimensionality()<vol2.dimensionality() ? *this: vol2;
    return tmp;
  }

  template <class T>
  volume<T> volume<T>::operator/(const volume<T>& vol2) const
  {
    volume<T> tmp(this->dimensionality()>=vol2.dimensionality() ? *this: vol2);
    tmp/=this->dimensionality()<vol2.dimensionality() ? *this: vol2;
    return tmp;
  }


template <class T>
void volume<T>::copyROI(const volume<T>& source,const int64_t minx, const int64_t miny, const int64_t minz, const int64_t mint, const int64_t maxx, const int64_t maxy, const int64_t maxz, const int64_t maxt)
{
  if ( ( source.maxx() != maxx-minx ) || ( source.maxy() != maxy-miny )  || ( source.maxz() != maxz-minz )  || ( source.maxt() != maxt-mint ) )
      imthrow("copyROI source must match requested dimensions",3);

  // now copy only the appropriate data
  for ( int64_t t=mint;t<=maxt;t++)
    for (int64_t z=minz; z<=maxz; z++)
      for (int64_t y=miny; y<=maxy; y++)
	for (int64_t x=minx; x<=maxx; x++)
	  (*this)(x,y,z,t)=source(x - minx, y - miny, z - minz, t - mint);
}


template <class T>
volume<T> volume<T>::ROI(const int64_t minx, const int64_t miny, const int64_t minz, const int64_t mint, const int64_t maxx, const int64_t maxy, const int64_t maxz, const int64_t maxt) const
{
  volume<T> roivol;
  roivol.reinitialize(maxx-minx+1,maxy-miny+1,maxz-minz+1,maxt-mint+1);

  // now copy only the appropriate data
  for ( int64_t t=mint;t<=maxt;t++)
    for (int64_t z=minz; z<=maxz; z++)
      for (int64_t y=miny; y<=maxy; y++)
	for (int64_t x=minx; x<=maxx; x++)
	  roivol(x - minx, y - miny, z - minz, t - mint) = (*this)(x,y,z,t);

  roivol.copyproperties(*this);
  // set sform and qform matrices appropriately (if set)
  Matrix roi2vol= IdentityMatrix(4);
  roi2vol(1,4) = minx;
  roi2vol(2,4) = miny;
  roi2vol(3,4) = minz;
  if (this->sform_code()!=NIFTI_XFORM_UNKNOWN)
    roivol.set_sform(this->sform_code(),this->sform_mat() * roi2vol);
  if (this->qform_code()!=NIFTI_XFORM_UNKNOWN)
    roivol.set_qform(this->qform_code(),this->qform_mat() * roi2vol);
  return roivol;
}


  template <class T>
  ReturnMatrix volume<T>::voxelts(int64_t x, int64_t y, int64_t z) const
  {
    ColumnVector res;
    if (this->tsize()<1) return res;
    res.ReSize(this->tsize());
    for (int t=0; t<this->tsize(); t++)
      res(t + 1) = (NEWMAT::Real) this->value(x,y,z,t);
    res.Release();
    return res;
  }

  template <class T>
  void volume<T>::setvoxelts(const ColumnVector& ts, int64_t x, int64_t y, int64_t z)
  {
    if (ts.Nrows() != this->tsize())
      imthrow("setvoxelts - incorrectly sized vector",3);
    for (int t=0; t<this->tsize(); t++) {
      this->value(x,y,z,t) = (T) ts(t+1);
    }
  }


  template < class T,class V >
  vector<double> calculateSums(const volume<T>& inputVolume, const volume<V>& mask)
  {
    vector<double> sums(2,0); //sum/sumsq
    size_t n(0), nn(0), nlim( max( (long)sqrt( (double)inputVolume.totalElements() ) ,10000L) );
    double sum=0, sum2=0;
    for ( maskedIterator<T,V> it(inputVolume,mask,"calculateSums: ");it.isValid(); ++it ) {
        double value=(*it);
	sum += value;
	sum2 += value*value;
	n++;
	if (n>nlim) { nn++; sums[0]+=sum; sums[1]+=sum2; sum=sum2=n=0; }
    }
    if (n + nn == 0)
     cerr << "ERROR:: Empty mask image" << endl;
    sums[0]+=sum;
    sums[1]+=sum2;
    return sums;
  }

template < class T, class V>
vector<T> calculateExtrema(const volume<T>& inputVolume, vector<int64_t>& coordinates, const volume<V>& mask )
{
  maskedIterator<T,V> it(inputVolume,mask,"calculateExtrema: ");
  if ( !it.isValid() ) {
    if ( inputVolume.totalElements() ) {
      // throw runtime_error("calculateExtrema: mask is empty");
      cerr << "calculateExtrema: mask is empty" << endl;	 //cerr for compat with old newimage behaviour
    }
    return vector<T>(2,0);
  }
  vector<T>extrema(2,*it); //max/min
  size_t maxVoxel(it.offset()),minVoxel(it.offset());
  while((++it).isValid()) {
    if ( *it > extrema[0] ) {
      extrema[0]=*it;
      maxVoxel=it.offset();
    } else if ( *it < extrema[1] ) {
      extrema[1]=*it;
      minVoxel=it.offset();
    }
  }
  coordinates=inputVolume.ptrToCoord(maxVoxel);
  vector<int64_t> minCoordinates(inputVolume.ptrToCoord(minVoxel));
  coordinates.insert(coordinates.begin()+7,minCoordinates.begin(),minCoordinates.end());
  return extrema;
}


  // END OF CURRENT CHANGES




  // MATRIX <-> VOLUME CONVERSIONS

  template <class T>
  ReturnMatrix volume<T>::matrix(const volume<T>& mask, bool using_mask) const
  {
    Matrix matv;
    if (this->totalElements()==0) return matv;
    if (this->dimensionality()>4)
      imthrow("Cannot use matrix() with images higher than 4D",77);
    if (using_mask && !samesize(*this,mask,3))
      imthrow("Mask of different size used in matrix()",3);
    long cidx(1);
    long ncols=this->nvoxels();
    if (using_mask) ncols=no_mask_voxels(mask);
    matv.ReSize(this->tsize(), ncols);
    for (int64_t z = 0; z < zsize(); z++) {
      for (int64_t y = 0; y < ysize(); y++) {
	for (int64_t x = 0; x < xsize(); x++) {
	  if (!using_mask || (mask(x,y,z)>0)) {
	    for (int64_t t = 0; t < this->tsize(); t++) {
	      matv(t+1,cidx) = this->value(x,y,z,t);
	    }
	    cidx++;
	  }
	}
      }
    }
    matv.Release();
    return matv;
  }


  template <class T>
  ReturnMatrix volume<T>::matrix(const volume<T>& mask, vector<int64_t>& voxelLabels) const
  {
    //typename boost::enable_if<boost::is_integral<T>, T>::type foo;
    //BOOST_STATIC_ASSERT_MSG(boost::is_integral<S>::value, "matrix only specialised for long int and int64_t");


    voxelLabels.clear();
    Matrix matv;
    if (this->totalElements()==0) return matv;
    if (this->dimensionality()>4)
      imthrow("Cannot use matrix() with images higher than 4D",77);
    if (!samesize(*this,mask,SUBSET))
      imthrow("Mask of different size used in matrix()",3);
    long cidx (1);
    matv.ReSize(this->tsize(), no_mask_voxels(mask) );

    for (int64_t z=0; z<mask.zsize(); z++) {
      for (int64_t y=0; y<mask.ysize(); y++) {
	for (int64_t x=0; x<mask.xsize(); x++) {
	  if (mask(x,y,z)>maskThreshold()) {
	    voxelLabels.push_back(x+y*mask.xsize()+z*mask.xsize()*mask.ysize());
	    for (int64_t t=0; t<this->tsize(); t++) {
	      matv(t+1,cidx) = this->value(x,y,z,t);
	    }
	    cidx++;
	  }
	}
      }
    }
    matv.Release();
    return matv;
  }

  template <class T>
  ReturnMatrix volume<T>::matrix() const
  {
    volume<T> dummy;
    return matrix(dummy,false);
  }


  template <class T>
  void volume<T>::setmatrix(const Matrix& newmatrix, const volume<T>& mask,
			    const T pad, bool using_mask)
  {
    if (using_mask && (mask.dimensionality()>3))
      imthrow("Cannot use setmatrix() with masks with more than 3 dimensions",77);
    int tsz = this->tsize();
    if (using_mask) {
      if ( (tsz==0) || (tsz!=newmatrix.Nrows()) || (!samesize(mask,*this,3)) )
	this->reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),newmatrix.Nrows());
      this->copyproperties(mask);
    } else {
      if ( (tsz==0) || (tsz!=newmatrix.Nrows()) )
	imthrow("Mismatch in number of timepoints in setmatrix()",4);
    }

    this->operator=(pad);
    int64_t nvox=this->nvoxels();
    if (using_mask) { nvox=no_mask_voxels(mask); }
    if (using_mask && (newmatrix.Ncols()!=nvox)) {
      imthrow("Incompatible number of valid voxels and matrix columns in setmatrix()",4);
    }
    long cidx = 1;
    for (int64_t z=0; z<this->zsize(); z++) {
      for (int64_t y=0; y<this->ysize(); y++) {
	for (int64_t x=0; x<this->xsize(); x++) {
	  if (!using_mask || (mask(x,y,z)>0)) {
	    for (int64_t t=0; t<this->tsize(); t++) {
	      this->value(x,y,z,t) = (T) newmatrix(t+1,cidx);
	    }
	    cidx++;
	  }
	}
      }
    }

  }


  template <class T>
  void volume<T>::setmatrix(const Matrix& newmatrix)
  {
    volume<T> dummymask;
    this->setmatrix(newmatrix,dummymask,0,false);
  }

  template <class T>
  volume<int> volume<T>::vol2matrixkey(const volume<T>& mask)  // TODO - in future this should return volume<int64_t>, but will work for a while to come...
  {
    if (!samesize(*this,mask,3))
      imthrow("Mask of different size used in vol2matrixkey",3);
    int count=1;
    volume<int> tmp(this->xsize(),this->ysize(),this->zsize());
    copybasicproperties(*this,tmp);
    for(int64_t z=0;z< this->zsize();z++){
      for(int64_t y=0;y<this->ysize();y++){
	for(int64_t x=0;x<this->xsize();x++){
	  if (mask(x,y,z)>0){
	    tmp(x,y,z)=count;
	    count++;
	  }
	  else{
	    tmp(x,y,z)=0;
	  }
	}
      }
    }
    return tmp;
  }

  template <class T>
  ReturnMatrix volume<T>::matrix2volkey(volume<T>& mask){
    if (!samesize(*this,mask,3))
      imthrow("Mask of different size used in matrix2volkey",3);
    int64_t count=no_mask_voxels(mask);
    Matrix key(count,3);
    count=1;
    for(int64_t z=0;z< this->zsize();z++)
      for(int64_t y=0;y<this->ysize();y++)
	for(int64_t x=0;x<this->xsize();x++)
	  if(mask(x,y,z)>0){
	    key(count,1)=x;
	    key(count,2)=y;
	    key(count,3)=z;
	    count++;
	  }
    key.Release();
    return key;
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////


  // SPECIFIC INSTANTIATIONS

  // provide only these instances of the class
  //  NB: unsigned int is not included as it is not a valid AVW DT_TYPE
  template class volume<char>;
  template class volume<short>;
  template class volume<int>;
  template class volume<float>;
  template class volume<double>;
  // WORKING template ReturnMatrix volume<int>::matrix(const volume<int>& mask, vector<int64_t>& voxelLabels) const;



} // end namespace
