/*  General IO functions (images and transformation files)

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 2000-2012 University of Oxford  */

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


#if !defined(__newimageio_h)
#define __newimageio_h

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include "NewNifti/NewNifti.h"
#include "newmatio.h"
#include "newimage.h"
#include "complexvolume.h"
#include "miscmaths/miscmaths.h"

using namespace NEWMAT;

namespace NEWIMAGE {

bool FslIsSingleFileType(int filetype);
int FslNiftiVersionFileType(int filetype);
int FslFiletypeFromHeader(const NiftiHeader& header);
int fslFileType(string filename);
int FslGetEnvOutputType(void);
string outputExtension(const int filetype);

bool FslFileExists(const string& filename);
bool FslImageExists(const string& filename);
inline bool fsl_imageexists(const string& filename) { return FslImageExists(filename); }
bool FslIsCompressedFileType(int filetype);
string return_validimagefilename(const string& filename, const bool quiet=false, const bool strict=false);
string make_basename(string& filename);
string make_basename(const string& filename);
string appendFSLfilename(const string inputName, const string addendum);
int find_pathname(string& filename);
int fslFileType(string filename);
template <class T>
void ConvertAndScaleNewNiftiBuffer(char* buffer, T*& tbuffer, const NiftiHeader& niihdr, const size_t & imagesize);
  // read

template <class T>
int read_volume(volume<T>& target, const string& filename);
template <class T>
int read_volumeROI(volume<T>& target, const string& filename,
		   int64_t x0, int64_t y0, int64_t z0, int64_t x1, int64_t y1, int64_t z1);
template <class T>
int read_volumeROI(volume<T>& target, const string& filename,
		   int64_t x0, int64_t y0, int64_t z0, int64_t x1, int64_t y1, int64_t z1, 
		   int64_t xskip, int64_t yskip, int64_t zskip);
template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		     int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
		     int64_t x1, int64_t y1, int64_t z1, int64_t t1);

template <class T>
int read_volumeROI(volume<T>& target, const string& filename, short& dtype,
		   int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
		   int64_t x1, int64_t y1, int64_t z1, int64_t t1,
		   const bool swap2radiological, const bool readAs4D=false) {
  return read_volumeROI(target,filename,dtype,x0,y0,z0,t0,-1L,-1L,-1L,x1,y1,z1,t1,-1L,-1L,-1L,swap2radiological,readAs4D);
 }

template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		     int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
		     int64_t x1, int64_t y1, int64_t z1, int64_t t1,
		     const int64_t xskip, int64_t yskip, int64_t zskip, int64_t tskip);
 
int read_complexvolume(volume<float>& realvols, volume<float>& imagvols,
		       const string& filename,  bool read_img_data=true);
int read_complexvolume(complexvolume& vol, const string& filename);  

template <class T>
void set_volume_properties(const NiftiHeader& niihdr, volume<T>& target);


template <class T>
int read_volume_hdr_only(volume<T>& target, const string& filename) {
  target.destroy();
  NiftiIO Reader;
  NiftiHeader niihdr;
  try {
    niihdr = Reader.loadHeader(return_validimagefilename(filename));
  } catch ( exception& e ) { imthrow("Failed to read volume "+filename+"\nError : "+e.what(),22); }
  for (int n=1; n<=7; n++) {
    if (niihdr.dim[n]<1) niihdr.dim[n]=1;  // make it robust to dim[n]=0
  }
  T* tbuffer=new T[1];
  target.initialize(niihdr.dim[1],niihdr.dim[2],niihdr.dim[3],niihdr.dim[4],niihdr.dim[5],niihdr.dim[6],niihdr.dim[7],tbuffer,true);
  // copy info from file
  set_volume_properties(niihdr,target);
  return 0;
 }

template <class T>
int read_volume4D_hdr_only(volume<T>& target, const string& filename) {
  return read_volume_hdr_only(target,filename);
}


  // save

template <class T>
int save_volume(const volume<T>& source, const string& filename, const int filetype=-1);
  
int save_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename);
int save_complexvolume(const complexvolume& vol, const string& filename);
  




// Helper functions
short closestTemplatedType(const short inputType);

int read_volume_size(const string& filename, 
		     int64_t& sx, int64_t& sy, int64_t& sz, int64_t& st, int64_t& s5, int64_t& s6, int64_t& s7);

short dtype(const char* T);
short dtype(const short* T);
short dtype(const int* T);
short dtype(const float* T);
short dtype(const double* T);

short dtype(const volume<char>& vol);
short dtype(const volume<short>& vol);
short dtype(const volume<int>& vol);
short dtype(const volume<float>& vol);
short dtype(const volume<double>& vol);

short dtype(const string& filename);

// Boring overloads to enable different names (load and write)


// load

template <class T>
int load_volume(volume<T>& target, const string& filename);
  
template <class T>
int load_complexvolume(volume<float>& realvol, volume<float>& imagvol,
		       const string& filename);
int load_complexvolume(complexvolume& vol, const string& filename);


// write

template <class T>
int write_volume(const volume<T>& source, const string& filename);
 
int write_complexvolume(const volume<float>& realvol, 
		       const volume<float>& imagvol, const string& filename);
int write_complexvolume(const complexvolume& vol, const string& filename);
  
int write_complexvolume(const volume<float>& realvol, 
			 const volume<float>& imagvol, 
  			 const string& filename);
 
////////////////////////////////////////////////////////////////////////
///////////////////////// TEMPLATE DEFINITIONS /////////////////////////
////////////////////////////////////////////////////////////////////////


// External functions

// READ FUNCTIONS


template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		   int64_t x0, int64_t y0, int64_t z0, int64_t x1, int64_t y1, int64_t z1,
		   int64_t xskip, int64_t yskip, int64_t zskip)
{
  int retval=read_volumeROI(target,filename,x0,y0,z0,x1,y1,z1);
  if (retval==0) {
    if (xskip<1) xskip=1;    
    if (yskip<1) yskip=1;
    if (zskip<1) zskip=1;
    int64_t sx=(target.maxx()-target.minx())/xskip + 1;
    int64_t sy=(target.maxy()-target.miny())/yskip + 1;
    int64_t sz=(target.maxz()-target.minz())/zskip + 1;
    volume<T> tmpvol(sx,sy,sz);
    int64_t xx=0, yy=0, zz=0, x=0, y=0, z=0;
    for (z=target.minz(), zz=0; z<=target.maxz(); z+=zskip, zz++) {
      for (y=target.miny(), yy=0; y<=target.maxy(); y+=yskip, yy++) {
	for (x=target.minx(), xx=0; x<=target.maxx(); x+=xskip, xx++) {
	  tmpvol(xx,yy,zz) = target(x,y,z);
	}
      }
    }
    tmpvol.copyproperties(target);
    target = tmpvol;
  }
  return retval;
}






template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		     int64_t x0, int64_t y0, int64_t z0, 
		     int64_t x1, int64_t y1, int64_t z1)
{
  short dtype;
  return read_volumeROI(target,filename,dtype,
			x0,y0,z0,0,x1,y1,z1,0);
}


template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		     int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
		     int64_t x1, int64_t y1, int64_t z1, int64_t t1,
		     int64_t xskip, int64_t yskip, int64_t zskip, int64_t tskip)
{
  read_volumeROI(target,filename,x0,y0,z0,t0,x1,y1,z1,t1);
  if (xskip<1) xskip=1;    
  if (yskip<1) yskip=1;
  if (zskip<1) zskip=1;
  if (tskip<1) tskip=1;
  int64_t sx=(target.maxx()-target.minx())/xskip + 1;
  int64_t sy=(target.maxy()-target.miny())/yskip + 1;
  int64_t sz=(target.maxz()-target.minz())/zskip + 1;
  int64_t st=(target.maxt()-target.mint())/tskip + 1;
  volume<T> tmpvol(sx,sy,sz,st);
  int64_t xx=0, yy=0, zz=0, tt=0, x=0, y=0, z=0, t=0;
  for (t=target.mint(), tt=0; t<=target.maxt(); t+=tskip, tt++) {
    for (z=target.minz(), zz=0; z<=target.maxz(); z+=zskip, zz++) {
      for (y=target.miny(), yy=0; y<=target.maxy(); y+=yskip, yy++) {
	for (x=target.minx(), xx=0; x<=target.maxx(); x+=xskip, xx++) {
	  tmpvol(xx,yy,zz,tt) = target(x,y,z,t);
	}
      }
    }
  }
  tmpvol.copyproperties(target[0]);
  target = tmpvol;
  return 0;
}

template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		     int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
		     int64_t x1, int64_t y1, int64_t z1, int64_t t1)
{
  short dtype;
  int retval = read_volumeROI(target,filename,dtype,
			      x0,y0,z0,t0,x1,y1,z1,t1,true);
  return retval;
}

template <class T>
int read_volume(volume<T>& target, const string& filename, 
		  short& dtype)
{
  return read_volumeROI(target,filename,dtype,0,0,0,0,-1,-1,-1,-1,true);
}


template <class T>
int read_volume(volume<T>& target, const string& filename)
{
  short dtype;
  read_volume(target,filename,dtype);
  return 0;
}



// SAVE FUNCTIONS


mat44 newmat2mat44(const Matrix& nmat);

template <class T>
int set_fsl_hdr(const volume<T>& source, NiftiHeader& niihdr) 
{
  niihdr.dim[0]=source.dimensionality();
  niihdr.dim[1]=source.size1();
  niihdr.dim[2]=source.size2();
  niihdr.dim[3]=source.size3();
  niihdr.dim[4]=source.size4();
  niihdr.dim[5]=source.size5();
  niihdr.dim[6]=source.size6();
  niihdr.dim[7]=source.size7();

  niihdr.datatype=dtype(source);

  niihdr.pixdim[1]=source.pixdim1();
  niihdr.pixdim[2]=source.pixdim2();
  niihdr.pixdim[3]=source.pixdim3();
  niihdr.pixdim[4]=source.pixdim4();
  niihdr.pixdim[5]=source.pixdim5();
  niihdr.pixdim[6]=source.pixdim6();
  niihdr.pixdim[7]=source.pixdim7();


  niihdr.sformCode = source.sform_code();
  niihdr.qformCode = source.qform_code();
  niihdr.setSForm(newmat2mat44(source.sform_mat()));
  niihdr.setQForm(newmat2mat44(source.qform_mat()));

  niihdr.intentCode = source.intent_code();
  niihdr.intent_p1 = source.intent_param(1);
  niihdr.intent_p2 = source.intent_param(2);
  niihdr.intent_p3 = source.intent_param(3);

  niihdr.sclSlope = 1.0;
  niihdr.sclInter = 0.0;

  niihdr.cal_min = source.getDisplayMinimum();
  niihdr.cal_max = source.getDisplayMaximum();
  niihdr.auxillaryFile = source.getAuxFile();

  niihdr.units= NIFTI_UNITS_SEC | NIFTI_UNITS_MM; //NIFTI unit setting is formed by bitwise addition of defined values, in this case 10
  niihdr.sliceOrdering=0;

  return 0;
}

template <class T> 
int save_basic_volume(const volume<T>& source, const string& filename, 
		      int filetype, bool save_orig=false);

template <class T>
  int save_volume(const volume<T>& source, const string& filename, const int filetype)
{
  return save_basic_volume(source,filename,filetype,false);
}


template <class T>
int save_volume_dtype(const volume<T>& source, const string& filename,short datatype,const int filetype=-1)
{
  datatype=closestTemplatedType(datatype);
  if (dtype(source) == datatype) 
    return save_volume(source,filename,filetype);
  switch(datatype) {
    case DT_SIGNED_SHORT: {
      volume<short> svol;
      copyconvert(source,svol);
      return save_volume(svol,filename,filetype);
      break;
    }
    case DT_UNSIGNED_CHAR: {
      volume<char> svol;
      copyconvert(source,svol);
      return save_volume(svol,filename,filetype);
      break;
    }
    case DT_SIGNED_INT:{
      volume<int> svol;
      copyconvert(source,svol);
      return save_volume(svol,filename,filetype);
      break;
    }
    case DT_FLOAT: {
      volume<float> svol;
      copyconvert(source,svol);
      return save_volume(svol,filename,filetype);
      break;
    }
    case DT_DOUBLE: {
      volume<double> svol;
      copyconvert(source,svol);
      return save_volume(svol,filename,filetype);
      break;
    }
    default:
      ostringstream errmsg;
      errmsg << "NEWIMAGE::save_volume_dtype: DT " << datatype <<  " not supported";
      perror(errmsg.str().c_str());
  }
  return -1;  // should never get here
}
  
// functions to save without doing any swapping (i.e. just as passed in)

template <class T>
int save_orig_volume(const volume<T>& source, const string& filename, const int filetype = -1)
{
  return save_basic_volume(source,filename,filetype,true);
}


////////////////////////////////////////////////////////////////////////
///// Boring overloads to enable different names (load and write) //////
////////////////////////////////////////////////////////////////////////

// load
template <class T>
int load_volume(volume<T>& target, const string& filename)
{ return read_volume(target,filename); }

 // write
template <class T>
int write_volume(const volume<T>& source, const string& filename)
{ return save_volume(source,filename); }

// Basic I/O functions 
// read original storage order - do not swap to radiological
template <class T>
int read_orig_volume(volume<T>& target, const string& filename)
{
  short dtype;
  read_volumeROI(target,filename,dtype,0,0,0,0,-1,-1,-1,-1,false);
  return 0;
}

template <class T>
int load_orig_volume(volume<T>& target, const string& filename)
{ return read_orig_volume(target,filename); }

 	 


//TODO REMOVE after 4D is dead
template <class T>
int write_volume4D(const volume<T>& source, const string& filename)
  { return save_volume(source,filename); }
template <class T>
int save_volume4D(const volume<T>& source, const string& filename)
  { return save_volume(source,filename); }
template <class T>
int load_volume4D(volume<T>& source, const string& filename)
  { return load_volume(source,filename); }
template <class T>
  int read_volume4D(volume<T>& source, const string& filename)
  { return load_volume(source,filename); } 
template <class T>
int read_orig_volume4D(volume<T>& target, const string& filename)
{
  short dtype;
  read_volumeROI(target,filename,dtype,0,0,0,0,-1,-1,-1,-1,false);
  return 0;
 }
	 
}

#endif

