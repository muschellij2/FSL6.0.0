/*  General IO functions (images and transformation files)

    Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group

    Copyright (C) 1999-2008 University of Oxford  */

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
#include "newimageio.h"
//#include "boost/filesystem.hpp"
#include <sys/stat.h>

using namespace MISCMATHS;

namespace NEWIMAGE {

////////////////////////////////////////////////////////////////////////////

class imageExtensions {
public:
  static const int size=4;
  static const string extension[size];
};

const string imageExtensions::extension[imageExtensions::size]={".nii.gz",".nii",".hdr",".hdr.gz"};
  
bool FslIsSingleFileType(int filetype)
{
  if ( filetype % 100 >= 10 )
    return false;
  return true;
}


bool FslIsCompressedFileType(int filetype)
{
  if ( filetype >=100 ) return true;
  return false;
}

int FslNiftiVersionFileType(int filetype) //Note that this relies on the definitions in Newimage.h, if the numbering scheme
{                                         //changes this may become invalid
  return filetype %100 %10 %3;
}

int FslFiletypeFromHeader(const NiftiHeader& header) 
{
  int filetype(0);
  filetype+=header.niftiVersion();
  if ( !header.singleFile() )
    filetype+=10;
  return filetype;
}

int fslFileType(string filename)
{
  filename=return_validimagefilename(filename);
  NiftiIO reader;
  NiftiHeader header = reader.loadHeader(filename);
  int filetype(FslFiletypeFromHeader(header));
  if ( filename.substr(filename.size()-3,3) == ".gz" )
    filetype+=100;
  return filetype;
}

bool FslFileExists(const string& filename) {
  try {
    return_validimagefilename(filename,true);
  } catch(...) {
    return false;
  }
  return true;
}

bool FslImageExists(const string& filename) { 
  return FslFileExists(filename); 
}

int FslGetEnvOutputType(void)
{
  /* return type is one of FSL_TYPE_* or -1 to indicate error */
  char *otype;
  otype = getenv("FSLOUTPUTTYPE");
  if (otype == NULL) {
    fprintf(stderr,"ERROR:: Environment variable FSLOUTPUTTYPE is not set!\n");
    fprintf(stderr,"Please make sure that the appropriate configuration file is sourced by your shell (e.g. by putting it in .profile).\n");
    fprintf(stderr,"e.g. bash or sh users add the line \". ${FSLDIR}/etc/fslconf/fsl.sh\"\n");
    fprintf(stderr,"e.g. tcsh or csh users add the line \"source ${FSLDIR}/etc/fslconf/fsl.csh\"\n");
    exit(EXIT_FAILURE);
  }
  if (strcmp(otype,"NIFTI")==0) { return FSL_TYPE_NIFTI; }
  if (strcmp(otype,"NIFTI2")==0) { return FSL_TYPE_NIFTI2; }
  if (strcmp(otype,"NIFTI_GZ")==0) { return FSL_TYPE_NIFTI_GZ; }
  if (strcmp(otype,"NIFTI2_GZ")==0) { return FSL_TYPE_NIFTI2_GZ; }
  if (strcmp(otype,"NIFTI_PAIR")==0) { return FSL_TYPE_NIFTI_PAIR; }
  if (strcmp(otype,"NIFTI2_PAIR")==0) { return FSL_TYPE_NIFTI2_PAIR; }
  if (strcmp(otype,"NIFTI_PAIR_GZ")==0) { return FSL_TYPE_NIFTI_PAIR_GZ; }
  if (strcmp(otype,"NIFTI2_PAIR_GZ")==0) { return FSL_TYPE_NIFTI2_PAIR_GZ; }
  fprintf(stderr,"ERROR:: Unrecognised value (%s) of environment variable FSLOUTPUTTYPE\n",otype);
  fprintf(stderr,"Legal values are: NIFTI, NIFTI_PAIR, NIFTI_GZ, NIFTI_PAIR_GZ\n");
  exit(EXIT_FAILURE);
  return -1;
}

string outputExtension(const int filetype)
{
  if ( filetype == FSL_TYPE_NIFTI || filetype == FSL_TYPE_NIFTI2 )
    return ".nii";
  if ( filetype == FSL_TYPE_NIFTI_GZ || filetype == FSL_TYPE_NIFTI2_GZ )
    return ".nii.gz";
  if ( filetype == FSL_TYPE_ANALYZE || FSL_TYPE_NIFTI_PAIR || filetype == FSL_TYPE_NIFTI2_PAIR )
    return ".hdr";
  if ( filetype == FSL_TYPE_NIFTI_PAIR_GZ || filetype == FSL_TYPE_NIFTI2_PAIR_GZ )
    return ".hdr.gz";
  return "";
}


int NiftiGetLeftRightOrder(const NiftiHeader& niihdr)
{
  int order=FSL_RADIOLOGICAL, sform_code, qform_code;
  mat44 sform44, qform44;
  sform_code = niihdr.sformCode;
  qform_code = niihdr.qformCode;
  sform44 = niihdr.getSForm();
  qform44 = niihdr.getQForm();
  // Determines if the image is stored in neurological or radiological convention
  float dets=-1.0, detq=-1.0, det=-1.0;
  mat33 sform33, qform33;
  if (qform_code!=NIFTI_XFORM_UNKNOWN) { 
    qform33 = mat44_to_mat33(qform44);
    detq = nifti_mat33_determ(qform33);
    det = detq;
  }
  if (sform_code!=NIFTI_XFORM_UNKNOWN) { 
    sform33 = mat44_to_mat33(sform44);
    dets = nifti_mat33_determ(sform33);
    det = dets;
  }
  
  if (det<0.0) order=FSL_RADIOLOGICAL;
  else order=FSL_NEUROLOGICAL;
  // check for inconsistency if both are set
  if ( (sform_code!=NIFTI_XFORM_UNKNOWN) && 
       (qform_code!=NIFTI_XFORM_UNKNOWN) ) { 
    if (dets * detq < 0.0) order=FSL_INCONSISTENT;
    if (fabs(dets * detq)<1e-12)  order=FSL_ZERODET;
  }
  if (fabs(det)<1e-12) order=FSL_ZERODET;
  return order;
}


// VOLUME I/O
template <class T>
void set_volume_properties(const NiftiHeader& niihdr, volume<T>& target)
{
  target.setdims(niihdr.pixdim[1],niihdr.pixdim[2],niihdr.pixdim[3],niihdr.pixdim[4],
		 niihdr.pixdim[5],niihdr.pixdim[6],niihdr.pixdim[7]);
  int sform_code, qform_code;
  mat44 smat, qmat;
  qmat = niihdr.getQForm();
  qform_code = niihdr.qformCode;
  smat = niihdr.getSForm();
  sform_code = niihdr.sformCode;
  Matrix snewmat(4,4), qnewmat(4,4);
  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      snewmat(i,j) = smat.m[i-1][j-1];
      qnewmat(i,j) = qmat.m[i-1][j-1];
    }
  }
  target.set_sform(sform_code,snewmat);
  target.set_qform(qform_code,qnewmat);
  target.RadiologicalFile = (NiftiGetLeftRightOrder(niihdr)==FSL_RADIOLOGICAL);

  target.set_intent(niihdr.intentCode,niihdr.intent_p1,niihdr.intent_p2,niihdr.intent_p3);
  target.setDisplayMinimum(niihdr.cal_min);
  target.setDisplayMaximum(niihdr.cal_max);
  target.setAuxFile(niihdr.auxillaryFile);
}

template void set_volume_properties(const NiftiHeader& niihdr, volume<char>& target);
template void set_volume_properties(const NiftiHeader& niihdr, volume<short>& target);
template void set_volume_properties(const NiftiHeader& niihdr, volume<int>& target);
template void set_volume_properties(const NiftiHeader& niihdr, volume<float>& target);
template void set_volume_properties(const NiftiHeader& niihdr, volume<double>& target);


int read_volume_size(const string& filename, 
		     int64_t& sx, int64_t& sy, int64_t& sz, int64_t& st, int64_t& s5, int64_t& s6, int64_t& s7)
{
  // read in sizes only
  NiftiIO Reader;
  NiftiHeader niihdr = Reader.loadHeader(filename);
  sx=niihdr.dim[1];
  sy=niihdr.dim[2];
  sz=niihdr.dim[3];
  st=niihdr.dim[4];
  s5=niihdr.dim[5];
  s6=niihdr.dim[6];
  s7=niihdr.dim[7];

  return 0;
}

template <class T>
int read_volumeROI(volume<T>& target, const string& filename, 
		     short& dtype,
		     int64_t x0, int64_t y0, int64_t z0, int64_t t0, int64_t d50, int64_t d60, int64_t d70,
		     int64_t x1, int64_t y1, int64_t z1, int64_t t1, int64_t d51, int64_t d61, int64_t d71,
		     const bool swap2radiological, const bool readAs4D)
{
  // to get the whole volume use x0=y0=z0=t0=0 and x1=y1=z1=t1=-1
  // NB: coordinates are in "radiological" convention when swapping (i.e.
  ///    *not* the same as nifti/fslview), or untouched otherwise
  target.destroy();

  NiftiIO Reader;
  NiftiHeader niihdr;
  char *buffer;
  size_t newBufferElements;
  try {
    niihdr = Reader.loadImageROI(return_validimagefilename(filename),buffer,target.extensions,newBufferElements,x0,x1,y0,y1,z0,z1,t0,t1,d50,d51,d60,d61,d70,d71);
  } catch ( exception& e ) { imthrow("Failed to read volume "+filename+"\nError : "+e.what(),22); }
  // sanity check stuff (well, forcing sanity really)
  for (int n=1; n<=7; n++) {
    if (niihdr.dim[n]<1) niihdr.dim[n]=1;  // make it robust to dim[n]=0
  }
  // allocate and fill buffer with required data
  T* tbuffer;
  //replicated code from NewNifti as a temp fix
   if ( x0 == -1 )
     x0=0;
   if ( y0 == -1 )
     y0=0;
   if ( z0 == -1 )
     z0=0;
   if ( t0 == -1 )
     t0=0;
   if ( d50 == -1 )
     d50=0;
   if ( d60 == -1 )
     d60=0;
   if ( d70 == -1 )
     d70=0;
   if ( x1 == -1 )
     x1=niihdr.dim[1]-1;
   if ( y1 == -1 )
     y1=niihdr.dim[2]-1;
   if ( z1 == -1 )
     z1=niihdr.dim[3]-1;
   if ( t1 == -1 )
     t1=niihdr.dim[4]-1;
   if ( d51 == -1 )
     d51=niihdr.dim[5]-1;
   if ( d61 == -1 )
     d61=niihdr.dim[6]-1;
   if ( d71 == -1 )
     d71=niihdr.dim[7]-1;

   //cerr << x1 << y1 << z1 << t1 << d51 << d61 << d71 << endl;

   size_t actualElements((++x1-x0)*(++y1-y0)*(++z1-z0)*(++t1-t0)*(++d51-d50)*(++d61-d60)*(++d71-d70));
   //cerr << "foo" <<  actualElements << endl;
  ConvertAndScaleNewNiftiBuffer(buffer,tbuffer,niihdr,actualElements);  // buffer will get deleted inside (unless T=char)
  if (tbuffer==NULL)
    cout << "help" << endl;
  target.initialize(x1-x0,y1-y0,z1-z0,t1-t0,d51-d50,d61-d60,d71-d70,tbuffer,true);
  // copy info from file
  set_volume_properties(niihdr,target);
  // return value gives info about file datatype
  dtype = niihdr.datatype;
  // swap to radiological if necessary
  if (swap2radiological && !target.RadiologicalFile) target.makeradiological();
  //TODO if readAs4D use the 5Dto4D method that we _will_ write
  return 0;
}

template int read_volumeROI(volume<char>& target, const string& filename, 
			      short& dtype, 
			      int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
			      int64_t x1, int64_t y1, int64_t z1, int64_t t1,
			      const bool swap2radiological, const bool readAs4D);
template int read_volumeROI(volume<short>& target, const string& filename, 
			      short& dtype, 
			      int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
			      int64_t x1, int64_t y1, int64_t z1, int64_t t1,
			      const bool swap2radiological, const bool readAs4D);
template int read_volumeROI(volume<int>& target, const string& filename, 
			      short& dtype, 
			      int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
			      int64_t x1, int64_t y1, int64_t z1, int64_t t1,
			      const bool swap2radiological, const bool readAs4D);
template int read_volumeROI(volume<float>& target, const string& filename, 
			      short& dtype, 
			      int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
			      int64_t x1, int64_t y1, int64_t z1, int64_t t1,
			      const bool swap2radiological, const bool readAs4D);
template int read_volumeROI(volume<double>& target, const string& filename, 
			      short& dtype, 
			      int64_t x0, int64_t y0, int64_t z0, int64_t t0, 
			      int64_t x1, int64_t y1, int64_t z1, int64_t t1,
			      const bool swap2radiological, const bool readAs4D);

template <class V>
int save_unswapped_vol(const V& source, const string& filename, int filetype,int bitsPerVoxel)
{
  NiftiIO niiwrite;
  NiftiHeader niihdr;
  set_fsl_hdr(source,niihdr);
  if ( filetype<0 ) 
    filetype=FslGetEnvOutputType(); 
  niihdr.description=BUILDSTRING;
  niihdr.setNiftiVersion(FslNiftiVersionFileType(filetype),FslIsSingleFileType(filetype));
  niihdr.bitsPerVoxel=bitsPerVoxel;
  niiwrite.saveImage(make_basename(filename)+outputExtension(filetype), (const char *)source.fbegin(), source.extensions, niihdr, FslIsCompressedFileType(filetype));
  return 0;
}

template <class T>
int save_basic_volume(const volume<T>& source, const string& filename,
			int filetype, bool noSwapping)
{
  if (source.tsize()<1) return -1;
  bool currently_rad = source.left_right_order()==FSL_RADIOLOGICAL;
  bool savingAnalyze(FslNiftiVersionFileType(filetype)==0);
  if (!noSwapping && ( !source.RadiologicalFile || savingAnalyze ) && currently_rad)  const_cast< volume <T>& > (source).makeneurological();
  save_unswapped_vol(source,filename,filetype,sizeof(T)*8);
  if (!noSwapping && ( !source.RadiologicalFile || savingAnalyze ) && currently_rad)  const_cast< volume <T>& > (source).makeradiological();
  return 0;
}

template int save_basic_volume(const volume<char>& source, const string& filename,
				 int filetype, bool save_orig);
template int save_basic_volume(const volume<short>& source, const string& filename,
				 int filetype, bool save_orig);
template int save_basic_volume(const volume<int>& source, const string& filename,
				 int filetype, bool save_orig);
template int save_basic_volume(const volume<float>& source, const string& filename,
				 int filetype, bool save_orig);
template int save_basic_volume(const volume<double>& source, const string& filename,
				 int filetype, bool save_orig);

mat44 newmat2mat44(const Matrix& nmat)
{
  mat44 ret;
  for (int i=1; i<=4; i++) {
    for (int j=1; j<=4; j++) {
      ret.m[i-1][j-1] = nmat(i,j);
    }
  }
  return ret;
}


string make_basename(string& filename) //note passed as reference for compat with old API - use as a return type is preferred
{
  for ( int i=0; i < imageExtensions::size; ++i )
    if ( filename.rfind(imageExtensions::extension[i]) != string::npos ) 
      filename.erase(filename.rfind(imageExtensions::extension[i]),imageExtensions::extension[i].size());
  return(filename);
}

string make_basename(const string& inputName) 
{
  string filename(inputName);
  return(make_basename(filename));
}

string appendFSLfilename(const string inputName, const string addendum) {
  return make_basename(inputName)+addendum;
}
 

bool valid_imagefilename(const string& filename) 
{
  //return boost::filesystem::exists(filename);
  struct stat buf;
  return (stat(filename.c_str(), &buf) == 0);
}

  string return_validimagefilename(const string& filename, const bool quiet, const bool strict)
{
  string bname(make_basename(filename)), validname="";
  int validcount=0;
  if (bname != filename) {
    if (valid_imagefilename(filename))
      return filename; //If the original filename had a full image extension and exists, return it.
    else if ( strict )
      imthrow(filename+" is not an image.",62,quiet); //Strict doesn't search for next-closest match
  } 
  //Look at possible extensions ( if filename is full, this will effectively repeat the test above )
  for ( int i=0; i < imageExtensions::size && validcount < 2; ++i )
    if (valid_imagefilename(bname+imageExtensions::extension[i])) { validname=bname+imageExtensions::extension[i]; validcount++; }
  if (validcount>1)
    imthrow("Multiple possible filenames detected for basename: "+bname,61,quiet);
  if (validcount==0)
    imthrow("No image files match: "+bname,63,quiet);
  return validname;
}

int find_pathname(string& filename)
{
  if (filename.size() < 1) return -1;
  string pathname = filename;
  int fsize = pathname.length(), indx;

  // working backwards, find '/' and remove everything after it

  indx = fsize-1;
  while ((pathname[indx] != '/') && (indx != 0))
    indx--;
  
  if (indx<fsize-1)
    pathname.erase(indx+1);
  
  filename = pathname;
  return 0;
}

short closestTemplatedType(const short inputType)
{
  switch (inputType) {
  case DT_UNSIGNED_CHAR:
  case DT_INT8:
    return DT_UNSIGNED_CHAR;
  case DT_SIGNED_SHORT:
    return DT_SIGNED_SHORT;
  case DT_SIGNED_INT:
  case DT_UINT16:
    return DT_SIGNED_INT;
  case DT_FLOAT:
  case DT_UINT32:
  case DT_INT64:
  case DT_UINT64:
    return DT_FLOAT;
  case DT_DOUBLE:
  case DT_FLOAT128:
    return DT_DOUBLE;
  case DT_COMPLEX:
    cerr << "COMPLEX not supported as an independent type" << endl;
    return -1;
  default:
    cerr << "Datatype " << inputType << " is NOT supported - please check your image" << endl;
    return -1;
  }
}

short dtype(const char* T)   { return DT_UNSIGNED_CHAR; }
short dtype(const short* T)  { return DT_SIGNED_SHORT; }
short dtype(const int* T)    { return DT_SIGNED_INT; }
short dtype(const float* T)  { return DT_FLOAT; }
short dtype(const double* T) { return DT_DOUBLE; }

short dtype(const volume<char>& vol)   { return DT_UNSIGNED_CHAR; }
short dtype(const volume<short>& vol)  { return DT_SIGNED_SHORT; }
short dtype(const volume<int>& vol)    { return DT_SIGNED_INT; }
short dtype(const volume<float>& vol)  { return DT_FLOAT; }
short dtype(const volume<double>& vol) { return DT_DOUBLE; }

short dtype(const string& filename) 
{
  if ( filename.size()<1 ) return -1;
  NiftiIO niireader;
  NiftiHeader niihdr = niireader.loadHeader(return_validimagefilename(filename));
  return niihdr.datatype;
}


template <class T>
void ConvertAndScaleNewNiftiBuffer(char* buffer, T*& tbuffer, const NiftiHeader& niihdr, const size_t & imagesize)
{
  short type = niihdr.datatype;
  float slope = niihdr.sclSlope, intercept = niihdr.sclInter;
  bool doscaling = false;
  if (fabs(slope)<1e-30) {
    slope = 1.0;
    intercept = 0.0;
  }
  if ( (fabs(slope - 1.0)>1e-30) || (fabs(intercept)>1e-30) ) {
    doscaling=true;
  }
  // create buffer pointer of the desired type and allocate if necessary (scaling or not)
  if (dtype(tbuffer) != type) { tbuffer = new T[imagesize]; }
  else { tbuffer=(T*) buffer; }

  if ( (dtype(tbuffer) != type) || (doscaling) ) {
    switch(type)
      {
      case DT_SIGNED_SHORT:
	{
	  if (doscaling) convertbuffer((short *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((short *) buffer,tbuffer,imagesize);
	}
	break;
      case DT_UNSIGNED_CHAR:
	{
	  if (doscaling) convertbuffer((unsigned char *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((unsigned char *) buffer,tbuffer,imagesize);
	}
	break;
      case DT_SIGNED_INT:
	{
	  if (doscaling) convertbuffer((int *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((int *) buffer,tbuffer,imagesize);
	}
	break;
      case DT_FLOAT:
	{
	  if (doscaling) convertbuffer((float *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((float *) buffer,tbuffer,imagesize);
	}
	break;
      case DT_DOUBLE:
	{
	  if (doscaling) convertbuffer((double *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((double *) buffer,tbuffer,imagesize);
	}
	break;
	/*------------------- new codes for NIFTI ---*/
      case DT_INT8:
	{
	  if (doscaling) convertbuffer((signed char *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((signed char *) buffer,tbuffer,imagesize);
	}
	break;
      case DT_UINT16:
	{
	  if (doscaling) convertbuffer((unsigned short *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((unsigned short *) buffer,tbuffer,imagesize);
	}
	break;
      case DT_UINT32:
	{
	  if (doscaling) convertbuffer((unsigned int *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((unsigned int *) buffer,tbuffer,imagesize);
	}
	break;
      case DT_INT64:
	{
	  if (doscaling) convertbuffer((long signed int *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((long signed int *) buffer,tbuffer,imagesize);
	}
	break;
      case DT_UINT64:
	{
	  if (doscaling) convertbuffer((long unsigned int *) buffer,tbuffer,imagesize,slope,intercept);
	  else convertbuffer((long unsigned int *) buffer,tbuffer,imagesize);
	}
	break;
      default:
	  /* includes: DT_BINARY, DT_RGB, DT_ALL, DT_FLOAT128, DT_COMPLEX's */
	delete [] tbuffer;
	delete [] buffer;
	imthrow("Fslread: DT " + num2str(type) + " not supported",8);
      }
    // delete old buffer *ONLY* if a new buffer has been allocated
    if (dtype(tbuffer) != type)  delete[] buffer;
  }
}

//////////////////////////////////////////////////////////////////////////

// COMPLEX IMAGE I/O
int read_complexvolume(volume<float>& realvols, volume<float>& imagvols,
			 const string& filename, bool read_img_data)
{
  if ( filename.size()<1 ) return -1;
  NiftiIO Reader;
  NiftiHeader niihdr;
  vector<NiftiExtension> extensions; 

  try {
    niihdr = Reader.loadHeader(return_validimagefilename(filename));
  } catch ( exception& e ) { imthrow("Failed to read volume "+filename+"\nError : "+e.what(),22); }
  bool isComplex( niihdr.datatype == DT_COMPLEX );

  int64_t sx(niihdr.dim[1]),sy(niihdr.dim[2]),sz(niihdr.dim[3]),st(std::max(niihdr.dim[4],(int64_t)1));
  size_t volsize=sx*sy*sz*st;

  float *rbuffer(NULL), *ibuffer(NULL); 
  char *buffer;
  if ( read_img_data) { 
    rbuffer=new float[volsize];
    if (rbuffer==0) 
      imthrow("Out of memory",99);
    if ( isComplex ) {
      ibuffer=new float[volsize];
      if (ibuffer==0) 
	imthrow("Out of memory",99); 
    }
    size_t elements;
    niihdr = Reader.loadImageROI(return_validimagefilename(filename),buffer,extensions,elements);
    for ( size_t voxel=0;voxel<volsize && !isComplex;voxel++) //Only copy real buffer for non-complex
      rbuffer[voxel]=((float*)buffer)[voxel];
    for ( size_t voxel=0;voxel<volsize && isComplex;voxel++) {
      rbuffer[voxel]=((float *)buffer)[2*voxel];
    ibuffer[voxel]=((float *)buffer)[2*voxel+1];
    }
    realvols.reinitialize(sx,sy,sz,st,rbuffer,true);
    imagvols.reinitialize(sx,sy,sz,st,ibuffer,true);
  }

  realvols.setdims(niihdr.pixdim[1],niihdr.pixdim[2],niihdr.pixdim[3],niihdr.pixdim[4]);
  imagvols.setdims(niihdr.pixdim[1],niihdr.pixdim[2],niihdr.pixdim[3],niihdr.pixdim[4]);
  // swap to Radiological when necessary
  if ( NiftiGetLeftRightOrder(niihdr) != FSL_RADIOLOGICAL ) {
    realvols.RadiologicalFile = false;
    realvols.makeradiological();
    imagvols.RadiologicalFile = false;
    imagvols.makeradiological();
  } else {
    realvols.RadiologicalFile = true;
    imagvols.RadiologicalFile = true;
  }
  delete [] buffer;
  return 0;
}

int read_complexvolume(complexvolume& vol, const string& filename)
{ return read_complexvolume(vol.re(),vol.im(),filename,true); }

int save_complexvolume(const volume<float>& realvols, const volume<float>& imagvols, const string& filename)
{
  NiftiIO niiwrite;
  NiftiHeader niihdr;
  if (realvols.tsize()<=0) return -1;
  // convert back to Neurological if necessary
  if (!realvols.RadiologicalFile) { const_cast< volume <float>& > (realvols).makeneurological(); }
  if (!imagvols.RadiologicalFile) { const_cast< volume <float>& > (imagvols).makeneurological(); }
  set_fsl_hdr(realvols,niihdr);
  niihdr.datatype=DT_COMPLEX;
  int filetype=FslGetEnvOutputType(); 
  niihdr.setNiftiVersion(FslNiftiVersionFileType(filetype),FslIsSingleFileType(filetype));
  niihdr.bitsPerVoxel=sizeof(float)*8;
  float *buffer=new float[ 2*realvols.totalElements() ];
  volume<float>::fast_const_iterator rit(realvols.fbegin()), iit(imagvols.fbegin());
  for ( size_t voxel=0;voxel<(size_t)realvols.totalElements();++voxel) {
    buffer[2*voxel]=*rit++;
    buffer[2*voxel+1]=*iit++;
  }
  niihdr.description=BUILDSTRING;
  niiwrite.saveImage(make_basename(filename)+outputExtension(filetype), (const char *)buffer, realvols.extensions, niihdr, FslIsCompressedFileType(filetype)); 
  // restore to original ?
  if (!realvols.RadiologicalFile) { const_cast< volume <float>& > (realvols).makeradiological(); }
  if (!imagvols.RadiologicalFile) { const_cast< volume <float>& > (imagvols).makeradiological(); }
  return 0;
}

int save_complexvolume(const complexvolume& vol, const string& filename)
{ return save_complexvolume(vol.re(),vol.im(),filename); }

//////////////////////////////////////////////////////////////////////////
int load_complexvolume(complexvolume& vol, const string& filename)
  { return read_complexvolume(vol,filename); }
int load_complexvolume(volume<float>& realvol, volume<float>& imagvol,const string& filename)
  { return read_complexvolume(realvol,imagvol,filename); }
int write_complexvolume(const complexvolume& vol, const string& filename)
  { return save_complexvolume(vol,filename); }
int write_complexvolume(const volume<float>& realvol, const volume<float>& imagvol, const string& filename)
  { return save_complexvolume(realvol,imagvol,filename); }
}








