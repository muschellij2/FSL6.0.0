/*
  Matthew Webster (WIN@FMRIB)
  Copyright (C) 2018 University of Oxford  */
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

#include "NewNifti.h"

#define LSB_FIRST 1
#define MSB_FIRST 2

int systemByteOrder(void)   /* determine CPU's byte order */
{
union {
   unsigned char c[2] ;
   short         s    ;
 } testUnion ;

 testUnion.c[0] = 1 ;
 testUnion.c[1] = 0 ;

 return (testUnion.s == 1) ? LSB_FIRST : MSB_FIRST;
}

void NiftiHeader::initialise(void) {
  dim.resize(8,1);
  pixdim.resize(8,0);
  resetNonNiftiFields();
}

void NiftiHeader::sanitise(void) {
  for(char axis(1); axis<=dim[0]; axis++)
    pixdim[axis] = pixdim[axis] > 0 ? pixdim[axis] : 1;
}

NiftiHeader::NiftiHeader()
{
  initialise();
}

void NiftiHeader::resetNonNiftiFields() {
  pixdim[0]=-1;
  sclSlope=sclInter=sliceStart=sliceEnd=sliceDuration=sliceCode=units=toffset=0;
  qformCode=sformCode=qB=qC=qD=qX=qY=qZ=0;
  sX.resize(4,0);
  sY.resize(4,0);
  sZ.resize(4,0);
  intentName=string(16,'\0');
  magic=string(4,'\0');
}


int NiftiHeader::bpvOfDatatype(void) {
  switch( datatype ) {
    case DT_UNSIGNED_CHAR:
    case DT_INT8:
      return 8;
    case DT_SIGNED_SHORT:
    case DT_UINT16:
      return 16;
    case DT_SIGNED_INT:
    case DT_FLOAT:
    case DT_UINT32:
      return 32;
    case DT_INT64:
    case DT_UINT64:
    case DT_DOUBLE:
    case DT_COMPLEX:
      return 64;
    case DT_FLOAT128:
      return 128;
  }
  return 0;
}

string NiftiHeader::niftiOrientationString( const int orientation ) const
{
  switch( orientation ) {
    case NIFTI_L2R: return "Left-to-Right" ;
    case NIFTI_R2L: return "Right-to-Left" ;
    case NIFTI_P2A: return "Posterior-to-Anterior" ;
    case NIFTI_A2P: return "Anterior-to-Posterior" ;
    case NIFTI_I2S: return "Inferior-to-Superior" ;
    case NIFTI_S2I: return "Superior-to-Inferior" ;
  }
  return "Unknown" ;
}

string NiftiHeader::niftiSliceString() const {
   switch( sliceCode ){
     case NIFTI_SLICE_SEQ_INC:  return "sequential_increasing";
     case NIFTI_SLICE_SEQ_DEC:  return "sequential_decreasing";
     case NIFTI_SLICE_ALT_INC:  return "alternating_increasing";
     case NIFTI_SLICE_ALT_DEC:  return "alternating_decreasing";
     case NIFTI_SLICE_ALT_INC2: return "alternating_increasing_2";
     case NIFTI_SLICE_ALT_DEC2: return "alternating_decreasing_2";
   }
   return "Unknown";
}

string NiftiHeader::niftiIntentString() const {
  switch( intentCode ){
    case NIFTI_INTENT_CORREL:     return "Correlation statistic";
    case NIFTI_INTENT_TTEST:      return "T-statistic";
    case NIFTI_INTENT_FTEST:      return "F-statistic";
    case NIFTI_INTENT_ZSCORE:     return "Z-score";
    case NIFTI_INTENT_CHISQ:      return "Chi-squared distribution";
    case NIFTI_INTENT_BETA:       return "Beta distribution";
    case NIFTI_INTENT_BINOM:      return "Binomial distribution";
    case NIFTI_INTENT_GAMMA:      return "Gamma distribution";
    case NIFTI_INTENT_POISSON:    return "Poisson distribution";
    case NIFTI_INTENT_NORMAL:     return "Normal distribution";
    case NIFTI_INTENT_FTEST_NONC: return "F-statistic noncentral";
    case NIFTI_INTENT_CHISQ_NONC: return "Chi-squared noncentral";
    case NIFTI_INTENT_LOGISTIC:   return "Logistic distribution";
    case NIFTI_INTENT_LAPLACE:    return "Laplace distribution";
    case NIFTI_INTENT_UNIFORM:    return "Uniform distribition";
    case NIFTI_INTENT_TTEST_NONC: return "T-statistic noncentral";
    case NIFTI_INTENT_WEIBULL:    return "Weibull distribution";
    case NIFTI_INTENT_CHI:        return "Chi distribution";
    case NIFTI_INTENT_INVGAUSS:   return "Inverse Gaussian distribution";
    case NIFTI_INTENT_EXTVAL:     return "Extreme Value distribution";
    case NIFTI_INTENT_PVAL:       return "P-value";
    case NIFTI_INTENT_LOGPVAL:    return "Log P-value";
    case NIFTI_INTENT_LOG10PVAL:  return "Log10 P-value";
    case NIFTI_INTENT_ESTIMATE:   return "Estimate";
    case NIFTI_INTENT_LABEL:      return "Label index";
    case NIFTI_INTENT_NEURONAME:  return "NeuroNames index";
    case NIFTI_INTENT_GENMATRIX:  return "General matrix";
    case NIFTI_INTENT_SYMMATRIX:  return "Symmetric matrix";
    case NIFTI_INTENT_DISPVECT:   return "Displacement vector";
    case NIFTI_INTENT_VECTOR:     return "Vector";
    case NIFTI_INTENT_POINTSET:   return "Pointset";
    case NIFTI_INTENT_TRIANGLE:   return "Triangle";
    case NIFTI_INTENT_QUATERNION: return "Quaternion";
    case NIFTI_INTENT_DIMLESS:    return "Dimensionless number";
  }
  return "Unknown" ;
}

string NiftiHeader::niftiTransformString(const int transform) const {
  switch( transform ) {
    case NIFTI_XFORM_SCANNER_ANAT:  return "Scanner Anat";
    case NIFTI_XFORM_ALIGNED_ANAT:  return "Aligned Anat";
    case NIFTI_XFORM_TALAIRACH:     return "Talairach";
    case NIFTI_XFORM_MNI_152:       return "MNI_152";
   }
  return "Unknown" ;
}

string NiftiHeader::unitsString(const int units) const {
   switch( units ){
     case NIFTI_UNITS_METER:  return "m" ;
     case NIFTI_UNITS_MM:     return "mm" ;
     case NIFTI_UNITS_MICRON: return "um" ;
     case NIFTI_UNITS_SEC:    return "s" ;
     case NIFTI_UNITS_MSEC:   return "ms" ;
     case NIFTI_UNITS_USEC:   return "us" ;
     case NIFTI_UNITS_HZ:     return "Hz" ;
     case NIFTI_UNITS_PPM:    return "ppm" ;
     case NIFTI_UNITS_RADS:   return "rad/s" ;
   }
   return "Unknown";
}

string NiftiHeader::datatypeString() const {
   switch( datatype ){
     case DT_UNKNOWN:    return "UNKNOWN"    ;
     case DT_BINARY:     return "BINARY"     ;
     case DT_INT8:       return "INT8"       ;
     case DT_UINT8:      return "UINT8"      ;
     case DT_INT16:      return "INT16"      ;
     case DT_UINT16:     return "UINT16"     ;
     case DT_INT32:      return "INT32"      ;
     case DT_UINT32:     return "UINT32"     ;
     case DT_INT64:      return "INT64"      ;
     case DT_UINT64:     return "UINT64"     ;
     case DT_FLOAT32:    return "FLOAT32"    ;
     case DT_FLOAT64:    return "FLOAT64"    ;
     case DT_FLOAT128:   return "FLOAT128"   ;
     case DT_COMPLEX64:  return "COMPLEX64"  ;
     case DT_COMPLEX128: return "COMPLEX128" ;
     case DT_COMPLEX256: return "COMPLEX256" ;
     case DT_RGB24:      return "RGB24"      ;
   }
   return "**ILLEGAL**" ;
}

string NiftiHeader::originalOrder() const {
  if ( systemByteOrder() == LSB_FIRST && wasWrongEndian )
    return("MSB_FIRST");
  return("LSB_FIRST");
}

string NiftiHeader::fileType() const {
  if ( isAnalyze() )
    return "ANALYZE-7.5";
  string type("NIFTI-"+string(1,char(niftiVersion()+'0')) );
  if (singleFile())
    type+="+";
  return type;
}

void NiftiHeader::setNiftiVersion( const char niftiVersion, const bool isSingleFile ) {
  sizeof_hdr=348;
  if ( niftiVersion==0 )
    return;
  magic=string("ni1\0",4);
  if ( isSingleFile )
    magic[1]='+';
  magic[2]=niftiVersion+(int)'0';
  if ( niftiVersion == 2 ) {
    magic.append("\r\n\032\n",4);
    sizeof_hdr=540;
  }
}


mat44 NiftiHeader::getQForm() const {
  return ( nifti_quatern_to_mat44( qB, qC, qD, qX, qY, qZ, pixdim[1], pixdim[2], pixdim[3], leftHanded() ) );
}

void NiftiHeader::setQForm(const mat44& qForm) {
  double dx,dy,dz;
  nifti_mat44_to_quatern( qForm , qB, qC, qD, qX, qY, qZ, dx, dy, dz, pixdim[0] );
}

mat44 NiftiHeader::getSForm() const {
   mat44 R ;
   for ( int i = 0; i <= 3; i++ ) {
     R.m[0][i]=sX.at(i);
     R.m[1][i]=sY.at(i);
     R.m[2][i]=sZ.at(i);
   }
   R.m[3][0]=R.m[3][1]=R.m[3][2] = 0.0 ;
   R.m[3][3]= 1.0 ;
   return R;

}

void NiftiHeader::setSForm(const mat44& sForm) {
  sX.resize(4,0);
  sY.resize(4,0);
  sZ.resize(4,0);
  for ( int i = 0; i <= 3; i++ ) {
    sX.at(i)=sForm.m[0][i];
    sY.at(i)=sForm.m[1][i];
    sZ.at(i)=sForm.m[2][i];
   }
}

NiftiIO::NiftiIO(void)
{
  debug=false;
  file=NULL;
}

//Public Methods
//This loads in the section of data stored between the limits input ( a value of -1 will default to either 0 ( for minimum limits ) or dim[foo]-1 ( for maximum limits ).
NiftiHeader NiftiIO::loadImageROI(string filename, char*& buffer, vector<NiftiExtension>& extensions, size_t& bufferElements, int64_t xmin, int64_t xmax, int64_t ymin, int64_t ymax, int64_t zmin, int64_t zmax, int64_t tmin, int64_t tmax, int64_t d5min, int64_t d5max, int64_t d6min, int64_t d6max, int64_t d7min, int64_t d7max)
{
   openImage(filename, true);
   NiftiHeader header=readHeader();
   vector<int64_t> dims(header.dim); //set dims above dim[0] to 1 so 7D loops below work.
   for(int badDim=dims[0]+1;badDim<=7;badDim++)
     dims[badDim]=1;
   readExtensions(header, extensions );
   if ( xmin == -1 )
     xmin=0;
   if ( ymin == -1 )
     ymin=0;
   if ( zmin == -1 )
     zmin=0;
   if ( tmin == -1 )
     tmin=0;
   if ( d5min == -1 )
     d5min=0;
   if ( d6min == -1 )
     d6min=0;
   if ( d7min == -1 )
     d7min=0;
   if ( xmax == -1 )
     xmax=dims[1]-1;
   if ( ymax == -1 )
     ymax=dims[2]-1;
   if ( zmax == -1 )
     zmax=dims[3]-1;
   if ( tmax == -1 )
     tmax=dims[4]-1;
   if ( d5max == -1 )
     d5max=dims[5]-1;
   if ( d6max == -1 )
     d6max=dims[6]-1;
   if ( d7max == -1 )
     d7max=dims[7]-1;

   //cerr << xmin << " " << xmax << " "  << ymin << " " << ymax << " " << zmin << " " << zmax << " " << tmin << " " << tmax << " " << d5min << " " << d5max << " " << d6min << " " << d6max << " " << d7min << " " << d7max << endl;

   if ( xmin < 0 || xmax > ( dims[1]-1 ) || ymin < 0 || ymax > ( dims[2]-1 ) || zmin < 0 || zmax > ( dims[3]-1 ) ||  tmin < 0 || tmax > ( dims[4]-1 ) || d5min < 0 || d5max > ( dims[5]-1 ) || d6min < 0 || d6max > ( dims[6]-1 ) || d7min < 0 || d7max > ( dims[7]-1 ) )
     throw NiftiException("Error: ROI out of bounds for "+string(filename));
   if ( xmin > xmax || ymin > ymax || zmin > zmax || tmin > tmax || d5min > d5max || d6min > d6max || d7min > d7max )
         throw NiftiException("Error: Nonsensical ROI for "+string(filename));

   bufferElements=( xmax-xmin+1 ) * ( ymax-ymin+1 ) * ( zmax-zmin+1 ) * ( tmax-tmin+1 ) * ( d5max-d5min+1 ) * ( d6max-d6min+1 ) * ( d7max-d7min+1 );
   buffer = new char[bufferElements*(header.bitsPerVoxel/8)];
   char *movingBuffer(buffer);



   size_t voxelsToRead(0);
   size_t voxelsToSeek(0);
   bool wasInROI( 0 >= xmin && 0 <= xmax && 0 >= ymin && 0 <= ymax && 0 >= zmin && 0 <= zmax && 0 >= tmin && 0 <= tmax && 0 >= d5min && 0 <= d5max && 0 >= d6min && 0 <= d6max && 0 >= d7min && 0 <= d7max );
   if ( !header.singleFile() )
     openImage( filename.replace(filename.rfind(".hdr"),4,".img"), true );
   znzseek(file, header.vox_offset, SEEK_SET);
   for ( int64_t d7 = 0; d7 < dims[7]; d7++ )
     for ( int64_t d6 = 0; d6 < dims[6]; d6++ )
       for ( int64_t d5 = 0; d5 < dims[5]; d5++ )
	 for ( int64_t t = 0; t < dims[4]; t++ )
	   for ( int64_t z = 0; z < dims[3]; z++ )
	     for ( int64_t y = 0; y < dims[2]; y++ )
	       for ( int64_t x = 0; x < dims[1]; x++ ) {
		 bool currentlyInROI( x >= xmin && x <= xmax && y >= ymin && y <= ymax && z >= zmin && z <= zmax && t >= tmin && t <= tmax && d5 >= d5min && d5 <= d5max && d6 >= d6min && d6 <= d6max && d7 >= d7min && d7 <= d7max );
		 if ( currentlyInROI )
		   voxelsToRead++;
		 else
		   voxelsToSeek++;
		 if ( wasInROI != currentlyInROI ) { //We have switched so flush according to wasInROI
		   if ( wasInROI ) { //read the accumulated bytes
		     //cerr << "reading " << voxelsToRead << " voxels" << x << " " << y << " " << z << endl;
		     readRawBytes(movingBuffer, voxelsToRead*(header.bitsPerVoxel/8) );
		     movingBuffer+=voxelsToRead*(header.bitsPerVoxel/8);
		     voxelsToRead=0;
		   }
		   else {
		     //cerr << "seeking " << voxelsToSeek << " voxels" << x << " " << y << " " << z <<endl;
		     znzseek(file, voxelsToSeek*(header.bitsPerVoxel/8), SEEK_CUR );
		     voxelsToSeek=0;
		   }
		   wasInROI = currentlyInROI;

		 }
	       }
   if ( voxelsToRead ) { //We still have some bytes to read - maybe the whole file?
     //cerr << "reading final " << voxelsToRead  << " voxels" << endl;
     readRawBytes(movingBuffer, voxelsToRead*(header.bitsPerVoxel/8) );
   }
   if ( voxelsToSeek ) {
     //cerr << "skipping final " << voxelsToSeek  << " voxels" << endl;
   }
   znzclose(file);

   if ( header.wasWrongEndian )
     byteSwap( header.bitsPerVoxel/8, buffer, bufferElements );
   return header;
}


NiftiHeader NiftiIO::loadImage(string filename, char*& buffer, vector<NiftiExtension>& extensions, bool allocateBuffer)
{
   openImage(filename, true);
   NiftiHeader header=readHeader();
   readExtensions(header, extensions );
   if ( !header.singleFile() )
     openImage( filename.replace(filename.rfind(".hdr"),4,".img"), true );
   if ( allocateBuffer )
     buffer = new char[header.nElements()*(header.bitsPerVoxel/8)];
   readData( header, buffer);
   znzclose(file);
   return header;
}

NiftiHeader NiftiIO::loadHeader(const string filename)
{
   openImage(filename, true);
   NiftiHeader header=readHeader();
   znzclose(file);
   return header;
}

NiftiHeader NiftiIO::loadExtensions(const string filename, vector<NiftiExtension>& extensions)
{
   openImage(filename, true);
   NiftiHeader header=readHeader();
   readExtensions(header, extensions );
   znzclose(file);
   return header;
}

void NiftiIO::saveImage(string filename, const char* buffer, const vector<NiftiExtension>& extensions, NiftiHeader header, const bool useCompression, const bool sanitise)
{
  openImage(filename, false, useCompression);
  long sizeOfExtensions(4);
  if ( sanitise )
    header.sanitise();
  header.vox_offset=header.sizeof_hdr;
  for ( unsigned int i(0); i < extensions.size(); i++ )
    sizeOfExtensions+=extensions[0].extensionSize();
  header.vox_offset+=sizeOfExtensions;
  if ( !header.singleFile() )
    header.vox_offset=0;
  writeHeader(header);
  if ( !header.isAnalyze() )
    writeExtensions(header,extensions);
  if ( !header.singleFile() ) {
    znzclose(file);
    openImage(filename.replace(filename.rfind(".hdr"),4,".img"), false, useCompression);
  }
  writeData(header,buffer);
  znzclose(file);
}
//End of public
void NiftiIO::writeHeader(const NiftiHeader& header)
{
  if ( header.isAnalyze() )
    writeAnalyzeHeader(header);
  else if ( header.niftiVersion() == 1 )
     writeNifti1Header(header);
  else if ( header.niftiVersion() == 2 )
    writeNifti2Header(header);
}

void NiftiIO::writeAnalyzeHeader(const NiftiHeader& header)
{
nifti_1_header rawNiftiHeader;
  retrieveCommonImageHeader(rawNiftiHeader, header );
  rawNiftiHeader.pixdim[0]=0;
  short originator[5] = { 0,0,0,0,0 };
  if (header.sformCode != NIFTI_XFORM_UNKNOWN) {
    mat44 inverse=nifti_mat44_inverse(header.getSForm());
    originator[0]=(short) inverse.m[0][3] + 1;
    originator[1]=(short) inverse.m[1][3] + 1;
    originator[2]=(short) inverse.m[2][3] + 1;
  } else if (header.qformCode != NIFTI_XFORM_UNKNOWN) {
    mat44 inverse=nifti_mat44_inverse(header.getSForm());
    originator[0]=(short) inverse.m[0][3] + 1;
    originator[1]=(short) inverse.m[1][3] + 1;
    originator[2]=(short) inverse.m[2][3] + 1;
  }
  memcpy(rawNiftiHeader.aux_file+25,originator,10);
  writeRawBytes( &rawNiftiHeader, (size_t)sizeof(rawNiftiHeader) );
}


void NiftiIO::writeNifti1Header(const NiftiHeader& header)
{
nifti_1_header rawNiftiHeader;
   retrieveCommonNiftiHeader(rawNiftiHeader, header );
   rawNiftiHeader.regular='r'; //set non-zero legacy field
   writeRawBytes( &rawNiftiHeader, (size_t)sizeof(rawNiftiHeader) );
}

void NiftiIO::writeNifti2Header(const NiftiHeader& header)
{
nifti_2_header rawNiftiHeader;
   retrieveCommonNiftiHeader(rawNiftiHeader, header );
   writeRawBytes( &rawNiftiHeader, (size_t)sizeof(rawNiftiHeader) );
}

template<class T>
void NiftiIO::retrieveCommonImageHeader( T& rawNiftiHeader, const NiftiHeader& header )
{ //These are fields shared by NIFTI and ANALYZE
  memset(&rawNiftiHeader,0,sizeof(rawNiftiHeader)); //Zero header for safety
  rawNiftiHeader.sizeof_hdr = header.sizeof_hdr;
  for ( int i = 0; i <= 7; i++ )
    rawNiftiHeader.dim[i] = header.dim[i];
  rawNiftiHeader.datatype = header.datatype;
  rawNiftiHeader.bitpix = header.bitsPerVoxel;
  for ( int i = 0; i <= 7; i++ )
    rawNiftiHeader.pixdim[i] = header.pixdim[i];
  rawNiftiHeader.vox_offset = header.vox_offset;
  rawNiftiHeader.cal_max = header.cal_max;
  rawNiftiHeader.cal_min = header.cal_min;
}

template<class T>
void NiftiIO::retrieveCommonNiftiHeader( T& rawNiftiHeader, const NiftiHeader& header )
{
   retrieveCommonImageHeader(rawNiftiHeader,header);
   memcpy ( rawNiftiHeader.magic, header.magic.data(), header.magic.size() );
   rawNiftiHeader.intent_p1 = header.intent_p1;
   rawNiftiHeader.intent_p2 = header.intent_p2;
   rawNiftiHeader.intent_p3 = header.intent_p3;
   strncpy ( rawNiftiHeader.intent_name, header.intentName.c_str(), 16 );
   rawNiftiHeader.intent_code = header.intentCode;
   rawNiftiHeader.dim_info = header.sliceOrdering;
   rawNiftiHeader.slice_start = header.sliceStart;
   rawNiftiHeader.slice_end = header.sliceEnd;
   rawNiftiHeader.slice_code = header.sliceCode;
   rawNiftiHeader.slice_duration = header.sliceDuration;
   rawNiftiHeader.scl_slope = header.sclSlope;
   rawNiftiHeader.scl_inter = header.sclInter;
   for ( int i = 0; i <= 3; i++ ) {
     rawNiftiHeader.srow_x[i] = header.sX[i];
     rawNiftiHeader.srow_y[i] = header.sY[i];
     rawNiftiHeader.srow_z[i] = header.sZ[i];
   }
   rawNiftiHeader.quatern_b = header.qB;
   rawNiftiHeader.quatern_c = header.qC;
   rawNiftiHeader.quatern_d = header.qD;
   rawNiftiHeader.qoffset_x = header.qX;
   rawNiftiHeader.qoffset_y = header.qY;
   rawNiftiHeader.qoffset_z = header.qZ;
   strncpy ( rawNiftiHeader.aux_file, header.auxillaryFile.c_str(), 24 );
   strncpy ( rawNiftiHeader.descrip, header.description.c_str(), 80 );
   rawNiftiHeader.sform_code = header.sformCode;
   rawNiftiHeader.qform_code = header.qformCode;
   rawNiftiHeader.xyzt_units = header.units;
   rawNiftiHeader.toffset = header.toffset;
}

void NiftiIO::openImage(const string filename, const bool reading, const bool useCompression)
{
string options;
  if ( debug )
    cout << "filename\t" << filename << endl;
  if ( reading )
    options="rb";
  else
    options="wb";
  if ( ! znz_isnull(file) ) //close currently open file
    znzclose(file);

  if ( useCompression )
    file=znzopen( filename.c_str(), options.c_str(), 1 );
  else
    file=znzopen( filename.c_str(), options.c_str(), 0 );

  if( znz_isnull(file) ) {
    throw NiftiException("Error: cant open file "+string(filename));
  }
}

void NiftiIO::readVersion( NiftiHeader& header )
{
   //read first 4 bytes from file to determine file type and endian-ness
   readRawBytes( &header.sizeof_hdr, 4);
   if ( header.niftiVersion() != 1 && header.niftiVersion() != 2 ) {
     znzclose(file);
     throw NiftiException("Error: file does not appear to be a valid NIFTI image");
   }

}

void NiftiIO::writeExtensions(const NiftiHeader& header, const vector<NiftiExtension>& extensions )
{
  char temp[4];
  temp[0]=temp[1]=temp[2]=temp[3]=0;
  if ( extensions.size() > 0 ) //We have extensions
    temp[0]=1;
  writeRawBytes(temp, 4 );
  for ( unsigned int current=0; current < extensions.size(); current++ ) {
     if ( extensions[current].esize % 16 != 0) {
       znzclose(file);
       throw NiftiException("Error: extension size must be multiple of 16 bytes");
     }
    writeRawBytes(&extensions[current].esize,4);
    writeRawBytes(&extensions[current].ecode,4);
    writeRawBytes(extensions[current].edata.data(),extensions[current].esize-8);
  }
}

void NiftiIO::readExtensions(const NiftiHeader& header, vector<NiftiExtension>& extensions)
{
  znzseek(file, header.sizeof_hdr, SEEK_SET);
  //To determine if extensions are present, read next 4 bytes ( if present ) if byte[0] is non-zero we have extensions
  char buffer[4];
  try { readRawBytes( buffer, 4 ); }
  catch ( NiftiException ) { //Short read at this point implies NIFTI-pair with no extension
    return;
  }
  //Otherwise we now have 4 bytes - test the first
  if ( buffer[0] != 0 ) { //We have extensions - begin looping over 8 bytes ( 2 ints ) + buffer
    bool moreExtensionsToRead(true);
    int iBuffer[2];
    long bytesRead(header.sizeof_hdr+4);
    while ( moreExtensionsToRead ) {
    NiftiExtension currentExtension;
       readRawBytes( iBuffer, 8 );
       if ( header.wasWrongEndian )
	 byteSwap(1, iBuffer, 8);
       currentExtension.esize=iBuffer[0];
       currentExtension.ecode=iBuffer[1];
       currentExtension.edata.resize(currentExtension.esize-8);
       readRawBytes( currentExtension.edata.data(), currentExtension.esize-8 );
       extensions.push_back( currentExtension );
       bytesRead+=currentExtension.esize;

       if ( header.singleFile() && bytesRead < header.vox_offset ) //We must have at least one more extension
	 moreExtensionsToRead=true;
       else
	 moreExtensionsToRead=false; //for paired image need to check not at end of filesp
    }
    if ( header.singleFile() && bytesRead > header.vox_offset ) { //Have we gone past vox_offset ( for single files
     znzclose(file);
     throw NiftiException("Error: Extension read overflows start of data");
    }
  }
}

void NiftiIO::readData(const NiftiHeader& header,void* buffer)
{
  znzseek(file, header.vox_offset, SEEK_SET);
  readRawBytes(buffer, header.nElements()*(header.bitsPerVoxel/8) );
  if ( header.wasWrongEndian )
    byteSwap( header.bitsPerVoxel/8, buffer, header.nElements() );
}

void NiftiIO::writeData(const NiftiHeader& header,const void* buffer)
{
  writeRawBytes(buffer, header.nElements()*(header.bitsPerVoxel/8) );
}

NiftiHeader NiftiIO::readHeader()
{
NiftiHeader header;
   readVersion(header);
   if ( header.niftiVersion() == 1 )
     readNifti1Header(header);
   if ( header.niftiVersion() == 2 )
     readNifti2Header(header);
   return header;
}

void NiftiIO::readRawBytes(void* buffer, size_t length)
{
  size_t bytes = (size_t)znzread( buffer, 1, (size_t)length, file );
   if ( bytes != length ) {
     znzclose(file);
     throw NiftiException("Error: short read, file may be truncated");
   }
}

void NiftiIO::writeRawBytes(const void* buffer, size_t length )
{
  size_t bytes = (size_t)znzwrite( buffer, 1, (size_t)length, file );
  if ( bytes != length ) {
     znzclose(file);
     throw NiftiException("Error: short write, output file will be truncated");
  }
}

template<class T>
void NiftiIO::storeCommonNiftiHeader( const T& rawNiftiHeader, NiftiHeader& header )
{
   header.sizeof_hdr = rawNiftiHeader.sizeof_hdr;
   header.magic = string(rawNiftiHeader.magic, sizeof(rawNiftiHeader.magic) );
   header.vox_offset = rawNiftiHeader.vox_offset;
   header.datatype = rawNiftiHeader.datatype;
   header.bitsPerVoxel = rawNiftiHeader.bitpix;
   header.dim.resize(8);
   for ( int i = 0; i <= 7; i++ )
     header.dim[i] = rawNiftiHeader.dim[i];
   header.intent_p1 = rawNiftiHeader.intent_p1;
   header.intent_p2 = rawNiftiHeader.intent_p2;
   header.intent_p3 = rawNiftiHeader.intent_p3;
   header.intentName = string(rawNiftiHeader.intent_name);
   header.intentCode = rawNiftiHeader.intent_code;
   header.pixdim.resize(8);
   for ( int i = 0; i <= 7; i++ )
     header.pixdim[i] = rawNiftiHeader.pixdim[i];
   header.sliceOrdering = rawNiftiHeader.dim_info;
   header.sliceStart = rawNiftiHeader.slice_start;
   header.sliceEnd = rawNiftiHeader.slice_end;
   header.sliceCode = rawNiftiHeader.slice_code;
   header.sliceDuration = rawNiftiHeader.slice_duration;
   header.sclSlope = rawNiftiHeader.scl_slope;
   header.sclInter = rawNiftiHeader.scl_inter;

   header.sX.resize(4);
   header.sY.resize(4);
   header.sZ.resize(4);
   for ( int i = 0; i <= 3; i++ ) {
     header.sX[i] = rawNiftiHeader.srow_x[i];
     header.sY[i] = rawNiftiHeader.srow_y[i];
     header.sZ[i] = rawNiftiHeader.srow_z[i];
   }
   header.qB = rawNiftiHeader.quatern_b;
   header.qC = rawNiftiHeader.quatern_c;
   header.qD = rawNiftiHeader.quatern_d;
   header.qX = rawNiftiHeader.qoffset_x;
   header.qY = rawNiftiHeader.qoffset_y;
   header.qZ = rawNiftiHeader.qoffset_z;


   header.auxillaryFile = string(rawNiftiHeader.aux_file);
   header.description = string(rawNiftiHeader.descrip);
   header.sformCode = rawNiftiHeader.sform_code;
   header.qformCode = rawNiftiHeader.qform_code;
   header.cal_max = rawNiftiHeader.cal_max;
   header.cal_min = rawNiftiHeader.cal_min;
   header.units = rawNiftiHeader.xyzt_units;
   header.toffset = rawNiftiHeader.toffset;

}

void NiftiIO::readNifti1Header(NiftiHeader& header)
{
nifti_1_header rawNiftiHeader;
  znzrewind(file);
  readRawBytes( &rawNiftiHeader, (size_t)sizeof(rawNiftiHeader) );
  header.wasWrongEndian = NIFTI2_NEEDS_SWAP(rawNiftiHeader);
  //Store original orient+originator for later
  char analyzeOrient[11];
  memcpy(analyzeOrient,&(rawNiftiHeader.qform_code),11);
  short analyzeOriginator[5];
  memcpy(analyzeOriginator,analyzeOrient+1,10);
  if ( header.wasWrongEndian ) {
    byteSwap( rawNiftiHeader );
    byteSwapLegacy(rawNiftiHeader);
    byteSwap(sizeof(analyzeOriginator[0]),analyzeOriginator,5);
  }
  storeCommonNiftiHeader( rawNiftiHeader, header );
  if ( header.isAnalyze() ) {
    header.resetNonNiftiFields();
    header.sX[0]=header.pixdim[1];
    header.sY[1]=header.pixdim[2];
    header.sZ[2]=header.pixdim[3];
    header.sX[3]=-(analyzeOriginator[0]-1)*header.pixdim[1];
    header.sY[3]=-(analyzeOriginator[1]-1)*header.pixdim[2];
    header.sZ[3]=-(analyzeOriginator[2]-1)*header.pixdim[3];
    header.setQForm(header.getSForm());
    header.qformCode=header.sformCode=NIFTI_XFORM_ALIGNED_ANAT;
    for ( int i = 1; i <= 7; i++ ) //pixdims 1..7 must be +ve for NIFTI
      header.pixdim[i] = fabs(header.pixdim[i]);
  }
  if ( debug ) {
    if ( header.wasWrongEndian )
      cout << "  Byte-Swapped" << endl;
    reportHeader(rawNiftiHeader);
    reportLegacyHeader(rawNiftiHeader);
  }
}

void NiftiIO::readNifti2Header(NiftiHeader& header)
{
nifti_2_header rawNiftiHeader;
  znzrewind(file);
  readRawBytes( &rawNiftiHeader, (size_t)sizeof(rawNiftiHeader) );
  header.wasWrongEndian = NIFTI2_NEEDS_SWAP(rawNiftiHeader);
  if ( header.wasWrongEndian )
    byteSwap( rawNiftiHeader );
  storeCommonNiftiHeader( rawNiftiHeader, header );
  if ( debug ) {
    reportHeader(rawNiftiHeader);
  if ( header.magic.substr(4,4) != "\r\n\032\n" )
     throw NiftiException("Error: bad NIFTI2 signature");
  }
}

void NiftiIO::byteSwapLegacy(nifti_1_header& rawNiftiHeader)
{
    byteSwap(sizeof(rawNiftiHeader.extents),&rawNiftiHeader.extents);
    byteSwap(sizeof(rawNiftiHeader.session_error),&rawNiftiHeader.session_error);
    byteSwap(sizeof(rawNiftiHeader.glmax),&rawNiftiHeader.glmax);
    byteSwap(sizeof(rawNiftiHeader.glmin),&rawNiftiHeader.glmin);
}


template<class T>
void NiftiIO::byteSwap(T& rawNiftiHeader)
{
  //cerr << "Byte Swapping common fields of NIFTI1/2 Header" << endl;
  byteSwap(sizeof(rawNiftiHeader.sizeof_hdr),&rawNiftiHeader.sizeof_hdr);
  byteSwap(sizeof(rawNiftiHeader.dim[0]),rawNiftiHeader.dim,8);
  byteSwap(sizeof(rawNiftiHeader.intent_p1),&rawNiftiHeader.intent_p1);
  byteSwap(sizeof(rawNiftiHeader.intent_p2),&rawNiftiHeader.intent_p2);
  byteSwap(sizeof(rawNiftiHeader.intent_p3),&rawNiftiHeader.intent_p3);
  byteSwap(sizeof(rawNiftiHeader.intent_code),&rawNiftiHeader.intent_code);
  byteSwap(sizeof(rawNiftiHeader.datatype),&rawNiftiHeader.datatype);
  byteSwap(sizeof(rawNiftiHeader.bitpix),&rawNiftiHeader.bitpix);
  byteSwap(sizeof(rawNiftiHeader.slice_start),&rawNiftiHeader.slice_start);
  byteSwap(sizeof(rawNiftiHeader.pixdim[0]),rawNiftiHeader.pixdim,8);
  byteSwap(sizeof(rawNiftiHeader.vox_offset),&rawNiftiHeader.vox_offset);
  byteSwap(sizeof(rawNiftiHeader.scl_slope),&rawNiftiHeader.scl_slope);
  byteSwap(sizeof(rawNiftiHeader.scl_inter),&rawNiftiHeader.scl_inter);
  byteSwap(sizeof(rawNiftiHeader.slice_end),&rawNiftiHeader.slice_end);
  byteSwap(sizeof(rawNiftiHeader.slice_code),&rawNiftiHeader.slice_code); //char in Nifti1, int in Nifti2
  byteSwap(sizeof(rawNiftiHeader.xyzt_units),&rawNiftiHeader.xyzt_units); //char in Nifti1, int in Nifti2
  byteSwap(sizeof(rawNiftiHeader.cal_max),&rawNiftiHeader.cal_max);
  byteSwap(sizeof(rawNiftiHeader.cal_min),&rawNiftiHeader.cal_min);
  byteSwap(sizeof(rawNiftiHeader.slice_duration),&rawNiftiHeader.slice_duration);
  byteSwap(sizeof(rawNiftiHeader.toffset),&rawNiftiHeader.toffset);
  byteSwap(sizeof(rawNiftiHeader.sform_code),&rawNiftiHeader.sform_code);
  byteSwap(sizeof(rawNiftiHeader.qform_code),&rawNiftiHeader.qform_code);
  byteSwap(sizeof(rawNiftiHeader.quatern_b),&rawNiftiHeader.quatern_b);
  byteSwap(sizeof(rawNiftiHeader.quatern_c),&rawNiftiHeader.quatern_c);
  byteSwap(sizeof(rawNiftiHeader.quatern_d),&rawNiftiHeader.quatern_d);
  byteSwap(sizeof(rawNiftiHeader.qoffset_x),&rawNiftiHeader.qoffset_x);
  byteSwap(sizeof(rawNiftiHeader.qoffset_y),&rawNiftiHeader.qoffset_y);
  byteSwap(sizeof(rawNiftiHeader.qoffset_z),&rawNiftiHeader.qoffset_z);
  byteSwap(sizeof(rawNiftiHeader.srow_x[0]),rawNiftiHeader.srow_x,4);
  byteSwap(sizeof(rawNiftiHeader.srow_y[0]),rawNiftiHeader.srow_y,4);
  byteSwap(sizeof(rawNiftiHeader.srow_z[0]),rawNiftiHeader.srow_z,4);
}

void NiftiIO::byteSwap(const size_t elementLength, void* vBuffer,const unsigned long nElements)
{
  //cerr << "Low level byte swap: " << elementLength << " " << vBuffer << " " << nElements << endl;
  char *buffer(static_cast<char *>(vBuffer));
  for ( unsigned long current = 0; current < nElements; current ++ ) {
    reverse(buffer,buffer+elementLength);
    buffer+=elementLength;
  }
}

template<class T>
void NiftiIO::reportHeader(const T& header)
{
  cout << "dim_info " << (int)header.dim_info << endl;
  for ( int i=0; i<=7; i++ )
    cout << "dim" << i << "\t\t" << header.dim[i] << endl;
  for ( int i=0; i<=7; i++ )
    cout << "pixdim" << i << "\t\t" << header.pixdim[i] << endl;

  cout << "intentp[1..3]\t" << header.intent_p1 << " " << header.intent_p2 << " " << header.intent_p3 << endl;
  cout << "intent_code\t" << header.intent_code << endl;
  cout << "datatype\t" << header.datatype << endl;
  cout << "bitpix\t" << header.bitpix << endl;
  cout << "slice_start\t" << header.slice_start << endl;
  cout << "size of header\t" << header.sizeof_hdr << endl;
  cout << "vox_offset\t" << header.vox_offset << endl;
  cout << "scl_slope\t" << header.scl_slope << endl;
  cout << "scl_inter\t" << header.scl_inter << endl;
  cout << "slice_end\t" << header.slice_end << endl;
  cout << "slice_code\t" << (int)header.slice_code << endl;
  cout << "xyzt_units\t" << (int)header.xyzt_units << endl;
  cout << "cal_max\t" << header.cal_max << endl;
  cout << "cal_min\t" << header.cal_min << endl;
  cout << "slice_duration\t" << header.slice_duration << endl;
  cout << "toffset\t" << header.toffset << endl;
  cout << "descrip\t" << string(header.descrip) << endl;
  cout << "aux_file\t" << string(header.aux_file) << endl;
  cout << "sform_code\t" << header.sform_code << endl;
  cout << "sform:1\t\t" << header.srow_x[0] << " "  << header.srow_x[1] << " " << header.srow_x[2] << " " << header.srow_x[3] << endl;
  cout << "sform:2\t\t" << header.srow_y[0] << " "  << header.srow_y[1] << " " << header.srow_y[2] << " " << header.srow_y[3] << endl;
  cout << "sform:3\t\t"  << header.srow_z[0] << " "  << header.srow_z[1] << " " << header.srow_z[2] << " " << header.srow_z[3] << endl;
    cout << "qform_code\t" << header.qform_code << endl;
    cout << "quatern_b\t" << header.quatern_b << endl;
    cout << "quatern_c\t" << header.quatern_c << endl;
    cout << "quatern_d\t" << header.quatern_d << endl;
    cout << "qoffset_x\t" << header.qoffset_x << endl;
    cout << "qoffset_y\t" << header.qoffset_y << endl;
    cout << "qoffset_z\t" << header.qoffset_z << endl;
  cout << "intent_name\t" << string(header.intent_name) << endl;
  cout << "magic\t" << string(header.magic) << endl;
}

void NiftiIO::reportLegacyHeader(const nifti_1_header& header)
{
  cout << "<legacy> data_type " << string(header.data_type) << endl;
  cout << "<legacy> db_name " << string(header.db_name) << endl;
  cout << "<legacy> extents " << header.extents << endl;
  cout << "<legacy> session_error " << header.session_error << endl;
  cout << "<legacy> regular " << header.regular << endl;
  cout << "<legacy> glmax " << header.glmax << endl;
  cout << "<legacy> glmin " << header.glmin << endl;
}

void NiftiHeader::report() const
{
  cout << "size of header\t" << sizeof_hdr << endl;
  cout << "data_type\t" << datatypeString() << endl;
  for ( int i=0; i<=7; i++ )
    cout << "dim" << i << "\t\t" << dim[i] << endl;
  cout << "vox_units\t" << unitsString(XYZT_TO_SPACE(units)) << endl;
  cout << "time_units\t" << unitsString(XYZT_TO_TIME(units)) << endl;
  cout << "datatype\t" << datatype << endl;
  cout << "nbyper\t\t" << bitsPerVoxel/8 << endl;
  cout << "bitpix\t\t" << bitsPerVoxel << endl;
  cout.setf(ios::fixed);
  for ( int i=0; i<=7; i++ )
    cout << "pixdim" << setprecision(6) << i << "\t\t" << pixdim[i] << endl;
  cout << "vox_offset\t" << vox_offset << endl;
  cout << "cal_max\t\t" << cal_max << endl;
  cout << "cal_min\t\t" << cal_min << endl;
  cout << "scl_slope\t" << sclSlope << endl;
  cout << "scl_inter\t" << sclInter << endl;
  cout << "phase_dim\t" << (int)phaseDim() << endl;
  cout << "freq_dim\t" << (int)freqDim() << endl;
  cout << "slice_dim\t" << (int)sliceDim() << endl;
  cout << "slice_name\t" << niftiSliceString() << endl;
  cout << "slice_code\t" << sliceCode << endl;
  cout << "slice_start\t" << sliceStart << endl;
  cout << "slice_end\t" << sliceEnd << endl;
  cout << "slice_duration\t" << sliceDuration << endl;
  cout << "toffset\t\t" << toffset << endl;
  cout << "intent\t\t" << niftiIntentString() << endl;
  cout << "intent_code\t" << intentCode << endl;
  cout << "intent_name\t" << intentName << endl;
  cout << "intent_p1\t" << intent_p1 << endl;
  cout << "intent_p2\t" << intent_p2 << endl;
  cout << "intent_p3\t" << intent_p3 << endl;
  cout << "qform_name\t" << qFormName() << endl;
  cout << "qform_code\t" << qformCode << endl;
  mat44 output(getQForm());
  for ( int i=0;i<4;i++ ) {
    cout << "qto_xyz:" << char(i+'1') << "\t";
    for ( int j=0;j<4;j++ )
      cout << output.m[i][j] << " ";
    cout << endl;
  }
  int icode,jcode,kcode;
  nifti_mat44_to_orientation(output,&icode,&jcode,&kcode);
  cout << "qform_xorient\t" << niftiOrientationString(icode) << endl;
  cout << "qform_yorient\t" << niftiOrientationString(jcode) << endl;
  cout << "qform_zorient\t" << niftiOrientationString(kcode) << endl;
  cout << "sform_name\t" << sFormName() << endl;
  cout << "sform_code\t" << sformCode << endl;
  output=getSForm();
  for ( int i=0;i<4;i++ ) {
    cout << "sto_xyz:" << char(i+'1') << "\t";
    for ( int j=0;j<4;j++ )
      cout << output.m[i][j] << " ";
    cout << endl;
  }
  nifti_mat44_to_orientation(output,&icode,&jcode,&kcode);
  cout << "sform_xorient\t" << niftiOrientationString(icode) << endl;
  cout << "sform_yorient\t" << niftiOrientationString(jcode) << endl;
  cout << "sform_zorient\t" << niftiOrientationString(kcode) << endl;
  cout << "file_type\t" << fileType() << endl;
  cout << "file_code\t" << (int)niftiVersion() << endl;
  cout << "descrip\t\t" << description << endl;
  cout << "aux_file\t" << auxillaryFile << endl;
}




NiftiHeader::NiftiHeader(vector<string> XMLreport)
{
  initialise();
  int xyzUnits(0);
  int tUnits(0);
  int freqDim(0);
  int phaseDim(0);
  int sliceDim(0);
  for(unsigned int setting=1;setting<XMLreport.size()-1;setting++) {
      //Tokenise
      size_t pos1=XMLreport[setting].find_first_of("  ")+2;
      size_t pos2=XMLreport[setting].find_first_of(" ",pos1);
      size_t pos3=XMLreport[setting].find_first_of('\'')+1;
      size_t pos4=XMLreport[setting].find_last_of('\'');
      string field(XMLreport[setting].substr(pos1,pos2-pos1));
      string value(XMLreport[setting].substr(pos3,pos4-pos3));
      istringstream values(value);
      if ( field == "datatype" ) {
	values >> datatype;
      } else if ( field == "image_offset" ) {
	values >> vox_offset;
      } else if ( field == "sto_xyz_matrix" ) {
	mat44 sForm;
	for ( int i=0;i<4;i++ )
	  for ( int j=0;j<4;j++ )
	    values >> sForm.m[i][j];
	setSForm(sForm);
      } else if ( field == "ndim" ) {
	values >> dim[0];
      } else if ( field == "nx" ) {
	values >> dim[1];
      } else if ( field == "ny" ) {
	values >> dim[2];
      } else if ( field == "nz" ) {
	values >> dim[3];
      } else if ( field == "nt" ) {
	values >> dim[4];
      } else if ( field == "nu" ) {
	values >> dim[5];
      } else if ( field == "nv" ) {
	values >> dim[6];
      } else if ( field == "nw" ) {
	values >> dim[7];
      }else if ( field == "qfac" ) {
	values >> pixdim[0];
      } else if ( field == "dx" ) {
	values >> pixdim[1];
      } else if ( field == "dy" ) {
	values >> pixdim[2];
      } else if ( field == "dz" ) {
	values >> pixdim[3];
      } else if ( field == "dt" ) {
	values >> pixdim[4];
      } else if ( field == "du" ) {
	values >> pixdim[5];
      } else if ( field == "dv" ) {
	values >> pixdim[6];
      } else if ( field == "dw" ) {
	values >> pixdim[7];
      } else if ( field == "cal_min" ) {
	values >> cal_min;
      } else if ( field == "cal_max" ) {
	values >> cal_max;
      } else if ( field == "scl_slope" ) {
	values >> sclSlope;
      } else if ( field == "scl_inter" ) {
	values >> sclInter;
      } else if ( field == "intentCode" ) {
	values >> sclSlope;
      } else if ( field == "intent_p1" ) {
	values >> sclInter;
      } else if ( field == "intent_p2" ) {
	values >> cal_min;
      } else if ( field == "intent_p3" ) {
	values >> cal_max;
      } else if ( field == "intent_name" ) {
	values >> intentName;
      } else if ( field == "toffset" ) {
	values >> toffset;
      } else if ( field == "xyz_units" ) {
	values >> xyzUnits;
      } else if ( field == "time_units" ) {
	values >> tUnits;
      } else if ( field == "descrip" ) {
	values >> description;
      } else if ( field == "aux_file" ) {
	values >> auxillaryFile;
      } else if ( field == "qformcode" ) {
	values >> qformCode;
      } else if ( field == "sformcode" ) {
	values >> sformCode;
      } else if ( field == "quatern_b" ) {
	values >> qB;
      } else if ( field == "quatern_c" ) {
	values >> qC;
      } else if ( field == "quatern_d" ) {
	values >> qD;
      } else if ( field == "quatern_x" ) {
	values >> qX;
      } else if ( field == "quatern_y" ) {
	values >> qY;
      } else if ( field == "quatern_z" ) {
	values >> qZ;
      } else if ( field == "freq_dim" ) {
	values >> freqDim;
      } else if ( field == "phase_dim" ) {
	values >> phaseDim;
      } else if ( field == "slice_dim" ) {
	values >> sliceDim;
      } else if ( field == "slice_code" ) {
	values >> sliceCode;
      } else if ( field == "slice_start" ) {
	values >> sliceStart;
      } else if ( field == "slice_end" ) {
	values >> sliceEnd;
      } else if ( field == "slice_duration" ) {
	values >> sliceDuration;
      } else {
	//cerr << XMLreport[setting] << endl;
	//cerr << "Unknown" << endl;
      }
  }
  units=xyzUnits+tUnits;
  sliceOrdering=freqDim | ( phaseDim << 2 ) | ( sliceDim << 4 );
  bitsPerVoxel=bpvOfDatatype();
}
