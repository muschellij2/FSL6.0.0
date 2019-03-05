/*  giftiInterface.h
    Matthew Webster (FMRIB)
    Copyright (C) 2013 University of Oxford  */
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

#include <stdlib.h>
#include <map>
#include <vector>
#include <iostream>
extern "C" {
#include <giftiio/gifti_io.h>
}

 struct GiftiException : public std::exception
{
  std::string errorMessage;
  GiftiException(const std::string& error) : errorMessage(error) {}
  ~GiftiException() throw() {}
  const char* what() const throw() { return errorMessage.c_str(); }
};

struct GIFTImeta {
  GIFTImeta(const char* inputName, const char* inputValue);
  GIFTImeta(const std::string inputName, const std::string inputValue);

  std::string name;
  std::string value;
};

struct GIFTIcoordinateSystem {
  GIFTIcoordinateSystem(const std::string& inputDataSpace, const std::string& inputTransformSpace, const std::vector<double>& inputTransform) : dataSpace(inputDataSpace), transformSpace(inputTransformSpace), transform(inputTransform){};
  std::string dataSpace;
  std::string transformSpace;
  std::vector<double> transform;
};

class GIFTIfield { //NOT templated by field type as I want to store all fields in single container
  friend class GIFTIwrapper;
  int intent;
  std::vector<int> dims; // size of this is number of dims
  int dataType;
  int ordering;
  std::vector<char> bdata; //Stores data for byte fields
  std::vector<int>  idata; //Stores data for int fields
  std::vector<float> fdata; //Stores data for float fields;
  void swapOrdering();

  std::vector<GIFTIcoordinateSystem> coordSystems;
  std::vector<GIFTImeta> metaData;
  std::vector<GIFTImeta> extraAttributes; //Any other attributes that the tag has - note these are _not_ saved by the GIFTI library despite being in the struct!!!
  long long externalOffset;
  std::string externalFilename;

public:

  GIFTIfield(const int intent, const int datatype, const int nDim,const int* dims,const void* data, const int ordering, const std::vector<GIFTIcoordinateSystem>& inputCoordSystems=std::vector<GIFTIcoordinateSystem>(), const std::vector<GIFTImeta>& inputMeta=std::vector<GIFTImeta>(), const std::vector<GIFTImeta>& inputExtraAttribures=std::vector<GIFTImeta>(), const long long inputExternalOffset=0, const std::string& inputExternalFilename=std::string());
  GIFTIfield(const giiDataArray* fieldPointer );
  void report() const;
  void printAsSixTensor() const;
  int getDim(const unsigned char dim) const;
  int getDataType(void) const;
  int getIntent(void) const;

  std::vector<GIFTIcoordinateSystem> getCoordSystems() const {return coordSystems;};
  std::vector<GIFTImeta> getMetaData() const {return metaData;};
  std::vector<GIFTImeta> getExtraAttributes() const {return extraAttributes;};
  long long getExternalOffset() const {return externalOffset;};
  std::string getExternalFilename() const {return externalFilename;};


  float fScalar(const size_t location) const; //Returns a uni-dimensional float fields value at location, throws exception for other types
  int iScalar(const size_t location) const; //Returns a uni-dimensional  int fields value at location, throws exception for other types
  char bScalar(const size_t location) const; //Returns a uni-dimensional byte fields value at location, throws exception for other types
  float asFscalar(const size_t location) const; //Returns a uni-dimensional float fields value at location, throws exception for other types
  int asIscalar(const size_t location) const; //Returns a uni-dimensional  int fields value at location, throws exception for other types
  char asBscalar(const size_t location) const; //Returns a uni-dimensional byte fields value at location, throws exception for other types
  void setFScalar(const size_t location, const float value);
  void setIScalar(const size_t location, const int value);
  void setBScalar(const size_t location, const char value );
  std::vector<float> fVector(const size_t location) const;  //Returns 2-dimensional vector field at location, throws exception for other types
  std::vector<int> iVector(const size_t location) const;  //Returns 2-dimensional vector field at location, throws exception for other types
  std::vector<char> bVector(const size_t location) const;  //Returns 2-dimensional vector field at location, throws exception for other types
};

class GIFTIlabel {
public:
  std::vector<float> RGBA;
  std::string name;
};

class GIFTIwrapper {
public:
  std::vector<GIFTIfield> allFields;
  std::vector<GIFTImeta> metaData;
  std::vector<GIFTImeta> extraAttributes; //Any other attributes that the tag has
  std::map<int,GIFTIlabel> GIFTIlabels;
  int readGIFTI(const std::string filename);
  int writeGIFTI(const std::string filename, int encoding=GIFTI_ENCODING_ASCII) const;
  std::vector<GIFTIfield> returnSurfaceFields();
  std::vector<GIFTIfield> returnNonSurfaceFields();
  void report() const;
};
