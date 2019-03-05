/*  Copyright (c) 2016 University of Oxford  */
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

#include "gsa.h"

using namespace cifti;
using namespace NEWIMAGE;

namespace GSA {

  template<class T>
  ReturnMatrix GeneralSpatialAbstractor<T>::abstractFile(const string&filename,const string& mode) {
    Matrix tmpData;
    sourceType=mode;
    if ( sourceType.compare("CIFTI") == 0 ) {
      cifti::CiftiFile inputCifti;
      inputCifti.openFile(make_basename(filename)+".nii");
      ciftiExemplar=inputCifti.getCiftiXML();
      cerr << "ndim " << ciftiExemplar.getNumberOfDimensions() << endl;
      cerr << "type1 " << ciftiExemplar.getMappingType(0) << endl;
      cerr << "type2 " << ciftiExemplar.getMappingType(1) << endl;
      CiftiBrainModelsMap spatialModels(ciftiExemplar.getBrainModelsMap(1));
      std::vector<CiftiBrainModelsMap::ModelInfo> models(spatialModels.getModelInfo());
      std::vector<StructureEnum::Enum> nSurfaces(spatialModels.getSurfaceStructureList());
      cerr << "Num surfaces " << nSurfaces.size() << endl;
      cerr << "Has volume " << spatialModels.hasVolumeData() << endl;
      //for(unsigned int currentModel=0;currentModel<models.size();currentModel++)
      //cerr << StructureEnum::toName(models[currentModel].m_structure) << endl;
      const vector<int64_t>& dims = inputCifti.getDimensions();
      tmpData.ReSize(dims[0],dims[1]); //swapped compared to cifti
      vector<float> scratchRow(dims[0]);//read/write a row at a time
      for (int64_t row=0;row<dims[1];row++) {
	inputCifti.getRow(scratchRow.data(),row);
	for (int64_t col=0;col<dims[0];col++) 
	  tmpData(col+1,row+1)=scratchRow[col];
      } 
    }
    if ( sourceType.compare("NIFTI") == 0 ) {
      volume4D<T> inputNifti;
      NEWIMAGE::read_volume4D(inputNifti,filename);
      niftiGeometry=inputNifti[0];
      storageMap.resize(niftiGeometry.xsize()*niftiGeometry.ysize()*niftiGeometry.zsize()+1);
      list<uniqueNode> allVoxels;
      int64_t index(1);
      for (int64_t z=0; z<niftiGeometry.zsize(); z++)
	for (int64_t y=0; y<niftiGeometry.ysize(); y++) 
	  for (int64_t x=0; x<niftiGeometry.xsize(); x++) {
	    storageMap[index].key.push_back(x);
	    storageMap[index].key.push_back(y);
	    storageMap[index].key.push_back(z);
	    storageMap[index].ID=index;
	    allVoxels.push_back(storageMap[index++]);
	}
      //VolumeGeometry<T> *myGeometry=new VolumeGeometry<T>(niftiGeometry,storageMap,allVoxels);
      //tmpData.ReSize(inputNifti.tsize(),myGeometry->nNodes());
      //myGeometry->updateMatrix(inputNifti,tmpData);
    }
    tmpData.Release();
    return tmpData; 
}


  template<class T>
  void VolumeGeometry<T>::maskMap(const NEWIMAGE::volume<T>& mask) {
    //Remove any known nodes which are outside of mask and tag in storageMap
    for (std::list<uniqueNode>::iterator iNode=this->myNodes.begin();iNode!=this->myNodes.end();iNode++) 
      if (mask(iNode->key[0],iNode->key[1],iNode->key[2])<mask.maskThreshold()) {
	this->dataMap[iNode->ID].ID=-1;
	iNode=this->myNodes.erase(iNode);
      }
  }
  
  template<class T>
  ReturnMatrix VolumeGeometry<T>::updateMapAndMatrix(const NEWIMAGE::volume<T>& mask, const Matrix& data) {
    volume4D<T> before=asVolume(data);
    updateMap(mask);
  }
  
  template<class T>
  NEWIMAGE::volume4D<T> VolumeGeometry<T>::asVolume(const Matrix& data) {
    volume4D<T> output(geometry.xsize(),geometry.ysize(),geometry.zsize(),data.Nrows());
    output.copyproperties(geometry);                
    output=0;
    for (std::list<uniqueNode>::iterator iMap=this->dataMap.begin();iMap!=this->dataMap.end();iMap++)
      if (this->dataMap[iMap->ID]>0)
	for(int64_t t=0;t<data.Nrows();t++)
	  output(iMap->key[0],iMap->key[1],iMap->key[2],t)=data(t+1,this->dataMap[iMap->ID]);
  }


  template<class T>
  void VolumeGeometry<T>::updateMatrix(const NEWIMAGE::volume4D<T>& source, NEWMAT::Matrix& target) {
    for (std::list<uniqueNode>::iterator iMap=this->myNodes.begin();iMap!=this->myNodes.end();iMap++)
      if (this->dataMap[iMap->ID]>0)
	for(int64_t t=0;t<source.tsize();t++)
	  target(t+1,this->dataMap[iMap->ID])=source(iMap->key[0],iMap->key[1],iMap->key[2],t);
  }
  
  template<class T>
  void GeneralSpatialAbstractor<T>::remapData(Matrix& data) {
    //remap data(after masking)
    int nValid(0);
    for(vector<uniqueNode>::iterator iNode=storageMap.begin();iNode!=storageMap.end();iNode++)
      if(iNode->ID>0)
	nValid++;
    if(nValid==data.Ncols()) //No masking has occurred do not modify data
      return;
    Matrix newData(data.Nrows(),nValid);
    int newColumn(1);
    int oldColumn(1);
    for(vector<uniqueNode>::iterator iNode=storageMap.begin();iNode!=storageMap.end();iNode++) {
      if(iNode->ID>0) //Copy
	newData.Column(newColumn++)=data.Column(oldColumn++);
      else if (iNode->ID==-1) {//skip next Column in old data
	oldColumn++;
	iNode->ID=0;
      }
    }
    data=newData;
  }
  
  template<class T>
  void GeneralSpatialAbstractor<T>::saveMatrix(const Matrix& data, const string& filename) {
    if ( sourceType.compare("CIFTI") == 0 ) {
      char intentNameOut[16];
      int intent=ciftiExemplar.getIntentInfo(ciftiExemplar.getParsedVersion(),intentNameOut);
      string extension("");
      CiftiMappingType::MappingType type=ciftiExemplar.getMappingType(0);
      if ( intent == 3001 || intent == 3002 )
	extension=".dscalar";
      if ( intent == 3003 || intent == 3004 )
	extension=".pscalar";
      cifti::CiftiScalarsMap scalarsMap;
      scalarsMap.setLength(data.Nrows());
      ciftiExemplar.setMap(0, scalarsMap);
      CiftiFile outputFile;
      outputFile.setWritingFile(make_basename(filename)+extension+".nii");//sets up on-disk writing with default writing version
      outputFile.setCiftiXML(ciftiExemplar,false);	    
      vector<float> scratchRow(data.Nrows());//read/write a row at a time
      for (int64_t row=0;row<data.Ncols();row++) {	    
	for (int64_t col=0;col<data.Nrows();col++) 
	  scratchRow[col]=data(col+1,row+1);
	outputFile.setRow(scratchRow.data(),row);
      }
    }
    else
      cerr << "GeneralSpatialAbstractor cannot save mode: " << sourceType << endl;
  }
  
  template class GeneralSpatialAbstractor<float>;

}
