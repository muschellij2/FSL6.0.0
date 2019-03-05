/*  Copyright (C) 2004 University of Oxford  */

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

#include <iostream>
#include <fstream>
#include "newimage/newimageall.h"
#include "csv_mesh.h"
#include <vector>
using namespace std;
using namespace NEWIMAGE;

void proj_thresh_volumes(const vector<string>& innames,const float& thresh){
  vector<volume<float> > tmpvec(innames.size());
  volume<float> tmp;
  cout<<"number of inputs "<<innames.size()<<endl;
  for(unsigned int i=1;i<=innames.size();i++){
    cout<<i<<" "<<innames[i-1]<<endl;
    read_volume(tmp,innames[i-1]);
    tmpvec[i-1]=tmp;
  }
  cerr<<"threshold "<<thresh<<endl;

  volume<float> total;
  volume<float> total_thresh;
  total.reinitialize(tmp.xsize(),tmp.ysize(),tmp.zsize());
  total_thresh.reinitialize(tmp.xsize(),tmp.ysize(),tmp.zsize());
  copybasicproperties(tmp,total_thresh);
  copybasicproperties(tmp,total);
  total=total*0;
  for(unsigned int i=0;i<tmpvec.size();i++){
    total+=tmpvec[i];
  }
  
  total_thresh=binarise(total,thresh);
  total.setDisplayMaximumMinimum(total.max(),total.min());
  save_volume(total,"total");
  for(unsigned int i=0;i<tmpvec.size();i++){
    tmp=divide(tmpvec[i],total,total_thresh);
    string outname =innames[i];
    make_basename(outname);
    string thrname="_thr_"+num2str(thresh);
    total.setDisplayMaximumMinimum(1,0);
    save_volume(tmp,outname+"_proj_seg"+thrname);
  }
}
void proj_thresh_surfaces(const vector<string>& innames,const float& thresh){
  vector<CsvMesh> meshes;
  CsvMesh m;
  cout<<"number of inputs "<<innames.size()<<endl;
  for (unsigned int i=0;i<innames.size();i++){
    m.load(innames[i]);
    meshes.push_back(m);
  }
  cerr<<"threshold "<<thresh<<endl;
  CsvMesh tot;
  tot=m;
  for (int i=0;i<m.nvertices();i++){
    float totval=0;
    for (unsigned int j=0;j<innames.size();j++){
      totval += (meshes[j].get_pvalue(i)>thresh?meshes[j].get_pvalue(i):0);
    }
    tot.set_pvalue(i,totval);
    for (unsigned int j=0;j<innames.size();j++){
      if(meshes[j].get_pvalue(i)<=thresh || totval==0)
	meshes[j].set_pvalue(i,0);
      else{
	meshes[j].set_pvalue(i,meshes[j].get_pvalue(i)/totval);
      }
    }
  }
  tot.save("total",meshFileType(innames[0]));
  for(unsigned int i=0;i<innames.size();i++){    
    string outname =innames[i];
    make_basename(outname);
    string thrname="_thr_"+num2str(thresh);
    meshes[i].save(outname+"_proj_seg"+thrname,meshFileType(innames[i]));
  }
}
bool test_input(const vector<string>& filenames){
  int nsurfs=0,nvols=0;
  for(unsigned int i=0;i<filenames.size();i++){
    if(meshExists(filenames[i])){nsurfs++;}
    else if(fsl_imageexists(filenames[i])){nvols++;}
    else{
      cerr<<"Cannot open "<<filenames[i]<<" for reading"<<endl;
      exit(1);
    }
  }
  if(nvols>0 && nsurfs>0){
    cerr<<"Please use list of EITHER volumes OR surfaces"<<endl;
    exit(1);
  }
  return (nvols>0?true:false);
}
int main ( int argc, char **argv ){
  if(argc<3){
    cerr<<""<<endl;
    cerr<<"usage:proj_thresh <lots of volumes/surfaces> threshold"<<endl;
    cerr<<"  please use EITHER volumes OR surfaces"<<endl;
    cerr<<""<<endl;
    exit(1);
  }

  vector<string> innames;
  innames.clear();
  for(int i=1;i<=argc-2;i++){
    innames.push_back(argv[i]);
  }
  bool isvols=test_input(innames);

  float thresh=atof(argv[argc-1]);
  if(!isvols){
    proj_thresh_surfaces(innames,thresh);
  }
  else{
    proj_thresh_volumes(innames,thresh);
  }

  
}










