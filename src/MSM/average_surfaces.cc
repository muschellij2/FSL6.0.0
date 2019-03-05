/*  average_surfaces.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

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
/* this program is designed to downsample freesurfer label files to be used in combination with the SPH6.vtk or other downsampled meshes*/

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include "utils/options.h"
#include "newmesh/meshfns.h"
#include "MeshReg/meshreg.h"


using namespace std;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "average_surfaces  <asciii list of data> <output> <-options>   " << endl;
  cout << "options:  " << endl;
  cout << "-target " << endl;
  cout << "-abs " << endl;
  cout << "-weighted distances" << endl;
  cout << "-simga " << endl;

}


int main(int argc, char **argv){


  boost::shared_ptr<NEWMESH::newmesh> EXCL_ref;
  SpMat<double> *DATA =new SpMat<double>();
  int ok;
  int _avdim=0;
  double newval;
  float _sigma;
  
  NEWMESH::newmesh totalcounts;
 
  NEWMESH::newmesh IN, TARG, HISTTARG,AVERAGE,VARIANCE,Zscore;
  
  string output,excl_in,excl_ref;
  vector<string> distances;

  vector<string> INlist,INsurflist;
  vector<newmesh> INPUTS;
  bool _isdatamatrix=false;
  string subtype;
  bool _exclude=false;
  bool _abs=false,_weighted=false;
  ColumnVector exclusion;

  if(argc < 2){

    Usage();
    exit(0);
  }

  INlist=read_ascii_list(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;
  newval=0.0;
  excl_in="";excl_ref="";
  EXCL_ref=boost::shared_ptr<NEWMESH::newmesh >(new NEWMESH::newmesh ()) ;
 

  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-target") == 0)){
      argc--;
      argv++;
      _isdatamatrix=1;
      cout << " normalize " << argv[1] << endl;
      TARG.load(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-inputsurfaces") == 0)){
      argc--;
      argv++;
      INsurflist=read_ascii_list(argv[1]);
      if(INsurflist.size()!=INlist.size()){ cout << " surf list is not the same length as data list" << endl; exit(1); }
      argc--;
      argv++;
      ok = 1;
    }
	else if((ok == 0) && (strcmp(argv[1], "-excl_in") == 0)){
      argc--;
      argv++;
      _exclude=1;
      excl_in=argv[1];
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-excl_ref") == 0)){
      argc--;
      argv++;
      _exclude=1;
      excl_ref=argv[1];
      EXCL_ref->load(excl_ref,false,false);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-abs") == 0)){
      argc--;
      argv++;
      _abs=1;
      ok = 1;
    }
     else if((ok == 0) && (strcmp(argv[1], "-averagedim") == 0)){
      argc--;
      argv++;
      _avdim=atoi(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
   else if((ok == 0) && (strcmp(argv[1], "-weighted") == 0)){
      argc--;
      argv++;
      _weighted=true;
      distances=read_ascii_list(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-sigma") == 0)){
      argc--;
      argv++;
      _sigma=atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else{cout << " option doesn't exist " << endl; exit(1);}
  }
  if(!_isdatamatrix)TARG.load(INlist[0]);
  if(excl_ref==""){
    *EXCL_ref=TARG;

    for (int i=0;i<EXCL_ref->nvertices();i++){
      EXCL_ref->set_pvalue(i,1);
    }
  }

  exclusion.ReSize(TARG.nvertices()); exclusion=1;

  /////////////// INITIALIZE //////////////////////
  int dim,vert;
  subtype = INlist[0].substr(INlist[0].size()-8, 4);
  if(_isdatamatrix){  /// set average either to dimensions of target mesh or, assuming all data have been resampled already to the size of first input mesh
    AVERAGE=TARG;
    vert=TARG.nvertices();
    if(INsurflist.size()==0) IN=AVERAGE;
  }
  else{
    AVERAGE.load(INlist[0]);
    vert=AVERAGE.nvertices();
  }
  
  if(subtype=="func" || subtype=="hape" ){ // define dimension of data
	IN.load(INlist[0],false,false);
	dim=IN.get_dimension();
  }else {
     DATA=new SpMat<double>(INlist[0]);
     if(DATA->Ncols()>DATA->Nrows()) dim=DATA->Nrows();
     else dim=DATA->Ncols();
  }
  

  /// initialise average and variances to zero 
  Matrix tmpData(dim,vert);  tmpData=0;
  AVERAGE.set_pvalues(tmpData);
  VARIANCE=AVERAGE;
  Zscore=AVERAGE;   
  totalcounts=AVERAGE;
  
  int atlasdim=AVERAGE.get_dimension();
  double sumweighted=0;
  vector<float> weight(INlist.size(),1);
  /////////////////////////////////////////////////////
  for (unsigned int i=0;i<INlist.size();i++){
	  cout << i << " " << INlist[i] << endl;
    if(subtype=="func"|| subtype=="hape"){
	  if(INsurflist.size()==INlist.size()) IN.load(INsurflist[i]);
      IN.load(INlist[i],false,false);
    }else if(_isdatamatrix){
      DATA=new SpMat<double>(INlist[i]);
      IN.set_pvalues(DATA->AsNEWMAT());
    }else
      IN.load(INlist[i]);

    if(IN.nvertices()!=TARG.nvertices() && TARG.nvertices()>0){
      //// resample data onto target
      resampler R; R.set_method("ADAP_BARY");
      Matrix tmpDATA=IN.get_pvalues();
      R.resampledata(IN,TARG,tmpDATA,1);
      IN=TARG;
      IN.set_pvalues(tmpDATA);
    }

	INPUTS.push_back(IN);
    boost::shared_ptr<NEWMESH::newmesh> EXCL;
    int dim=IN.get_dimension();
    if(dim!=atlasdim){cout << i << " dim " << dim << " atlasdim " << endl; throw MeshReg_error("Feature dimensions do not agree");}
    
    if(_weighted){
      if(_sigma==0){ cout << "must set sigma for weighted" <<endl; exit(1);}
      weight[i]=exp(-0.5*((atof(distances[i].c_str()))/_sigma)*((atof(distances[i].c_str()))/_sigma));
     // sumweighted+=weight[i];
      cout << "distance " <<  distances[i] << " weight " << weight[i] << endl;
    }//else{ sumweighted=1;}
      cout << i << " weight " << weight[i] << endl;
    // cout << " add to total ... " << IN.nvertices() << endl;
    for(int d=0;d<atlasdim;d++){
      cout << d << " atlasdim " <<atlasdim << " dim " << dim << "EXCL_ref->nvertices() " << EXCL_ref->nvertices() << endl; 
      for (int j=0;j<AVERAGE.nvertices();j++){
	if(EXCL_ref->get_pvalue(j) && abs(IN.get_pvalue(j,d)) > EPSILON){
	  if(_abs) newval=AVERAGE.get_pvalue(j,d)+weight[i]*abs(IN.get_pvalue(j,d));
	  else newval=AVERAGE.get_pvalue(j,d)+weight[i]*IN.get_pvalue(j,d);
	  totalcounts.set_pvalue(j,totalcounts.get_pvalue(j,d)+weight[i]*1,d);
	  AVERAGE.set_pvalue(j,newval,d);
	}else{
	  if(_exclude)
	    exclusion(j+1)=0;
	}

      }
     
    }
    // cout <<" here "<< endl;
  }

  double sumav=0;
 
  for(int d=0;d<atlasdim;d++){  
    for (int i=0;i<AVERAGE.nvertices();i++){
      if(totalcounts.get_pvalue(i,d)==0) newval=0;
      else newval=AVERAGE.get_pvalue(i,d)/(totalcounts.get_pvalue(i,d));
      AVERAGE.set_pvalue(i,newval,d);
      if(_avdim==-1 || _avdim==d) sumav+=newval;
      //   cout << i << " newval " << newval << " AVERAGE.get_pvalue(d,i) " <<  AVERAGE.get_pvalue(d,i) << " totalcounts.get_pvalue(d,i) " << totalcounts.get_pvalue(d,i) << " sumav " << sumav << endl;
    }
  }
  

  for (unsigned int i=0;i<INlist.size();i++){
   
    for(int d=0;d<atlasdim;d++){
      for (int j=0;j<AVERAGE.nvertices();j++){
	if(EXCL_ref->get_pvalue(j) && exclusion(j+1)){
	  if(_abs) newval=weight[i]*(abs(INPUTS[i].get_pvalue(j,d))-AVERAGE.get_pvalue(j,d));
	  else  newval=weight[i]*(INPUTS[i].get_pvalue(j,d)-AVERAGE.get_pvalue(j,d));
	  newval=VARIANCE.get_pvalue(j,d)+newval*newval;
	  VARIANCE.set_pvalue(j,newval,d);
	}
      }

    }
  }
  double meanvar=0;
  for(int d=0;d<atlasdim;d++){
    for (int i=0;i<VARIANCE.nvertices();i++){
      if(EXCL_ref->get_pvalue(i) && totalcounts.get_pvalue(i,d) > 0 && exclusion(i+1)){
	newval=sqrt(VARIANCE.get_pvalue(i,d)/(totalcounts.get_pvalue(i,d)));
	VARIANCE.set_pvalue(i,newval,d);
	if(_avdim==-1 || _avdim==d) meanvar+=newval;
	//cout << i << " newval " << newval << " VARIANCE.get_pvalue(d,i) " <<  VARIANCE.get_pvalue(d,i) << " AVERAGE.get_pvalue(d,i) " <<  AVERAGE.get_pvalue(d,i) << " totalcounts.get_pvalue(d,i) " << totalcounts.get_pvalue(d,i) << " meanvar " << meanvar << endl;

      }
      else VARIANCE.set_pvalue(i,0,d);
    }
  }

  for(int d=0;d<atlasdim;d++){
    for (int i=0;i<VARIANCE.nvertices();i++){
      if(VARIANCE.get_pvalue(i,d)>0 && exclusion(i+1))
	newval=(abs(AVERAGE.get_pvalue(i,d)))/VARIANCE.get_pvalue(i,d);
      else
	newval=0;

      Zscore.set_pvalue(i,newval,d);
    }
  }

int totdim;
  if(_avdim==-1) totdim=atlasdim;
  else totdim=1;
  cout << " meanval " << sumav/(VARIANCE.nvertices()*totdim) << endl;

  cout << " meanvar " << meanvar/(VARIANCE.nvertices()*totdim) << endl;

  AVERAGE.save(output+"AVERAGE.func");
  VARIANCE.save(output+"STD.func");
  Zscore.save(output+"ZSCORE.func");
  delete DATA;
}
