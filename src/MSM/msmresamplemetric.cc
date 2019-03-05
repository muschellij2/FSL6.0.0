/*  msmresamplemetric.cc

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

#include "newmesh/newmesh.h"
#include "MeshReg/meshreg.h"
#include "miscmaths/SpMat.h"

using namespace MESHREG;

void Usage()
{
  cout << "msmresamplemetric  <inputmesh> <output> <-option>   " << endl;
  cout << "-option " << endl;
  cout << " -project - project final result back down onto surface (requires argument) " << endl;
  cout << " -labels - load labels file for input (inc. .func. and .shape, requires argument)" << endl;
  cout << " -datamat - multivariate data matrix (requires argument)" << endl; 
  cout << " -barycentric use barycentric interpolation" << endl;
  cout << " -adap_bary use adaptive barycentric interpolation" << endl;
  cout << " -linear - linear interpolation kernel with kernel size X (requires argument)" << endl;
  cout << " -gaussian - gaussian interpolation kernel with std deviation X (requires argument) " << endl;
  cout << " -normalize  - normlalize intensity range to target (requires argument)" << endl;
  cout << " -excl exclude the area of the cut form contributing to the resampling" << endl;
  cout << " -arealdist estimate areal distortion" << endl;
}


int main(int argc, char **argv){

  Matrix M;
  int ok,getwm;
  bool barycentric,adap_barycentric,_normalize;
  resampler resample;
  //MeshReg MR; // only used for read list -> move to resample?
  NEWMESH::newmesh in,SPH,REG,ICO,wm,HISTTARG;
  bool _outputlabel,_transform,_project,_exclude;
  bool _save_rel,_datamatrix,getdist, _datamatsp;
  string output;
  double thr;
  double sigma,rad;
  boost::shared_ptr<RELATIONS>  _rel;
  string rel_out;
  boost::shared_ptr<BFMatrix > datamat;
  bool labels=false;
  // BFMatrix *datamat;
  

  if(argc < 3){

    Usage();
    exit(0);
  }
 
  
  in.load(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;

  recentre(in);
  
  getwm=0; thr=0.0; sigma=1.0;rad=100.0;
  _transform=false;_project=false;_outputlabel=false;
  _datamatrix=false; _save_rel=false;  _datamatsp=false;
  _exclude=false; resample.set_method("NN");
  barycentric=false; adap_barycentric=false;  _normalize=false;
  
  // cout << " here 2 " << endl;
  while (argc > 1) {
    ok = 0;
    cout << argv[1] << endl;
   
    if((ok == 0) && (strcmp(argv[1], "-project") == 0)){
      argc--;
      argv++;
      _project=true;
      SPH.load(argv[1]);
      cout << " recentre target " << endl;
      recentre(SPH);
      argc--;
      argv++;
      ok = 1;
    }else if((ok == 0) && (strcmp(argv[1], "-labels") == 0)){
      argc--;
      argv++;
      labels=true;
      cout << " labels " << argv[1] << endl;
      in.load(argv[1],false,false);
      Matrix M=in.get_pvalues();
      datamat = boost::shared_ptr<BFMatrix >(new FullBFMatrix (M));

      argc--;
      argv++;
      ok = 1;
    }
  else if((ok == 0) && (strcmp(argv[1], "-datamat") == 0)){
      argc--;
      argv++;
      _datamatrix=true;
      Matrix M=read_ascii_matrix(argv[1]);
      datamat = boost::shared_ptr<BFMatrix >(new FullBFMatrix (M));
      in.set_pvalues(M);
      argc--;
      argv++;
      ok = 1;
    }
  else if((ok == 0) && (strcmp(argv[1], "-datamatsp") == 0)){
      argc--;
      argv++;
      _datamatsp=true;
 cout << " load sparsemat " << endl;
      SpMat<double> M(argv[1]);
      //   tmp=new SpMat<double> (argv[1]);
      cout << " after load sparsemat " << endl;
      datamat = boost::shared_ptr<BFMatrix >(new SparseBFMatrix<double> (M));
      in.set_pvalues(datamat->AsMatrix());
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-normalize") == 0)){
      argc--;
      argv++;
      _normalize=true;
      cout <<" normalize " << endl;
      HISTTARG.load_gifti(argv[1],false,false);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-wmmesh") == 0)){
      argc--;
      argv++;
      getwm=1;
      wm.load(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }   
    else if((ok == 0) && (strcmp(argv[1], "-excl") == 0)){
      argc--;
      argv++;
      _exclude=true;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-barycentric") == 0)){
      argc--;
      argv++;
      barycentric=true;
      resample.set_method("BARYCENTRIC");  
      cout << " METHOD IS BARCENTRIC " << endl;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-adap_bary") == 0)){
      argc--;
      argv++;
      adap_barycentric=true;
      resample.set_method("ADAP_BARY");  
      cout << " METHOD IS ADAPTIVE BARCENTRIC " << endl;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-linear") == 0)){
      argc--;
      argv++;
      resample.set_method("LINEAR");
      sigma=atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-gaussian") == 0)){
      argc--;
      argv++;
      resample.set_method("GAUSSIAN");
      sigma=atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    } 
    else{cout << " option doesn't exist " << endl; exit(1);}

   
  }

  
  double meannorm=0;
  double varnorm=0;
 
  boost::shared_ptr<NEWMESH::newmesh> EXCL_IN, EXCL_REF;
  Matrix DATAIN;
  if(_exclude){
    DATAIN=in.get_pvalues();
    EXCL_IN= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(in,DATAIN,0,0))); 
    // EXCL_IN->save("IN_EXCL.func");
    //DATAIN->Print("DATAIN");
   
  }
 
  if(_normalize)  {
    Matrix DATAREF;
    boost::shared_ptr<BFMatrix > BFIN,BFREF;
    DATAREF=HISTTARG.get_pvalues();
    if(!datamat.get()){
      Matrix M=in.get_pvalues();
      datamat = boost::shared_ptr<BFMatrix >(new FullBFMatrix (M));

    }

    if(_exclude){
      cout << "exclude " << endl;
      EXCL_REF= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(SPH,DATAREF,0,0))); 
      // EXCL_REF->save("REF_EXCL.func");
    }
    
    BFREF=boost::shared_ptr<BFMatrix > (new FullBFMatrix (DATAREF)); //
    multivariate_histogram_normalization(*datamat,*BFREF,EXCL_IN,EXCL_REF,false);
    in.set_pvalues(datamat->AsMatrix());

  }

 


  bool _scalar=false;

  if(_project){
    for (int i=0;i<SPH.nvertices();i++)
      SPH.set_pvalue(i,0); 
    
    if(datamat.get()){
      resample.resampledata(in,SPH,EXCL_IN,datamat,sigma); 
      //SPH.set_pvalues(datamat->AsMatrix());
           
    }
    else{
      resample.resample_scalar(in,SPH,sigma,EXCL_IN);
      //resample.downsample(in,SPH,sigma,EXCL);
    }
  

    if(labels){
      SPH.set_pvalues(datamat->AsMatrix());
      SPH.save(output+".func");
    }else{
      datamat->Print(output);
   
    }
    
   

    if(getwm){
      
      for (int i = 0; i < SPH.nvertices(); i++)
	wm.set_pvalue(i,SPH.get_pvalue(i));
	
    }

 

    if(getwm)
      wm.save(output +"targetmesh");
  
  }
  //delete datamatsp;
}
