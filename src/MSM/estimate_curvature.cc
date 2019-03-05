/*  estimate_curvature.cc

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
#include "newmat.h"
#include "newmatio.h"
#include "utils/options.h"
#include "newimage/newimageall.h"


#include "MeshReg/meshreg.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "estimate_curvature  <anatomical mesh>  <output base> <-option>   " << endl;
  cout << " estimates mean curvature by fiting a least-squares quadratic patch to the local neighborhood of a vertex f(x,y) = ax^2 + by^2 + cxy + dx + ey + f (inspired by patchcuvature.m the matlab function " << endl;
  cout << " options: " << endl;
  cout << " -fit-radius default 2 " << endl; 
  cout << " " << endl;
}


int main(int argc, char **argv){

  
  newmesh ORIG,SPHERE;
  string outdir;
  double fit_radius=2, smoothval=2.0;
  boost::shared_ptr<RELATIONS> REL;
  bool _calcrel=false,_smooth=false;
  bool _outputrho=false;
  int ok;

  if(argc < 3){

    Usage();
    exit(0);
  }

 
  ORIG.load(argv[1]);
  argc--; 
  argv++;
  outdir=argv[1];
  argc--; 
  argv++;


  while (argc > 1) {
    ok = 0;
   
    if((ok == 0) && (strcmp(argv[1], "-fit-radius") == 0)){
       argc--; 
       argv++;
       fit_radius=atof(argv[1]);
       argc--; 
       argv++;
       ok=1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-sphere") == 0)){
       argc--; 
       argv++;
       SPHERE.load(argv[1]);
       _calcrel=true;
       argc--; 
       argv++;
       ok=1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-smooth") == 0)){
       argc--; 
       argv++;
       _smooth=true;
       smoothval=atof(argv[1]);
       argc--; 
       argv++;
       ok=1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-rho") == 0)){
       argc--; 
       argv++;
       _outputrho=true;
       ok=1;
    } else{cout << " option doesn't exist " << endl; exit(1);}
  }

 if(_calcrel){
    REL=boost::shared_ptr<RELATIONS >(new RELATIONS (SPHERE,SPHERE,3*asin(fit_radius/RAD))); 
    REL->update_RELATIONS(SPHERE);
  }
  cout << " estimate curvature " << endl;
  mean_curvature(fit_radius,ORIG,REL);

  if(_smooth){
    if(!_calcrel) {cout << " need sphere for smoothing currently " << endl; exit(1);}
    resampler R;
    boost::shared_ptr<BFMatrix> scalardata;
    boost::shared_ptr<newmesh> EXCL;
    scalardata =boost::shared_ptr<BFMatrix> (new FullBFMatrix (1,SPHERE.nvertices()));
    R.set_method("GAUSSIAN");

    for (int i=0;i<ORIG.nvertices();i++)
      scalardata->Set(1,i+1,ORIG.get_pvalue(i));

    R.resampledata(SPHERE,SPHERE,EXCL,scalardata,smoothval,REL)  ;

    for (int i=0;i<ORIG.nvertices();i++)
      ORIG.set_pvalue(i,scalardata->Peek(1,i+1));

  }
  char filename[1000];
  sprintf(filename,"%s-meancurvature.func",outdir.c_str());
  ORIG.save(filename);
  if(_outputrho){

    for (int i=0;i<ORIG.nvertices();i++)
      ORIG.set_pvalue(i,1/ORIG.get_pvalue(i));

    sprintf(filename,"%s-rho.func",outdir.c_str());
    ORIG.save(filename);
  }

 
}
