/*  estimate_metric_distortion.cc

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
  cout << "estimate_metric_distortion  <original mesh> <transformed mesh > < output> <-option>   " << endl;
  cout << " Compares triangle areas before and after registration - can be used to assess strength of regularisation" << endl;
  cout << " options: " << endl;
  cout << " -target project result onto a target image (requires argument)" << endl;
  cout << " -abs take absolute value " << endl;
  
}


int main(int argc, char **argv){

  
  newmesh ORIG, TRANSFORMED,TARG,DISTMAP;
  resampler R;
  string outdir;
  boost::shared_ptr<NEWMESH::newmesh> EXCL;
  int ok;
  bool _targ=false,_abs=false;
  if(argc < 4){

    Usage();
    exit(0);
  }

 
  ORIG.load(argv[1]);
  argc--; 
  argv++;
  TRANSFORMED.load(argv[1]);
  argc--; 
  argv++;
  outdir=argv[1];
  argc--; 
  argv++;

  TARG=TRANSFORMED;

  while (argc > 1) {
    ok = 0;
   
    if((ok == 0) && (strcmp(argv[1], "-target") == 0)){
      argc--;
      argv++;
      _targ=true;
      TARG.load(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }if((ok == 0) && (strcmp(argv[1], "-abs") == 0)){
      argc--;
      argv++;
      _abs=true;
      ok = 1;
    }
    else{cout << " option doesn't exist " << endl; exit(1);}
  }


  //  true_rescale(ORIG,RAD);
  //true_rescale(TRANSFORMED,RAD);
  DISTMAP=TRANSFORMED;
  ColumnVector distortion=getarealdistortion(ORIG,TRANSFORMED);
  double meandist=0;
  double val=0;
  for(int i=0;i<DISTMAP.nvertices();i++){
    if(_abs==true)val=abs(distortion(i+1));
    else  val=distortion(i+1);
    DISTMAP.set_pvalue(i,val);
    meandist+=val;
  }
  meandist/=DISTMAP.nvertices();
  cout << " mean=" << meandist << endl;
  DISTMAP.save(outdir + "distortion.func");

  if(_targ){ double MVD=Calculate_MVD(TRANSFORMED);  
  R.resample_scalar(DISTMAP,TARG,asin(MVD/(2*RAD)),EXCL);
  DISTMAP.save(outdir + "distortion_on_target.func");
  }
}
