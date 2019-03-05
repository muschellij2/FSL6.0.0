/*  Copyright (C) 2009 University of Oxford  */

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
#include <cmath>
#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "newmat.h"
#include "newimage/newimageall.h"
#include "diffmodels.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;



string title="basgen - generate diffusion data using the ball and stick model";
string examples="baspred --dir=<bpx_directory> -o data.nii -b bvals -r bvecs -m mask.nii";

Option<bool> help(string("-h,--help"),false,
		       string("display this message"),
		       false,no_argument);
Option<string> dir(string("--dir"),string(""),
		   string("bedpostX directory"),
		   true,requires_argument);
Option<string> odata(string("-o,--output"),string(""),
		     string("output data"),
		     true,requires_argument);
Option<string> bvalsfile(string("-b,--bvals"),string(""),
			 string("ascii bvals file"),
			 true,requires_argument);
Option<string> bvecsfile(string("-r,--bvecs"),string(""),
			 string("ascii bvecs file"),
			 true,requires_argument);
Option<string> maskfile(string("-m,--mask"),string(""),
			string("brain mask"),
			true,requires_argument);



int do_bpxgen(){
  volume<float> s0vol,dvol,d_stdvol,f0vol;
  vector< volume<float> > fvol;
  vector< volume4D<float> > dyadsvol;
  volume<float> mask;

  volume4D<float> data;



  read_volume(mask,maskfile.value());
  read_volume(s0vol,dir.value()+"/mean_S0samples");
  read_volume(dvol,dir.value()+"/mean_dsamples");

  int model=1;
  if( fsl_imageexists(dir.value()+"/mean_d_stdsamples") ){
    read_volume(d_stdvol,dir.value()+"/mean_d_stdsamples");
    model=2;
  }
  bool usef0=false;
  if( fsl_imageexists(dir.value()+"/mean_f0samples") ){
    read_volume(f0vol,dir.value()+"/mean_f0samples");
    usef0=true;
  }  

  // fibres
  int nfibres=1;
  volume<float> tmpf;volume4D<float> tmpdyad;
  while(fsl_imageexists(dir.value()+"/mean_f"+num2str(nfibres)+"samples")){
    read_volume(tmpf,dir.value()+"/mean_f"+num2str(nfibres)+"samples");
    fvol.push_back(tmpf);
    read_volume4D(tmpdyad,dir.value()+"/dyads"+num2str(nfibres));
    dyadsvol.push_back(tmpdyad);
    nfibres++;
  }
  nfibres=(int)fvol.size();

  // bvals/bvecs  
  Matrix bvecs = read_ascii_matrix(bvecsfile.value());
  Matrix bvals = read_ascii_matrix(bvalsfile.value());
 
  data.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),bvals.Ncols());
  copybasicproperties(fvol[0],data);
  data=0;

  
  float s0,d,d_std=0;
  ColumnVector pvf(nfibres);
  Matrix dyads(nfibres,3);
  ColumnVector Y(bvals.Ncols());
  Y=0;

  cout << "generate data" << endl << endl;;
  for(int z=0;z<mask.zsize();z++){
    cout << "processing slice " << z << endl;
    for(int y=0;y<mask.ysize();y++){
      for(int x=0;x<mask.xsize();x++){
	if(mask(x,y,z)==0)continue;
	
	s0=s0vol(x,y,z);
	d=dvol(x,y,z);
	if(model==2)
	  d_std=d_stdvol(x,y,z);
	for(int f=1;f<=nfibres;f++){
	  pvf(f)=fvol[f-1](x,y,z);
	  dyads.Row(f) << dyadsvol[f-1](x,y,z,0)
		       << dyadsvol[f-1](x,y,z,1)
		       << dyadsvol[f-1](x,y,z,2);
	}

	if(model==1){
	  PVM_single pvm1(Y,bvecs,bvals,nfibres,usef0);
	  pvm1.set_s0(s0);
	  pvm1.set_d(d);
	  pvm1.set_f(pvf);
	  pvm1.set_th_ph(dyads);
	  if(usef0){
	    pvm1.set_f0(f0vol(x,y,z));
	  }

	  Y = pvm1.get_prediction();
	 
	  for(int t=1;t<=bvals.Ncols();t++)
	    data(x,y,z,t-1)=Y(t);
	}
	else{
	  PVM_multi pvm2(Y,bvecs,bvals,nfibres,usef0);
	  pvm2.set_s0(s0);
	  pvm2.set_d(d);
	  pvm2.set_d_std(d_std);
	  pvm2.set_f(pvf);
	  pvm2.set_th_ph(dyads);
	  if(usef0){
	    pvm2.set_f0(f0vol(x,y,z));
	  }

	  Y = pvm2.get_prediction();
	  
	  for(int t=1;t<=bvals.Ncols();t++)
	    data(x,y,z,t-1)=Y(t);
	}
      }
    }
  }

  //cout<<"saving results" << endl;
  data.setDisplayMaximumMinimum(data.max(),0);
  save_volume4D(data,odata.value());
  
  cout<<endl<<"Done."<<endl;

  return 0;
}

int main(int argc,char *argv[]){

  Tracer tr("main");
  OptionParser options(title,examples);

  try{
    options.add(help);
    options.add(dir);
    options.add(odata);
    options.add(bvalsfile);
    options.add(bvecsfile);
    options.add(maskfile);

    options.parse_command_line(argc,argv);

    
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
      options.usage();
      exit(EXIT_FAILURE);
    }
  }
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 
  
  return do_bpxgen();
  
  
}
