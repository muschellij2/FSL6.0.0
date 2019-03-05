/*  msm.cc

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
#include "MeshReg/meshreg.h"
#include "msmOptions.h"
//#include "QPBO/QPBO.h"
//#include "HOCR/HOCR.h"

//#include "HOCR/QPBO.h"


#include "newran.h"
#include <time.h>

using namespace std;
using namespace NEWMESH;
using namespace NEWRAN;
using namespace MESHREG;

msmOptions* msmOptions::gopt = NULL;


 
int main(int argc, char *argv[]) {

  msmOptions& opts = msmOptions::getInstance();
  Log& logger = LogSingleton::getInstance();
  string CMpathin;

  opts.parse_command_line(argc,argv,logger);
  //opts.parse_command_line(argc,argv);
  MeshReg MR;

  //// INITIALISE VARIABLES //////////////////// 
  if(opts.printoptions.value()){ MR.print_config_options(); return 0; }
  if(opts.version.value()){ cout << " MSM_XS version 0.0" << endl; return 0; }
  if(opts.verbose.value())    MR.set_verbosity(opts.verbose.value());
  if(opts.debug.value())    MR.set_debug(opts.debug.value());

  MR.set_input(opts.inputmesh.value());
  if(opts.referencemesh.value()=="")  MR.set_reference(opts.inputmesh.value());
  else MR.set_reference(opts.referencemesh.value());

  ////// add anatomical meshes for stress and strain ////////////
  if(opts.inputanatmesh.value()!=""){
    if(opts.referenceanatmesh.value()=="") { cout << " Error: must supply both anatomical meshes or none "<< endl; exit(1);}
    MR.set_anatomical(opts.inputanatmesh.value(),opts.referenceanatmesh.value());
  }

  MR.set_outdir(opts.outbase.value());
  MR.set_output_format(opts.outformat.value());
 

  if(opts.transformed_sphere.set()) MR.set_transformed(opts.transformed_sphere.value());  
  if(opts.cfweight_in.set())  MR.set_input_cfweighting(opts.cfweight_in.value());
  if(opts.cfweight_ref.set()) MR.set_reference_cfweighting(opts.cfweight_ref.value());

  if(opts.in_register.set()){
    //// if data is supplied at a different resolution to the input or transformed mesh it is necessary to resample i.e. HCP 32K to native
    //// the in_register sphere HAS to be the sphere on which the data is supplied and the input or transformed mesh MUST be in alignment with it.
    newmesh in_register, target;
    resampler R; R.set_method("ADAP_BARY");
    boost::shared_ptr<BFMatrix > DATA;
    char filename[1000];
    if(opts.transformed_sphere.set()) target.load(opts.transformed_sphere.value());
    else target.load(opts.inputmesh.value());

    in_register.load(opts.in_register.value());
    set_data(opts.CMmatrixin.value(),DATA,in_register);
    R.resampledata(in_register,target,DATA,0.1);
    sprintf(filename,"%sinput_in_register.func.gii",opts.outbase.value().c_str());
    target.set_pvalues(DATA->AsMatrix());
    target.save(filename);
    CMpathin=filename;
   
  }else CMpathin=opts.CMmatrixin.value();
  cout <<  CMpathin << endl;
  MR.set_CMpathin(CMpathin);
  MR.set_CMpathref(opts.CMmatrixref.value());
  if(opts.training_data.set()){ MR.set_training(opts.training_data.value(), opts.trainingconcat.value()); 
    if(!opts.L1matlabpath.set()) { cout << "Must set matlab path " << endl; exit(1);}
    else{ MR.set_matlab(opts.L1matlabpath.value()); }
  }

  MR.run_multiresolutions(opts.multiresolutionlevels.value(),opts.smoothoutput.value(),opts.parameters.value());
 
  if(opts.in_register.set()){
    /// if data is supplied at a lower mesh resolution, then resample final warp accordingly. This is typically a HCP formating issue
    newmesh ORIG, FINAL, in_register;
    char filename[1000];

    if(opts.transformed_sphere.set()) ORIG.load(opts.transformed_sphere.value());
    else ORIG.load(opts.inputmesh.value());

    in_register.load(opts.in_register.value());
    sprintf(filename,"%ssphere.reg%s",opts.outbase.value().c_str(),MR.get_surf_format().c_str());
    cout << " filename " << filename << endl;
    FINAL.load(filename);
    barycentric_mesh_interpolation(in_register,ORIG,FINAL); 
    sprintf(filename,"%ssphere.in_register.reg.%s",opts.outbase.value().c_str(),MR.get_surf_format().c_str());
    in_register.save(filename);

  }
}
