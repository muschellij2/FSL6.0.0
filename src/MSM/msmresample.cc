/*  msmresample.cc

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
#include "newmat.h"
#include "newmesh/meshfns.h"
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



using namespace NEWMAT;
using namespace NEWMESH;

void Usage()
{ cout << " msmresample <in_anat> <in_sphere> <ico resolution> <output base> -options " << endl;
  cout << " -data X supply data for resampling " << endl;
  cout << " -adap resample adaptively " << endl;

}


int main(int argc, char **argv){

  
  newmesh in_anat,in_sphere,ico;
  newmesh ANAT_res;
  string output;
  double res;
  char filename[1000];
  resampler R; R.set_method("ADAP_BARY");
  int ok;
  bool _resampledata=false;
  bool _adap=false;
  if(argc < 3){

    Usage();
    exit(0);
  }

 
  in_anat.load(argv[1]);
  argc--; 
  argv++;
  in_sphere.load(argv[1]);
  argc--; 
  argv++;
  res=atoi(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;
  
  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-data") == 0)){
      argc--;
      argv++;    
      _resampledata=true;
      in_sphere.load(argv[1],false,false);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-adap") == 0)){
      argc--;
      argv++;    
      _adap=true;
      ok = 1;
    }else{cout << " option doesn't exist " << endl; exit(1);}
  }
  cout << " Make ico " << endl;
  ico.make_mesh_from_icosa(res); true_rescale(ico,RAD); 
  
 
  Matrix TRANSLATE=recentre(in_sphere);
  true_rescale(in_sphere,RAD);
  //recentre anat ///
Pt mean;
  if(TRANSLATE.NormFrobenius() > EPSILON){
  
 // ColumnVector P_in(4);
  
  cout << " translate anat " << endl;
  for (int i=0;i< in_anat.nvertices();i++)
      mean+=in_anat.get_coord(i);

    mean/=in_anat.nvertices();
    cout << "inant  before mean" << mean.X << " mean.Y " << mean.Y << " mean.Z " << mean.Z <<  endl;
    TRANSLATE=recentre(in_anat);
    mean*=0;
  /*for ( vector<boost::shared_ptr<NEWMESH::Mpoint> >::const_iterator i= in_anat.vbegin(); i!=in_anat.vend(); i++){
       NEWMESH::Pt p = (*i)->get_coord();
       NEWMESH::Pt p2;
       if(p.norm()){
	 P_in(1) = p.X; P_in(2) = p.Y; P_in(3) = p.Z; P_in(4) = 1;  
	 P_in = TRANSLATE * P_in;
	 p2.X=P_in(1); p2.Y=P_in(2);  p2.Z=P_in(3);
	 (*i)->set_coord(p2);
       }
      
     }	
     * */ 
for (int i=0;i< in_anat.nvertices();i++)
      mean+=in_anat.get_coord(i);

    mean/=in_anat.nvertices();
    cout << "inant  mean" << mean.X << " mean.Y " << mean.Y << " mean.Z " << mean.Z <<  endl;
  
  }
 
  ANAT_res=mesh_resample(in_anat,in_sphere,ico,_adap);
  sprintf(filename,"%s-anat.surf",output.c_str());
  mean*=0;
  for (int i=0;i< ANAT_res.nvertices();i++)
      mean+=ANAT_res.get_coord(i);

    mean/=ANAT_res.nvertices();
    cout << "anat res  mean" << mean.X << " mean.Y " << mean.Y << " mean.Z " << mean.Z <<  endl;
  ANAT_res.save(filename);
  sprintf(filename,"%s-regular_sphere.surf",output.c_str());
  ico.save(filename);

  if(_resampledata){
    R.resample_scalar(in_sphere,ico,1);
    sprintf(filename,"%s-resampled_data.func",output.c_str());
    in_sphere.save(filename);

  }
}
