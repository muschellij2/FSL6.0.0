/*  msmapplywarp.cc

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

#include "newmat.h"
#include "newmesh/meshfns.h"
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



using namespace NEWMAT;
using namespace NEWMESH;

void Usage()
{ cout << " msmapplywarp <to-be-transformedmesh> <outputname> -options  " << endl;
  cout << " Projects the to-be-transformed mesh though a transformation defined by an original mesh and its deformed counterpart" << endl;
  cout << " It is optional to supply the undeformed (original) mesh, but if no original mesh is supplied the algorithm will assume that the warp is prescribed by an icospheric template" << endl;
  cout << " -options " << endl;
  cout << " -original X.surf.gii " << endl;  
  cout << " -deformed DEFORMED.surf.gii (MUST be supplied in order to warp data)" << endl;
  cout << " -anat TARGET_SPHERE.surf.gii TARGET_ANAT.surf (2 inputs!). This will effectively project the INPUT anatomical mesh through the spherical warp." << endl;   
  cout << " -nospherical don't save spherical warp " << endl;  
  cout << " -affine (estimate affine transformation between input (-original) and deformed (-deformed) meshes and apply this to the to-be-transformed-mesh " << endl;  
  cout << " -readaffine X; where X  is an affine transformation matrix " << endl;  
  cout << " -writeaffine write out affine transformation matrix" << endl;  


}

void get_areas(newmesh &M){
  double val;
  ColumnVector meanarea(M.nvertices()); meanarea=0;

  for ( vector<NEWMESH::Triangle>::const_iterator i=M.tbegin() ; i!=M.tend(); i++){
    NEWMESH::Pt v1 = (*i).get_vertex_coord(0),  v2 = (*i).get_vertex_coord(1), v3 = (*i).get_vertex_coord(2);
	  
    val=computeArea(v1, v2, v3);
    meanarea((*i).get_vertex_no(0)+1)+=(1/3)*val; 
    meanarea((*i).get_vertex_no(1)+1)+=(1/3)*val; 
    meanarea((*i).get_vertex_no(2)+1)+=(1/3)*val; 
	
  }
  for (int i=0;i<M.nvertices();i++){
    M.set_pvalue(i,meanarea(i+1));
  }


}
int main(int argc, char **argv){

  
  newmesh to_be_deformed,to_be_deformed2, initial,final;
  newmesh TARGETSPHERE, TARGETANAT;
  newmesh ANAT_TRANS;
  Matrix _affinemat;
  string output;
  double MVD;
  boost::shared_ptr<RELATIONS > _rel;
  bool _isregular=true, _projectanat=false;
  bool _deformed=false, _areas=false;
  bool _affine=false, _supplyaffine=false, _writeaffine=false;
  char filename[1000];
  int ok;

  if(argc < 3){

    Usage();
    exit(0);
  }

 
  to_be_deformed.load(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;
  
  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-original") == 0)){
      argc--;
      argv++;    
      initial.load(argv[1]);
      _isregular=false;
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-deformed") == 0)){
      argc--;
      argv++;    
      final.load(argv[1]);
      _deformed=true;
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-anat") == 0)){
      argc--;
      argv++;   
      TARGETSPHERE.load(argv[1]); 
      argc--;
      argv++; 
      TARGETANAT.load(argv[1]);
      argc--;
      argv++; 
      _projectanat=true;
     
      ok = 1;
    } 
    else if((ok == 0) && (strcmp(argv[1], "-affine") == 0)){
      argc--;
      argv++;    
      _affine=true;   
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-outputareas") == 0)){
      argc--;
      argv++;    
      _areas=true;   
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-readaffine") == 0)){
      argc--;
      argv++;
      _supplyaffine=true;
      _affinemat=read_ascii_matrix(argv[1]);    
      argc--;
      argv++;
      ok = 1;
    }
     else if((ok == 0) && (strcmp(argv[1], "-writeaffine") == 0)){
      argc--;
      argv++;
      _writeaffine=true;
      ok = 1;
    }
    else{ cout <<argv[1] << "option npot present " << endl; exit(1);} 
  }
  
  if(!_deformed && !_supplyaffine) {
    if(!_projectanat){ cout <<_supplyaffine <<  "Cannot resample without a target and anatomical mesh. Please supply a target and anatomical mesh (-anat option) for resampling or a deformed mesh (-deformed option) for warping" << endl; exit(1);}
    ANAT_TRANS=projectmesh(to_be_deformed,TARGETSPHERE,TARGETANAT);
    sprintf(filename,"%s_anatresampled.surf.gii",output.c_str());
    ANAT_TRANS.save(filename);
    if(_areas){
      get_areas(ANAT_TRANS);
	
	
      sprintf(filename,"%s_anatresampled_areas.func.gii",output.c_str());
      ANAT_TRANS.save(filename);

    }

  }
  else{
 
	if(_supplyaffine){
		cout << _affinemat << endl;
		affine_transform(to_be_deformed,_affinemat);
        
	}else{
		newmesh ICO;

		if(_isregular){
		int ico=final.get_ico_resolution();
    
		// char buffer[1000];
		//SPH_in=SPH_orig;
		ICO.make_mesh_from_icosa(ico); true_rescale(ICO,RAD);
 
      //  double MVDinput=Calculate_MVD(ICO);  
      cout << " transform " << final.nvertices() << " ICO.nvertices " << ICO.nvertices() << " ico " << ico << endl;
    }else ICO=initial;

    if(_affine){
      cout << " affine" << endl;
      Matrix affinetrans=affine_transform(to_be_deformed,ICO,final);
      if(_writeaffine){
		sprintf(filename,"%s_affinewarp.txt",output.c_str());
		write_ascii_matrix(filename,affinetrans);
	}	
    }
    else
      barycentric_mesh_interpolation(to_be_deformed,ICO,final);
    
   
    if(_projectanat){
      ANAT_TRANS=projectmesh(to_be_deformed,TARGETSPHERE,TARGETANAT);
      sprintf(filename,"%s_projectedwarp.surf.gii",output.c_str());
      ANAT_TRANS.save(filename);

      if(_areas){
	get_areas(ANAT_TRANS);

	
	sprintf(filename,"%s_projectedwarp_areas.func.gii",output.c_str());
	ANAT_TRANS.save(filename);

      }

    }
	
  }
   sprintf(filename,"%s_warp.surf.gii",output.c_str());
   to_be_deformed.save(filename);
  }
}
