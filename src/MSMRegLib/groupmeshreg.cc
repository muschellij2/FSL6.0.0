/*  groupmeshreg.cc

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
#include "groupmeshreg.h"

#include <time.h>

namespace MESHREG {

 void GroupMeshReg::Initialize_level(int level)
  {
    cout << " In Initialize level " <<  DATAlist.size() <<endl;
    MeshModify::Initialize();
    vector<double>  sigma(MESHES.size(),_sigma_in[level]);
    NEWMESH::newmesh CONTROL; // original low res icosphere mesh
    CONTROL.make_mesh_from_icosa(_gridres[level] ); true_rescale(CONTROL,RAD);
    FEAT=boost::shared_ptr<featurespace>(new featurespace(DATAlist));
    FEAT->set_smoothing_parameters(sigma);

    FEAT->set_cutthreshold(_threshold); // will also generate exclusion masks at the same mesh resolution as datagrid
    FEAT->logtransform(_logtransform);// if true logtransforms AND normalises
    FEAT->varnorm(_varnorm);// variance normalises
    FEAT->intensitynormalize(_IN, _scale); // matches the intensities of the source to to the target (will rescale all to the top feature of the target if scale is true)
    FEAT->resamplingmethod(_dataInterpolator);
    FEAT->is_sparse(_issparse);
    SPH_orig=FEAT->Initialize(_genesis[level],MESHES,_exclude);  /// downsamples and smooths data, creates and exclusion mask if exclude is true
    cout << " after FEAT INIT " << endl;

    if(cost[level]=="AFFINE"){
      throw  MeshReg_error("GroupMeshREG ERROR:: affine not currently available");
    }else{

      cout << " 1 " << endl;
      bool multivariate=true;
      if(FEAT->get_dim()==1){multivariate=false; if (_simval[level]==4) {throw  MeshReg_error("MeshREG ERROR:: simval option 4 (alphaMI) is not suitable for univariate costfunctions");}}

      if(multivariate==true)  throw  MeshReg_error("GroupMeshREG ERROR:: multivariate not currently available");
     
      PARAMETERS.insert(parameterPair("multivariate",multivariate));  // check whether data is multivariate (and NOT using patch based alignment) and add this to the parameter set

      cout << " create MODEL " << _debug << endl;
      MODEL=boost::shared_ptr<SRegDiscreteModel>(new GroupDiscreteModel(PARAMETERS));
      cout << " after create Model " << endl;
      if(_debug) MODEL->set_debug();
      // if(_L1path=="" && _quartet==true)  throw  MeshReg_error("GroupMeshREG ERROR::quartet version requires matlab path");
      //else MODEL->set_L1path(_L1path);
      cout << " set feat " << endl;
      MODEL->set_featurespace(FEAT);
      cout << " set meshspace " << endl;

      MODEL->set_meshspace(SPH_orig,SPH_orig,MESHES.size());
      MODEL->Initialize(CONTROL);
    
      
    }

  }

  void GroupMeshReg::Evaluate()
  {

    cout << " in Evaluate " <<_level<< endl;
    newmesh OLDREG;

    for(int n=0;n<MESHES.size();n++){
     
      if(_level==1) ALL_SPH_REG.push_back(project_CPgrid(SPH_orig,OLDREG)); // first project data grid through any predefined transformation or, from transformation from previous resolution level
      else {
	OLDREG=ALL_SPH_REG[n];
	ALL_SPH_REG[n]=project_CPgrid(SPH_orig,OLDREG,n);
      }
    }
   
    cout << "  run discrete opt " <<  ALL_SPH_REG.size() << endl;
      // sprintf(filename,"SPHreg-EVal%d.surf", _level); SPH_reg.save(filename);
    run_discrete_opt(ALL_SPH_REG);
    
    if(_verbose) cout << "exit main algorithm     " << endl;
  
  }

  void GroupMeshReg::Transform(const string &filename){
    char fullpath[1000];
    // char buffer[1000];

    if(_verbose) cout << " Transform Group" << endl;
    
    for(int n=0;n<MESHES.size();n++){
      barycentric_mesh_interpolation(MESHES[n],SPH_orig,ALL_SPH_REG[n]); 
   
      sprintf(fullpath,"%ssphere-%d.reg%s",filename.c_str(), n,_surfformat.c_str());
      cout <<n << " " << SPH_orig.nvertices() << " " << ALL_SPH_REG[n].nvertices() << " fullpath " << fullpath <<  endl;

      MESHES[n].save(fullpath);
    }

  }

  void GroupMeshReg::saveTransformedData(const double &sigma, const string &filename){
    if(_verbose) cout << " save transformed data " << endl;

    resampler R; R.set_method("ADAP_BARY"); 
    boost::shared_ptr<RELATIONS > REL;
    double ang; ang=R.guess_angular_spacing(TEMPLATE.nvertices());
    char fullpath[1000];

 
    for(int n=0;n<MESHES.size();n++){
      sprintf(fullpath,"%stransformed_and_reprojected-%d%s",filename.c_str(),n,_dataformat.c_str());
	
      boost::shared_ptr<BFMatrix> DATA;
   

      REL = boost::shared_ptr<RELATIONS > ( new RELATIONS(MESHES[n],TEMPLATE,ang));
      REL->update_RELATIONS(MESHES[n]);
  
      set_data(DATAlist[n],DATA,MESHES[n]);
    
      R.resampledata(MESHES[n],TEMPLATE,DATA,0.0,REL);

      boost::shared_ptr<FullBFMatrix > pin =boost::dynamic_pointer_cast<FullBFMatrix>(DATA);
      TEMPLATE.set_pvalues(DATA->AsMatrix());
      TEMPLATE.save(fullpath);
    }
  }

  ////// iterates over discrete optimisations/////
    void GroupMeshReg::run_discrete_opt(vector<NEWMESH::newmesh> &source){
    resampler R;  R.set_method("ADAP_BARY");
    int iter=1;
    
    NEWMESH::newmesh transformed_controlgrid,targetmesh; 
    vector<newmesh> controlgrid;
    boost::shared_ptr<RELATIONS> controlgrid_neighbourhood;
    myparam::iterator it;
    int res,_itersforlevel;
    int numNodes;
    double energy=0,newenergy=0;
    char filename[1000];
    ofstream Energyout; 

    it=PARAMETERS.find("CPres");res=boost::get<int>(it->second);
    it=PARAMETERS.find("iters");_itersforlevel=boost::get<int>(it->second);
   
    numNodes=MODEL->getNumNodes();
    targetmesh=MODEL->get_TARGET();
      
  
    while(iter <= _itersforlevel){
      controlgrid.clear();
      cout << " reset meshspace " << endl;
      for(int n=0;n<source.size();n++){
	cout << n << " " << source[n].nvertices() << endl;
	MODEL->reset_meshspace(source[n],n); // source mesh is updated and control point grids are reset
	//	sprintf(filename,"sourcebeforereg-level%d-it%d-m%d.surf.gii",_level, iter,n); source[n].save(filename);
	controlgrid.push_back(MODEL->get_CPgrid(n)); 
	//sprintf(filename,"controlbeforereg-it%d-m%d.surf.gii",iter,n); controlgrid[n].save(filename);

      }
   
      cout << " MODEL->setupCostFunction(); " <<endl;
      MODEL->setupCostFunction();
     
      int *Labels=MODEL->getLabeling();
     
#ifdef HAS_HOCR

	Reduction HOCR_mode;

	if(_discreteOPT.compare(0,4,"HOCR")==0)
	  HOCR_mode=HOCR;
	else if(_discreteOPT=="ELC"){
	  HOCR_mode=ELC_HOCR;
	}else if(_discreteOPT=="ELC_approx"){
	  HOCR_mode=ELC_APPROX;
	}else {throw  MeshReg_error("discrete optimisation mode is not available");}
	
	//	tick_count start = tick_count::now();

	cout << _discreteOPT << " optimise 2 " << endl;
	srand (time(NULL));
	newenergy=Fusion::optimize(MODEL,HOCR_mode,_verbose);

#endif	//	MODEL->saveSTRAINS(iter);
    

	if(iter>1 && ((iter-1) % 2==0) && newenergy>=energy) { cout << iter << " level has converged" << (iter-1) %2  << endl; break;}

 
	for (int i = 0; i <numNodes; i++){
	  if(_verbose)
	    cout << i << " _iter " << iter <<" _iters " << _itersforlevel <<  " label " << Labels[i] <<  endl;
	}
   
	MODEL->applyLabeling();
	controlgrid_neighbourhood=MODEL->get_cp_neighbourhood();

	// apply these choices in order to deform the CP grid 
	for(int n=0;n<source.size();n++){
	  transformed_controlgrid=MODEL->get_CPgrid(n);     
	  cout << n << " transformed_controlgrid.nvertices() " << transformed_controlgrid.nvertices() << " controlgrid[n].nvertices() " << controlgrid[n].nvertices() << " source vertices " << " " << source[n].nvertices() << " " <<  controlgrid_neighbourhood->Ncols() << endl;
	  barycentric_mesh_interpolation(source[n],controlgrid[n],transformed_controlgrid,controlgrid_neighbourhood);
      
	  unfold(transformed_controlgrid);
	  MODEL->reset_CPgrid(transformed_controlgrid,n); // source mesh is updated and control point grids are reset   
	//}

      //sprintf(filename,"controlgridtrans-%d-%d.surf",iter,res);
      //transformed_controlgrid.save(filename);
      //sprintf(filename,"sourcetrans-%d-%d.surf",iter,res);
      //source.save(filename);

	  unfold(source[n]);
	}
	energy=newenergy;
	iter++;
    }

      

  }
}
   

 
