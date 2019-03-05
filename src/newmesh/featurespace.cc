/*  featurespace.cc

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
#include "featurespace.h"

namespace NEWMESH {


  ////////////////////////// RESAMPLING OF DATA DIMENSION//////////////////////////////

  NEWMESH::newmesh featurespace::Initialize(const int &ico, vector<NEWMESH::newmesh> &IN,const bool & exclude){
   
    NEWMESH::newmesh icotmp;
    bool isfunc;

    if(IN.size()!=CMfile_in.size()){ throw  NEWMESHException(" NEWMESH::featurespace::Initialize do not have the same number of datasets and surface meshes");	}
    else {DATA.resize(IN.size(),boost::shared_ptr<BFMatrix > ()); }
    
    if(ico>0){icotmp.make_mesh_from_icosa(ico); true_rescale(icotmp,RAD); }
     

    
    for (unsigned int i=0;i<IN.size();i++){
      boost::shared_ptr<BFMatrix> tmp=boost::shared_ptr<BFMatrix > ();
      isfunc=set_data(CMfile_in[i],DATA[i],IN[i],_issparse); 
      IN[i].set_pvalues(DATA[i]->AsMatrix());
    
      if(ico==0){icotmp=IN[i];} 
      // create exclusion mask 

 
      if(exclude){
	NEWMESH::newmesh excl_tmp=create_exclusion(IN[i],DATA[i]->AsMatrix(),_fthreshold[0],_fthreshold[1]);
	EXCL.push_back(boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(excl_tmp))); ///mask data according to some min (_fthreshold[0]) and max (_fthreshold[1]) threshold
		     
      }else{EXCL.push_back(boost::shared_ptr<NEWMESH::newmesh>()); }

      ///// downsample data to regular grid of resolution defined by "ico"
    
      resample(_sigma_in[i],DATA[i],icotmp,IN[i],EXCL[i]);
      
      

    }
    icotmp.set_pvalues(DATA[0]->AsMatrix());
    
    /// intensity normalise using histogram matching 
    if(_intensitynorm){ 
     
      for (unsigned int i=1;i<IN.size();i++){
	multivariate_histogram_normalization(*DATA[i],*DATA[0],EXCL[i],EXCL[0],_scale);  // match input data feature distributions to equivalent in ref, rescale all to first feature in reference if _scale is 

      }
   
    }
    
    
    
    if(_logtransform){
      for (unsigned int i=0;i<IN.size();i++)
	log_transform_and_normalise(*DATA[i]);
    }

    if(_varnorm){
      for (unsigned int i=0;i<IN.size();i++){
	varnorm(DATA[i],EXCL[i]); 
       
      
      }
    }
   
    return icotmp;
 
  }

  ////////////////////////// RESAMPLING OF DATA DIMENSION//////////////////////////////


  void featurespace::resample(const double &sigma, boost::shared_ptr<BFMatrix> &DATAtmp, NEWMESH::newmesh &icotemp,const NEWMESH::newmesh &M, boost::shared_ptr<NEWMESH::newmesh> &EXCLtmp){
    
    resampler R;
    newmesh tmp=M;

    if(_resamplingmethod=="ADAP_BARY"){
      R.set_method("ADAP_BARY");  
      R.resampledata(tmp,icotemp,EXCLtmp,DATAtmp,sigma);

 
      if(sigma>0){
	R.set_method("GAUSSIAN");
	R.smooth_data(sigma,DATAtmp,icotemp);
	
      }   
    }
    else{
      R.set_method("GAUSSIAN");     
      R.resampledata(tmp,icotemp,EXCLtmp,DATAtmp,sigma);
    }
  }


  void featurespace::smooth(NEWMESH::newmesh &IN, boost::shared_ptr<BFMatrix> &tmp, const double &_sigma)
  {
    
    resampler R;
    R.set_method("GAUSSIAN");
    R.smooth_data(_sigma,tmp,IN);
   
  }

  void featurespace::varnorm(boost::shared_ptr<BFMatrix> & DATA, boost::shared_ptr<NEWMESH::newmesh> &EXCL){ /// need to condense code of next two functions into one
   

    vector<vector<double> > _data;
    vector<double> var1;
    vector<double> tmp;
    int ind;
    ind=0;
    
    
    for (unsigned int i=1;i<=DATA->Ncols();i++){
      if(!EXCL.get() || EXCL->get_pvalue(i-1)>0){
	_data.push_back(tmp);
	for(unsigned int k=1;k<=DATA->Nrows();k++){
	  _data[ind].push_back(DATA->Peek(k,i));
	}
	ind++;
      }
    } 
    ind=0;
  
    var1=online_variance_normalize(_data); /// can replace with normalise(_sourcedata)
   
    ind=0;

    for (unsigned int i=1;i<=DATA->Ncols();i++){
      if(!EXCL.get() ||EXCL->get_pvalue(i-1)>0){
	for(unsigned int k=1;k<=DATA->Nrows();k++){
	  DATA->Set(k,i,_data[ind][k-1])  ;  
	 
	}
	ind++;
      }
    }
       
  }


  /// calculate mean and variance simultaneously to speed up estimation
  vector<double> featurespace::online_variance_normalize( vector<vector<double> > &M){   // variance normalisation is susceptible to precision error especially when data points are close to the mean
  
    RowVector mean(M[0].size());
    vector<double> var(M[0].size(),0);
    double delta,maxvar;
    mean=0; maxvar=0;

    for (unsigned int j=0;j<M[0].size();j++){
      var[j]=0; mean(j+1)=0;
      for (unsigned int i=0;i<M.size();i++){
	delta=M[i][j]-mean(j+1);
	mean(j+1)=mean(j+1)+delta/(i+1);
	var[j]= var[j] +delta*(M[i][j]-mean(j+1));
      }
      var[j]= var[j]/(M.size()-1);

      for (unsigned int i=0;i<M.size();i++){
	M[i][j]-=mean(j+1);
	if(var[j]>0) M[i][j]=M[i][j]/sqrt(var[j]);
	if(var[j]>maxvar) maxvar=var[j];
      }

     }

    for (unsigned int i=0;i<M[0].size();i++){
	var[i]/=maxvar;
      }

    return var;
   
  }


 
    
}

  
  
