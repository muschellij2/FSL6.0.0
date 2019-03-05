/*  featurespace.h

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
/*  CLASS FOR PROCESSING DATA PRIOR TO MSM REGISTRATION */
/* ideally this class should probably be split up and the major components should go into newmesh for manipulation the mesh data for each mesh individually */

#if !defined(featurespace_h)
#define featurespace_h

#include <fstream>
#include <stdio.h>

#include "meshfns.h"
#include "utils/options.h"

using namespace Utilities;

namespace NEWMESH{

   enum StringValue { evNotDefined, 
		      evStringValue1, 
		      evStringValue2, 
		      evStringValue3, 
		      evEnd };
  // Map to associate the strings with the enum values
   static std::map<std::string, StringValue> s_mapStringValues;
  
   /* features - should this also inherit from newmat?  */
  class featurespace{
    
  private:
   
    NEWMESH::newmesh source;

    vector<boost::shared_ptr<BFMatrix> > DATA; // holds generic BFMATRIX data which can be sparse or full matrices
    vector<boost::shared_ptr<NEWMESH::newmesh> >  EXCL;  // exclusion masks for binary weighting of the data during resampling
   
    vector<string >  CMfile_in;  // path to data
    vector<double> _sigma_in;  // smoothing parameters for input and reference
    vector<float> _fthreshold;

    string inorig;   
    string reforig;
    string _resamplingmethod;
 
    bool _logtransform;  // will log transform and normalise
    bool _intensitynorm; // will histogram match
    bool _scale;  /// will rescale each feature to crudley match the distribution of the first in a multivariate distribution
    bool _issparse;  /// notes that data is sparse
    bool _varnorm;  // performs online variance normalisation - maybe replace with non online version called during logtransformandnormalise()

  public:

    featurespace(){};
    ~featurespace(){};
    featurespace(const string &datain, const string &dataref){
      _sigma_in.resize(2,5.0);_logtransform=false; _issparse=false;  _scale=false; _intensitynorm=false; _fthreshold.resize(2,0.0); 
      CMfile_in.push_back(datain); CMfile_in.push_back(dataref);
    };

    featurespace(const string &datain, const vector<string> &datareflist){
      _sigma_in.push_back(5);_logtransform=false; _issparse=false;  _scale=false; _intensitynorm=false; _fthreshold.resize(2,0.0); 
       CMfile_in.push_back(datain); 
       for(int i=0;i<(int) datareflist.size();i++){
	 CMfile_in.push_back(datareflist[i]);
	 _sigma_in.push_back(5);
       }

    };

    featurespace(const vector<string> &datalist){
      _logtransform=false; _issparse=false;  _scale=false; _intensitynorm=false; _fthreshold.resize(2,0.0); 
      for(int i=0;i<(int) datalist.size();i++){
	 CMfile_in.push_back(datalist[i]);
	 _sigma_in.push_back(5);
       }

    };
    /////////////////// INITIALIZE //////////////////////////////
    void set_smoothing_parameters(const vector<double> s){ 
      _sigma_in.clear();
      if(s.size()!=CMfile_in.size()){ 
	if(s.size()==1){ 
	  for (int i=0;i<(int) CMfile_in.size();i++)
	    _sigma_in.push_back(s[0]);	  
	}else throw  NEWMESHException("Mewmesh::featurespace:: smoothing sigma size incompatible with data dimensions");
      }else _sigma_in=s;};

    void set_cutthreshold(vector<float> & thr){_fthreshold=thr;};
    void logtransform(const bool &log){_logtransform=log;}
    void varnorm(const bool &norm){_varnorm=norm;}

    void is_sparse(const bool &sp){_issparse=sp;}
    void intensitynormalize(const bool & norm, const bool &scale){_intensitynorm=norm;};
    void resamplingmethod(string method){_resamplingmethod=method;}

    NEWMESH::newmesh Initialize(const int &, vector<NEWMESH::newmesh> &,const bool &);
    /////////////////// RESAMPLE DATA /////////////////////////
    void resample(const double &, boost::shared_ptr<BFMatrix> &,NEWMESH::newmesh &,const NEWMESH::newmesh &, boost::shared_ptr<NEWMESH::newmesh> &);
    void smooth(NEWMESH::newmesh &, boost::shared_ptr<BFMatrix> &, const double &);
    void varnorm(boost::shared_ptr<BFMatrix> &, boost::shared_ptr<NEWMESH::newmesh> &); // combine this and next function
    vector<double> online_variance_normalize(vector<vector<double> >  &); // mean centre and v ariance normalize

    //////////////////// ACCESS //////////////////////////////////////////////
    string get_path(const int i)const{return CMfile_in[i];}
    int get_dim()const{return DATA[0]->Nrows();}

    /////// PAIRWISE ACCESSS /////////////////////////////////////////////////////////////////////
    string get_input_path()const{return CMfile_in[0];}
    string get_reference_path()const{return CMfile_in[1];}
    double get_input_val(const int &i,const int &j)const{ return DATA[0]->Peek(i,j);}
    double get_ref_val(const int &i,const int &j)const{ return DATA[1]->Peek(i,j);}
    double get_data_val(const int &i,const int &j,const int n)const{ return DATA[n]->Peek(i,j);}
    void set_input_val(const int &i,const int &j, const double &val)const{ return DATA[0]->Set(i,j,val);}
    void set_ref_val(const int &i,const int &j, const double &val)const{ return DATA[1]->Set(i,j,val);}
    boost::shared_ptr<NEWMESH::newmesh>  get_input_excl()const{return EXCL[0];}
    boost::shared_ptr<NEWMESH::newmesh>  get_reference_excl()const{return EXCL[1];}

    boost::shared_ptr<BFMatrix>  get_input_data()const{
      return DATA[0];
    }

    boost::shared_ptr<BFMatrix>  get_data(const int i)const{
      return DATA[i];
    }
    Matrix  get_data_matrix(const int i)const{
       return DATA[i]->AsMatrix();
    }
    boost::shared_ptr<BFMatrix> get_reference_data()const{
      return DATA[1];
    }
  };
  
}
#endif
  
  
