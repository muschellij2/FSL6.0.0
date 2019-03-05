/*  meshmodify.h

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
/* this class sets up the framework for mesh registration */
#ifndef meshmodify_h
#define meshmodify_h

#include <time.h>
#include "ContinuosOpt.h"


using namespace std;
using namespace NEWMAT;
using namespace NEWMESH;
using namespace MISCMATHS;
using namespace FPD;
using namespace DISCRETEOPT;

namespace MESHREG {

  class  MeshReg_error: public std::exception
    {
    public:
      MeshReg_error(const string& pmsg) throw() : msg(pmsg) {}
      const char *what() const throw() {return(string("MESHREG::" + msg).c_str());}
      ~MeshReg_error() throw() {}
    private:
      string msg;
    };
  class MeshModify{
    
  protected:
    
    Tangent Tang;  // holds tangent plane basis at each point
    vector<int> datarange; // for processing only a range of data points - primarily only used for segmentation, likely temporary
    
    // use smart pointers allows for checking presence of these data objects before code is run
    string CMfile_in; //location of connectivity/data matrix for input
    string CMfile_ref; //location of connectivity/data matrix for reference
    string _outdir;
    string _surfformat;
    string _dataformat;
    
    vector<NEWMESH::newmesh> MESHES;  // original input (moving) mesh
    NEWMESH::newmesh TEMPLATE;
    NEWMESH::newmesh transformed_mesh;  // original input (moving) mesh in transformed position (i.e. from previous alignment step)
    NEWMESH::newmesh in_anat;
    NEWMESH::newmesh ref_anat; 
    NEWMESH::newmesh SPH_orig; // original low res icosphere mesh
    NEWMESH::newmesh SPH_reg;  // transformed low res icosphere mesh
    NEWMESH::newmesh ANAT_orig; // original low res icosphere mesh
    // NEWMESH::newmesh ANAT_ref;  // transformed low res icosphere mesh
    
    boost::shared_ptr<NEWMESH::newmesh> IN_CFWEIGHTING; // cost function weights for high resolution meshes
    boost::shared_ptr<NEWMESH::newmesh> REF_CFWEIGHTING;
    
    Matrix SPHin_CFWEIGHTING; // downsampled costfunction weightings
    Matrix SPHref_CFWEIGHTING;
    
    boost::shared_ptr<featurespace> FEAT;
    
    vector<string> DATAlist;
    //////////////DATA OPTIONS //////////////////
    
    bool _logtransform; // logtransform and rescale the data 
    bool _IN;
    bool _issparse;
    string _meshInterpolator; // TPS by default should be barycentric?
    string _dataInterpolator; // ADAP_BARY by default was Gaussian
    
    string _discreteOPT;
    string _L1path;
    //// REGISTRATION PARAMETERS /////////////
    myparam PARAMETERS;
    vector<string> cost;   // controls registration method i.e. affine, discrete, gradient descent
    vector<int>   _genesis;           // ico mesh resolution at this level
    
    vector<float> _sigma_in;        // smoothing of input 
    vector<float> _sigma_ref;      // smoothing of reference 
    vector<int> _simval; // code determines how similarity is assessed 1 is aleks' correlation measure 2 is conventional correlation 3=SSD 4=NMI 5 alpha entropy
    
    vector<float> _lambda;         // controls regularisation
    vector<float> _threshold;         // controls cut exclusion (2D upper and lower thresholds for defining cut vertices)
    vector<int> _iters; // total per resolution level
    vector<int> _gridres; // control point grid resolution (for discrete reg)
    vector<int> _anatres; // control point grid resolution (for discrete reg)
    
    vector<int> _sampres; // sampling grid for discrete reg (should be higher than control point)
    vector<int>  _alpha_kNN;   // for alpha mutual information
    bool _scale; // rescale all features to have distribution of the first feature in a multivariate set
    
    
    double     MVD;         //mean inter-vertex distance
    bool _verbose;
    bool _exclude;  // exclusion zone
    bool _varnorm; // variance normalise
    bool _initialise;
    bool _tricliquelikeihood;
    bool _debug;
    bool _concattraining;
    bool _usetraining;
    bool _anat;
    bool _weight;
    bool _quartet;
    bool _set_group_lambda;

    float _potts;
    int _resolutionlevels;  /// default 1 
    int _regmode; //regulariser option
    int _numthreads;
    /////////// REGULARISER OPTIONS //////////////
       
    float _regexp; // choice of exponent
    float _regscaling; // choice of exponent scaling for options 2 and 3
    float _pairwiselambda; ///scaling for group alignment
    float _maxdist; // max areal disortion allowed before scaling is applied
    float _shearmod; // for strain regulariser
    float _bulkmod;    /// for strain regulariser
    float _cprange;
    int _featlength;
   
    //////////////// AFFINE PARAMETERS //////////////////////
    float _affinestepsize;
    float _affinegradsampling;

   
    
  public:
    
    // Constructors
    
    inline MeshModify(){
      MESHES.resize(2,newmesh());
      _verbose=false;
      _resolutionlevels=0; 
      _issparse=false;
      _debug=false;
      _set_group_lambda=false;
      _usetraining=false;
      _tricliquelikeihood=false;
      _quartet=false;
      _anat=false;
      _concattraining=false;
      _potts=0.0;
      _surfformat=".surf";
      _dataformat=".func";
      _discreteOPT="FastPD";
      _L1path="";
      _cprange=1;
      _bulkmod=10;
      _shearmod=2.5;
      _numthreads=1;
      FEAT=boost::shared_ptr<featurespace>(new featurespace());
    };
    
    // Destructor
    virtual inline ~MeshModify(){};
    
    //// parses config file
    void parse_reg_options(const string &);
    
    ////// sets up parameter space for one resolution level of the registration
    void fix_parameters_for_level(const int & );

    inline NEWMESH::Pt setpoint( double X, double Y, double Z){Pt p; p.X=X;p.Y=Y;p.Z=Z; return p;}
    
    // Initialize 
    /// reads high resolution input mesh   
    inline void set_input(const NEWMESH::newmesh &M) {MESHES[0]=M; recentre(MESHES[0]);true_rescale(MESHES[0],RAD); };
    //// reads high resolution target mesh
    inline void set_reference(const NEWMESH::newmesh &M) {MESHES[1]= M; recentre(MESHES[1]);  true_rescale(MESHES[1],RAD); };

    /// reads high resolution input mesh   from path
    inline void set_input(const string &M) {MESHES[0].load(M); recentre(MESHES[0]);    true_rescale(MESHES[0],RAD); };

    inline void set_inputs(string s){  
      vector<string> meshlist=read_ascii_list(s); 
      newmesh tmp;
      MESHES.clear();
      for(unsigned int i=0;i<meshlist.size();i++) {
	if(_verbose)	cout << i << " " << meshlist[i] << endl;
	tmp.load(meshlist[i]);
	MESHES.push_back(tmp);

      }

    }

    /// reads high resolution target mesh  from path
    inline void set_reference(const string &M) {MESHES[1].load(M);  recentre(MESHES[1]);  true_rescale(MESHES[1],RAD);  };

    inline void set_template(const string &M) {TEMPLATE.load(M); recentre(TEMPLATE);  true_rescale(TEMPLATE,RAD);  };


    inline void set_anatomical(const string &M1,const string &M2) {_anat=true; in_anat.load(M1); ref_anat.load(M2);};

    /// transformed mesh allows registration to be initialised using a  transformation from a previous iteration
    inline void set_transformed(const string &M) {transformed_mesh.load(M);true_rescale(MESHES[1],RAD); _initialise=true; };
    //// reads costfunction weighting mask for input data    
    inline void set_input_cfweighting(string E) { IN_CFWEIGHTING= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(MESHES[0])); IN_CFWEIGHTING->load(E,false,false);  true_rescale(*IN_CFWEIGHTING,RAD); };
    /// reads costfunction weighting mask for target data
    inline void set_reference_cfweighting(string E) {REF_CFWEIGHTING= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(MESHES[1]));   REF_CFWEIGHTING->load(E,false,false); true_rescale(*REF_CFWEIGHTING,RAD);};
    /// sets output path   
    inline void set_outdir(string s) {_outdir=s;};
    /// sets output format   
    void set_output_format(string type);
   
    //// leads to additional outputs being saved out
    inline void set_debug(bool test) {_debug=test;};
    //////// sets path to input data
    inline void set_CMpathin(string s){CMfile_in = s;}; 
    /// sets path to reference data
    inline void set_CMpathref(string s){CMfile_ref=s;};
    /// gets all paths to training data
    inline void set_training(string s, bool concat){DATAlist=read_ascii_list(s);  _concattraining=concat; _usetraining=true;};

    inline void set_data_list(string s){DATAlist=read_ascii_list(s); };


    inline void set_matlab(string s){_L1path=s;}
    
    inline void set_datarange(const vector<int> range){datarange=range;};

    void Initialize();
    virtual void Initialize_level(int)=0; // for multires registration
     
    void print_config_options(){parse_reg_options("usage");}
    void set_verbosity(bool V){_verbose=V;}
    void is_sparse(bool sp){_issparse=sp;}; /// input data is sparse
    void check(); // checks you have all the necessary data
    
    ////////COMMON FUNCTIONS//////////////// 
  
    virtual void Evaluate()=0; 
    virtual void Transform(const string &)=0;

    ////// UPDATING /////////
  
    void update_similarity();
  
    ////// RETURNING FUNCTIONS /////////

    inline  NEWMESH::newmesh return_registered_input_mesh()const {return MESHES[0];};

    inline  string get_indata_path()const {return CMfile_in;};
    inline  string get_surf_format()const {return _surfformat;};
    inline  string get_refdata_path()const {return CMfile_ref;};
    /// save low resolution estimate of registration
    virtual inline void saveSPH_reg(const string &filename) const {
      char fullpath[1000];
      cout << " orig save SPH " << endl;
      sprintf(fullpath,"%ssphere.LR.reg%s",filename.c_str(),_surfformat.c_str());
      SPH_reg.save(fullpath);

    

    }
    /// saves transformed and reprojected data
    virtual void saveTransformedData(const double &,const string &filename)=0;
    
    Matrix combine_costfunction_weighting(const Matrix &, const Matrix &);
  };


}
 




#endif

