/*  ContinuosOpt.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  

    Some sections of code inspired by Alek Petrovic.*/
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
/// Class for gradient based optimisation - hopefully soon to be obselete just need to convert affine method to discrete ////

#ifndef ContinuosOpt_h
#define ContinuosOpt_h

#include "newmesh/featurespace.h"

#include <FastPD/FastPD.h>

#ifdef HAS_HOCR
#include <DiscreteOpt/Fusion.h>
#endif

using namespace DISCRETEOPT;


#define CORLIM  1E-10  // when calculating corrzero consider any correlations below this value to be ZERO- NOT USED ANY MORE?

using namespace NEWMESH;

namespace MESHREG{
  
  class MESHREGException: public std::exception
    {
    private:
      std::string m_msg;
    public:
      MESHREGException(const std::string& msg) throw(): m_msg(msg) {}
      
      virtual const char * what() const throw() {
	return string("Fnirt: msg=" + m_msg).c_str();
      }
      
      ~MESHREGException() throw() {}
    };


  class MeshCF {    /// basic mesh cost function - uses the gradient of the similarity function, estimated using weighted least squares
   

  protected:
    NEWMESH::newmesh _TARGET; // TARGET MESH
    NEWMESH::newmesh _SOURCE; // SOURCE MESH
   
    boost::shared_ptr<RELATIONS> _rel; // saves mesh neighbourhood information

    Matrix _inweight; // exclusion/weighting mask for cost function masking
    Matrix _refweight;

    const boost::shared_ptr<featurespace> FEAT; /// holds data
    sparsesimkernel<double> sim; // similarity matrix
  
  
    double MVD; // mean vertex distance

    ///////// user defined parameters /////////
    int _simmeasure; 
    int _iters; // total iterations
    float _stepsize;
    float _spacing;
       //////////////////////////////////////////

    
    string m_outdir;
  
    vector<Tangs> BASIS; /////////////  tangent vector bases for all vertices

    bool  _verbosity;
    ////// GRAD DESCENT PARAMETERS /////////////////
    mutable ColumnVector current_sim;

  
    double min_sigma;
    double CF; 
    double totJP;
    double totJD;

   
  public:
   
     MeshCF(const NEWMESH::newmesh &        target, 
	   const NEWMESH::newmesh &        source,
	   const Matrix T_cfweight, 
	   const Matrix S_cfweight,
	   const boost::shared_ptr<featurespace> features)      
      : _TARGET(target),_SOURCE(source),_inweight(S_cfweight),_refweight(T_cfweight), FEAT(features)
    {
    
      ////  DEFAULT PARAMETRISATION /////
      _simmeasure=3;
      _iters=10000;
      _verbosity=false;
    }

    MeshCF (const NEWMESH::newmesh target, const NEWMESH::newmesh &source,
	    const boost::shared_ptr<featurespace> features):_TARGET(target), _SOURCE(source),FEAT(features){  
      _simmeasure=3;
      _iters=20;
      _verbosity=false;
      _inweight.ReSize(1, _SOURCE.nvertices()); _refweight.ReSize(1, _TARGET.nvertices());  // if no cf weighting then set weighting to 1 by default
      _inweight=1; _refweight=1;
      
    } 
    
    void set_parameters(myparam PAR){
      myparam::iterator it;
      it=PAR.find("iters");_iters=boost::get<int>(it->second);
      it=PAR.find("simmeasure");_simmeasure=boost::get<int>(it->second); 
      it=PAR.find("verbosity");_verbosity=boost::get<bool>(it->second);
      it=PAR.find("stepsize");_stepsize=boost::get<float>(it->second);
      it=PAR.find("gradsampling");_spacing=boost::get<float>(it->second);
      it=PAR.find("outdir");m_outdir=boost::get<string>(it->second);
  
    }

    ///////////////////// INITIALIZE ///////////////////
    virtual void Initialize();  
    void set_simmeasure(const int & simval){_simmeasure=simval;}
     
    //////////////// MAKE UPDATES ////////////////////
    
    void update_similarity();  
    void update_similarity(const int &, vector<int> &); 
    void update_source(const NEWMESH::newmesh& M) {_SOURCE=M; };
    
    //////////////// SIMILARITY GRADIENT ESTIMATION /////////////////////  
    ///// similarity gradient via weighted regression -
    ColumnVector WLS_simgradient(const Tangs &T, int ,const vector<int> &);  
    //// prepares data for sim gradient calculation
   
    ColumnVector Evaluate_SIMGradient(int i,const Tangs &T);
    
    virtual  NEWMESH::newmesh run()=0;
    ////////////////// RETURNING FUNCTIONS ///////////////////////

    NEWMESH::newmesh get_REG()const{return _SOURCE;}  
  };


  class affineMeshCF: public MeshCF 
  {


  public:
    
  affineMeshCF(const NEWMESH::newmesh & target, const NEWMESH::newmesh &source,const boost::shared_ptr<featurespace> features ): MeshCF(target,source,features){}
    virtual ~ affineMeshCF() {};
    
    void Initialize(){MeshCF::Initialize();}
    /////////////// AFFINE FUNCTIONS - uses image gradients and euler rotations - taken from Alek Petrovic's implementation.
    void Rotate_IN_mesh(const double &, const double &,const double &);
    double Affine_cost_mesh(const double & , const double & ,const double &);
    NEWMESH::newmesh run();
    ColumnVector return_transformation(){return REC;} 
   private:

    ColumnVector REC;
    affineMeshCF();
    virtual  affineMeshCF& operator=(const  affineMeshCF& inf) {throw MESHREGException("LM_MeshCF:: Assignment explicitly disallowed"); return(*this);}

  };
  

  
}

#endif
