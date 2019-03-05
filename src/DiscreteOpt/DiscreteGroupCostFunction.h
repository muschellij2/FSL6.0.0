/*  DiscreteGroupCostFunction.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2014 University of Oxford  */

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
#include <vector>
#include "DiscreteCostFunction.h"
//#include <boost/thread.hpp>
#include <boost/variant/variant.hpp>
#include <boost/variant/get.hpp>
//#include "tbb/tick_count.h"
//using namespace tbb;

typedef map<string, boost::variant<int, string, double, float, bool> > myparam;
typedef pair<string, boost::variant<int, string, double, float, bool> > parameterPair;


namespace DISCRETEOPTHOCR{

  class GroupDiscreteCostFunction: public SRegDiscreteCostFunction
  {
    
  public:
    //
    // Constructor.
    //
    GroupDiscreteCostFunction();
    // ~GroupDiscreteCostFunction();
    //// SET UP //////////////////

    void set_parameters(myparam & ); 
    /// neighbouhood info
    void set_relations(const boost::shared_ptr<RELATIONS> &CONTROL,const boost::shared_ptr<RELATIONS> &TARG);
    void get_spacings();
    void set_meshes(const NEWMESH::newmesh & target,const NEWMESH::newmesh & source, const NEWMESH::newmesh & GRID, int num=1){_TEMPLATE=target; 

      _ORIG=source;
      num_subjects=num;
      VERTICES_PER_SUBJ=GRID.nvertices();
      TRIPLETS_PER_SUBJ=GRID.ntriangles();
      MVD_LR=Calculate_MVD(GRID);
      cout << " costfunction set meshes " << num << endl;
      for (int i=0;i<num;i++){
	_DATAMESHES.push_back(source);
	_CONTROLMESHES.push_back(GRID);
      }
    } 

    //// INITIALISATION //////////////////

    virtual void initialize(int numNodes, int numLabels, int numPairs, int numTriplets = 0,int numQuartets = 0); // quartets not used yet so no code for them below
    void define_template_patches();
    void resample_to_template(); /// resamples all data ontotemplate

    /* inline void set_matlab_path(string s){_matlabpath=s;   
      cout << " in group set matlab path " << _matlabpath.c_str() << endl;
      char filename[1000];
      sprintf(filename,"addpath('%s')",_matlabpath.c_str());
      engEvalString(eng,filename);
     
   }
    */
    //////////////// Updates ///////////////////////
    void reset_source(const NEWMESH::newmesh & source, int num=0){_DATAMESHES[num]=source;}
    virtual void reset_CPgrid(const NEWMESH::newmesh & grid,int num=0){_CONTROLMESHES[num]=grid;}


    double computeQuartetCost(int quartet, int labelA, int labelB, int labelC,int labelD);

    double computeTripletCost(int triplet, int labelA, int labelB, int labelC);
    double computePairwiseCost(int pair, int labelA, int labelB);

    void resample_patches();
    void resampler_worker_function(const int &, const int &,const vector<bool> &);
    map<int,float>  resample_onto_template(const int &,const int &,const Pt &, const vector<int> &);

  protected:

    //////////////////////// MESHES //////////////////////
    vector<NEWMESH::newmesh> _DATAMESHES; // TARGET MESH
    vector<NEWMESH::newmesh> _CONTROLMESHES; // TARGET MESH
    newmesh _TEMPLATE;
    /////////////////////// NEIGHBOURHOOD INFO ////////////////////
    vector<boost::shared_ptr<RELATIONS> > _CONTROLRELATIONS; // hold control grid neighbours of each source vertex
    vector<boost::shared_ptr<RELATIONS> > _TEMPLATERELATIONS; // hold target grid neighbours of each source vertex
    vector<vector<int> > TEMPLATEPTS; 
    vector<ColumnVector> SPACINGS;

    //////////////////////DATA//////////////////
    vector<Matrix>  RESAMPLEDDATA;
    vector<map<int,float> > PATCHDATA;
    //double **L1data;

    //Engine *eng; /// matlab engine

    double MVD_LR;

    float _lambdapairs;
    int num_subjects;
    int TRIPLETS_PER_SUBJ;
    int VERTICES_PER_SUBJ;

    bool _quadcost;
    bool _setpairs;

  };
  

}
