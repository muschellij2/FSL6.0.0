/*  DiscreteModel.h

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
#include "DiscreteCostFunction.h"
#include <vector>


namespace DISCRETEOPT{
  
  void sort_nodeids(int *, int );

  class DiscreteModel
  {
  public:
    
    /**
     * Constructor.
     */
  DiscreteModel()
    : m_num_nodes(0), m_num_labels(0), m_num_pairs(0), m_num_triplets(0), m_num_quartets(0), m_lambda(0.0), labeling(0), pairs(0), triplets(0),  quartets(0)
      {
	_numthreads=1;
      }
  
    /**
     * Constructor.
     */
  DiscreteModel(myparam P)
    : m_num_nodes(0), m_num_labels(0), m_num_pairs(0), m_num_triplets(0), m_num_quartets(0), m_lambda(0.0), labeling(0), pairs(0), triplets(0), quartets(0)
      {
	_numthreads=1;
      }

    //====================================================================================================================//
    /**
     * Destructor.
     */
    virtual ~DiscreteModel()
  {
    if(labeling) { delete[] labeling; labeling = 0; }
    if(pairs) { delete[] pairs; pairs = 0; }
    if(triplets) { delete[] triplets; triplets = 0; }
    if(quartets) { delete[] quartets; quartets= 0; }

    // if(costfct) { delete costfct; costfct = 0; }
  }
    
    /**
     * Sets the weighting between unary and pairwise costs.
     */
    void setLambda(double lambda) { m_lambda = lambda; }

    /**
     * Returns the number of points.
     */
    int getNumNodes() const { return m_num_nodes; }
    
    
    /**
     * Returns the number of labels.
     */
    int getNumLabels() const { return m_num_labels; }
    
    /**
     * Returns the number of pairs.
     */
    int getNumPairs() const { return m_num_pairs; }
    
    /**
     * Returns the number of triplets.
     */
    int getNumTriplets() const { return m_num_triplets; }
    
    /**
     * Returns the number of quartets.
     */
    int getNumQuartets() const { return m_num_quartets; }
    
    /**
     * Returns the weighting between unary and pairwise costs.
     */
    double getLambda() const { return m_lambda; }
    
    /**
     * Returns the labeling.
     */
    int* getLabeling() { return labeling; }
    const int* getLabeling() const { return labeling; }

    /**
     * Returns the pairs.
     */
    const int* getPairs() const { return pairs; }
    
    /**
     * Returns the triplets.
     */
    const int* getTriplets() const { return triplets; }
    
    /**
     * Returns the quartets.
     */
    const int* getQuartets() const { return quartets; }
   
    /**
     * Resets the labeling.
     */
    void resetLabeling();
    
    /**
     * Returns the cost function.
     */
    virtual boost::shared_ptr<DiscreteCostFunction> getCostFunction() =0;// this class is pure virtual
    
    //virtual void InitializeUnaryCost(int)=0;

    /**
     * Invokes the unary costs look-up table computation.
     */
    virtual void computeUnaryCosts(){};

    /**
     *  Invokes the  unary potential for a the given node
     */
    virtual double computeUnaryCost(int node, int label){ return 0; };

    /**
     * Invokes the pairwise costs look-up table computation.
     */
    virtual void computePairwiseCosts(){};
    /**
     *  Invokes the  pairwise potential for a the given pair and labels.
     */
    virtual double computePairwiseCost(int pair, int labelA, int labelB){ return 0;};
 
    /**
     * Invokes the triplet costs look-up table computation.
     */
    //virtual void computeTripletCosts(){};
    
    /**
     * Invokes the  triplet potential for a the given triplet and labels.
     */
    virtual double computeTripletCost(int triplet, int labelA, int labelB, int labelC) { return 0; }
   
    /**
     * Invokes the pair weights initialization.
     */
    virtual void initPairWeights(const double *pairweights = 0){};
    
    /**
     * Invokes the triplet weights initialization.
     */
    virtual void initTripletWeights(const double *tripletweights = 0){};
    
    /**
     *  Invokes the  quartet costs look-up table.
     */
    //  virtual void computeQuartetCosts() {};
    
    /**
     * Invokes the  quartet potential for a the given triplet and labels.
     */
    virtual double computeQuartetCost(int quartet, int labelA, int labelB, int labelC, int labelD) { return 0; }
    
    /**
     * Enables the memory for the pairwise potentials computation.
     */
    //void enablePairwiseMemory() { if(costfct) costfct->enablePairwiseMemory(); }
    
    /**
     * Enables the memory for the triplet potentials computation.
     */
    //void enableTripletMemory() { if(costfct) costfct->enableTripletMemory(); }
    
    /**
     * Evaluates the total cost for the zero labeling.
     */
    virtual double evaluateTotalCostSumZeroLabeling(){return 0;};
    
    /**
     * Evaluates the total cost w.r.t to the current labeling.
     */
    virtual double evaluateTotalCostSum(){return 0;}; 

    /**
     * Evaluates the unary cost sum w.r.t to the current labeling.
     */
    virtual double evaluateUnaryCostSum(){return 0;};
    
    /**
     * Evaluates the pairwise cost sum w.r.t to the current labeling.
	*/
    virtual double evaluatePairwiseCostSum(){return 0;};
    
    /**
     * Evaluates the triplet cost sum w.r.t to the current labeling.
     */
    virtual double evaluateTripletCostSum(){return 0;};
    
    /**
     * Evaluates the quartet cost sum w.r.t to the current labeling.
     */
    virtual double evaluateQuartetCostSum(){return 0;};
   
    /**
     * Computes the uncertainties.
     */
    // virtual void computeUncertainties(int level = 0, int iteration = 0) =0;
    
    /**
     * Applies the current labeling to the model configuration.
     */
    virtual void applyLabeling() =0;
    
    /**
     * Applies the given labeling to the model configuration.
     */
    virtual void applyLabeling(int *discreteLabeling) =0;
    
    /**
     * Frees memory which is not needed any more after optimization.
     */
    // virtual void freeMemory() { deleteCostFunction(); }
    
    /**
     * Prints the labeling.
     */

    virtual  void report()=0;

        int getNumthreads(){return _numthreads;}


    void printLabeling();
    //  virtual void Initialize()=0;
    virtual void Initialize(const newmesh &)=0;
    virtual void Initialize()=0;
    virtual void setupCostFunction()=0;
    virtual void set_parameters(myparam PAR)=0;
  protected:
   
    
    /**
     * Deletes the cost function.
     */
    //void deleteCostFunction() { if(costfct) { delete costfct; costfct = 0; } }
    
    /**
     * Initializes the labeling buffer.
     */
    void initLabeling();
    
    int m_num_nodes;						///< Number of model nodes (e.g. grid nodes, variables, etc).
    int m_num_labels;						///< Number of labels.
    int m_num_pairs;						///< Number of node pairs.
    int m_num_triplets;						///< Number of node triplets.
    int m_num_quartets;						///< Number of node quartets (currently none).


    float m_lambda;							///< Weighting between unary and pairwise costs.
    
    int *labeling;							///< Labeling array.
    int* pairs;								///< Node pairs array.
    int *triplets;							///< Node triplets array.
    int *quartets;							///< Node triplets array.
    int _numthreads;

    string m_outdir;
    bool m_verbosity;


    					///< Cost function.
  };
  
  class SRegDiscreteModel: public DiscreteModel
  {
    
  public:
    /**
     * Constructor.
     */
    SRegDiscreteModel(){ 
      m_CPres=2; m_SGres=4; m_simmeasure=2; m_multivariate=false; m_verbosity=false; m_outdir=""; m_debug=false;  m_regoption=2;
      _pairwise=false; _estquartet=false; m_triclique=false;
      m_inputrel=boost::shared_ptr<RELATIONS >(new RELATIONS()); 
      m_cp_neighbourhood=boost::shared_ptr<RELATIONS >(new RELATIONS ());
    };

    /**
     * Constructor.
     */
    SRegDiscreteModel(myparam& PAR){
      m_CPres=2; m_SGres=4; m_simmeasure=2; m_multivariate=false; m_verbosity=false; m_outdir="";m_debug=false;  m_regoption=2;
      _pairwise=false; _estquartet=false; m_triclique=false;
      set_parameters(PAR);
      initialize_cost_function(m_multivariate,m_simmeasure,PAR);
      m_inputrel=boost::shared_ptr<RELATIONS >(new RELATIONS());
      m_cp_neighbourhood=boost::shared_ptr<RELATIONS >(new RELATIONS ());
   
    }

    /**
     * Destructor.
     */
    ~SRegDiscreteModel(){};
    
    /**
     * Returns the cost function.
     */
    boost::shared_ptr<DiscreteCostFunction> getCostFunction() { boost::shared_ptr<DiscreteCostFunction> dcostfct=costfct; return dcostfct; }// upcast
    
    //void InitializeUnaryCost(int label){ if(costfct) costfct->InitializeUnaryCost(label); };

    /**
     * Invokes the unary costs look-up table computation.
     */
    void computeUnaryCosts(){ if(costfct) costfct->computeUnaryCosts(); }

    /**
     * Invokes the unary potential for a the given node
     */
    double computeUnaryCost(int node, int label) {return (costfct) ? costfct->computeUnaryCost(node,label): 0.0f; }

    /**
     * Invokes the pairwise costs look-up table computation.
     */
    void computePairwiseCosts() { if(costfct) costfct->computePairwiseCosts(pairs); }

    double computePairwiseCost(int pair, int labelA, int labelB){return (costfct) ? costfct->computePairwiseCost(pair,labelA,labelB): 0.0f; }
    /**
     * Invokes the pair weights initialization.
     */
    void initPairWeights(const double *pairweights = 0) { if(costfct) costfct->initPairWeights(pairweights); }
    
      /**
     * Invokes the pairwise costs look-up table computation.
     */
    //  void computeTripletCosts() { if(costfct) costfct->computeTripletCosts(); }

    double computeTripletCost(int triplet, int labelA, int labelB, int labelC){return (costfct) ?costfct->computeTripletCost(triplet,labelA,labelB,labelC): 0.0f; }

    /**
     * Computes the quartet costs look-up table.
     */
    //    void computeQuartetCosts(){ if(costfct) costfct->computeQuartetCosts(); }
    
    /**
     * Computes the quartet potential for a the given triplet and labels.
     */
    double computeQuartetCost(int quartet, int labelA, int labelB, int labelC, int labelD) { return (costfct) ? costfct->computeQuartetCost(quartet,labelA,labelB,labelC,labelD): 0.0f;}

    /**
     * Invokes the pair weights initialization.
     */
    //   void initTripletWeights(const double *tripletweights = 0) { if(costfct) costfct->initTripletWeights(tripletweights); }
    
 
    /**
     * Evaluates the total cost for the zero labeling.
     */
    double evaluateTotalCostSumZeroLabeling(){ return (costfct) ? costfct->evaluateTotalCostSumZeroLabeling() : 0.0f; }
    
    /**
     * Evaluates the total cost w.r.t to the current labeling.
     */
    double evaluateTotalCostSum(){ return (costfct) ? costfct->evaluateTotalCostSum(labeling,pairs,triplets,quartets) : 0.0f; }

    /**
     * Evaluates the unary cost sum w.r.t to the current labeling.
     */
    double evaluateUnaryCostSum(){ if(costfct) return costfct->evaluateUnaryCostSum(labeling); else return 0.0f; }
    
    /**
     * Evaluates the pairwise cost sum w.r.t to the current labeling.
	*/
    double evaluatePairwiseCostSum(){ if(costfct) return costfct->evaluatePairwiseCostSum(labeling,pairs); else return 0.0f; }
       

    /**
     * Evaluates the triplet cost sum w.r.t to the current labeling.
     */
    double evaluateTripletCostSum() { if(costfct) return costfct->evaluateTripletCostSum(labeling,triplets); else return 0.0f; }

    /**
     * Evaluates the quartet cost sum w.r.t to the current labeling.
     */
    double evaluateQuartetCostSum() { if(costfct) return costfct->evaluateQuartetCostSum(labeling,triplets); else return 0.0f; }

    virtual void set_parameters(myparam PAR);

    /*                                     INITIALIZE MODEL                         */

    virtual void set_meshspace(const NEWMESH::newmesh & target,const NEWMESH::newmesh & source, const int num=1){ 
m_TARGET=target; m_SOURCE=source; }
    // void set_anatomical_meshspace(const boost::shared_ptr<NEWMESH::newmesh> &ref_sphere, const boost::shared_ptr<NEWMESH::newmesh> &ref_anat,const NEWMESH::newmesh & source_anat){  if(costfct.get()){ costfct->set_anatomical(ref_sphere, ref_anat,source_anat);}}
    void set_anatomical_meshspace(const NEWMESH::newmesh &ref_sphere, const NEWMESH::newmesh &ref_anat,const NEWMESH::newmesh & source_sphere, const NEWMESH::newmesh & source_anat){  if(costfct.get()){ costfct->set_anatomical(ref_sphere, ref_anat,source_sphere,source_anat);}}
    void set_anatomical_neighbourhood(const vector<map<int,double> > & weights, const vector<vector<int> > neighbourhood){  if(costfct.get()){ costfct->set_anatomical_neighbourhood(weights,neighbourhood);}}

    void set_featurespace(const boost::shared_ptr<featurespace> &FEATURES, bool  _concatenate=false){  if(FEATURES->get_dim()>1) m_multivariate=true;
      if(costfct.get()){ costfct->set_featurespace(FEATURES,_concatenate); }     
      else {throw  DISCRETEOPTHOCRException("Discrete Model:: You have not initialised the discrete costfunction before setting the featurespace");}
    }

    void set_L1path(string s){if(costfct.get()){ costfct->set_matlab_path(s); }}
    //// costfunction weighting combines source and reference weightings at beginning of optimisation iteration - will not be 100% accurate but will remove any sensitivity of the label choices to weighting   
    void setupCostFunctionWeighting(const Matrix & Weight){costfct->set_dataaffintyweighting(Weight);};

    /// source needs to be reset after every iteration of discrete optimisation
    virtual void reset_meshspace(const NEWMESH::newmesh & source, int num=0){  
      
      m_SOURCE=source;
      if(costfct.get()){costfct->reset_source(source);}
      else{throw  DISCRETEOPTHOCRException("Discrete Model:: You cannot reset the source mesh withou initialising the discrete costfunction");} 
    }

    virtual void reset_CPgrid(const NEWMESH::newmesh & grid, int num=0){  
      m_CPgrid=grid;
    }

    virtual void warp_CPgrid(NEWMESH::newmesh & START, NEWMESH::newmesh & END, int num=0){ barycentric_mesh_interpolation(m_CPgrid,START,END); 	unfold(m_CPgrid);}
    
    void initialize_cost_function(const bool &MV, const int & sim, myparam &P);
    void Initialize_sampling_grid();
    void label_sampling_grid(const int &, const double &, NEWMESH::newmesh &);
 
    virtual void Initialize(const newmesh &);    

    virtual void set_debug(){m_debug=true; costfct->debug();} // for debuging 

    boost::shared_ptr<RELATIONS> get_cp_neighbourhood(){return m_cp_neighbourhood;};

    NEWMESH::newmesh get_SOURCE(){ if(costfct.get()) return costfct->get_SOURCE(); else{ cout << " no costfunction return inital source mesh " << endl;return m_SOURCE; } }

    NEWMESH::newmesh get_TARGET(){return m_TARGET;}
    virtual NEWMESH::newmesh get_CPgrid(int num=0){return m_CPgrid;}
    void report(){if(costfct.get()){ costfct->report();}}

  protected:

    NEWMESH::newmesh m_TARGET; // TARGET MESH
    NEWMESH::newmesh m_SOURCE; // SOURCE MESH   
    NEWMESH::newmesh m_CPgrid; ///// CONTROL POINT GRID
    NEWMESH::newmesh m_samplinggrid;

    boost::shared_ptr<RELATIONS> m_cp_neighbourhood; // hold control grid neighbours of each source vertex
    boost::shared_ptr<RELATIONS> m_inputrel; // hold target grid neighbours of each source vertex
   
    int m_CPres;//control grid resolution
    int m_SGres;// sampling grid resolution
    int m_iter;// iteration of the discrete optimisation
    int m_centroid; /// used for selecting which sampling grid vertex will form the center of the sampling grid
    int m_simmeasure; // sim measure i.e. correlation
    int m_regoption; // sim measure i.e. correlation

    double m_maxs_dist; // define maximum distance between the centre of the sampling grid and the furthest label
    double MVD;
    bool m_multivariate;
    bool m_debug;
    bool m_triclique;
    bool _pairwise;
    bool _estquartet;

    vector<Pt> m_samples;  // samples based on  vertices of sampling grid
    vector<Pt> m_barycentres; // samples based on barycentres of sampling grid
    vector<Pt> m_labels; // labels iterates between samples and barycnetres and is the label set used within cosfct

    vector<Matrix> m_ROT; // rotates sampling grid to each control point

    boost::shared_ptr<SRegDiscreteCostFunction> costfct;  // costfunction object
  };

  class AffineSRegDiscreteModel: public SRegDiscreteModel
  {
    
  public:
    /**
     * Constructor.
     */
    AffineSRegDiscreteModel(){};
  AffineSRegDiscreteModel(myparam & P):SRegDiscreteModel(P){};
    void setupCostFunction() ;
  };

  class NonLinearSRegDiscreteModel: public SRegDiscreteModel
  {     
  protected:

  public:
    /**
     * Constructor.
     */
    NonLinearSRegDiscreteModel(){};
    NonLinearSRegDiscreteModel(myparam & P):SRegDiscreteModel(P){  set_parameters(P); };
    
  
    void applyLabeling(){applyLabeling(labeling);} ;
    void applyLabeling(int *discreteLabeling);

    void estimate_pairs();
    void estimate_triplets();

    void Initialize(const newmesh &);
    void Initialize(){};

    void get_rotations(vector<Matrix> &);
    void setupCostFunction();
   
  };


  class MetricDistortionDiscreteModel: public NonLinearSRegDiscreteModel
  {     
  protected:


  public:
    /**
     * Constructor.
     */
    MetricDistortionDiscreteModel(newmesh ANAT, int SAMP, double shear, double bulk, double exp){ m_SGres=SAMP;  costfct=boost::shared_ptr<SRegDiscreteCostFunction>(new MetricDistortionDiscreteCostFunction(ANAT,shear,bulk,exp)); m_num_pairs=0;};
    
    void Initialize(const newmesh &);
    void Initialize(){};
    void setupCostFunction();
   

   
  };
  class GroupSegmentationModel: public DiscreteModel
  {     
    newmesh m_SOURCE; // constant surface for all label maps  - this enables us to have one RELATIONS matrix

    boost::shared_ptr<RELATIONS> m_neighbourhood; // hold control grid neighbours of each source vertex
   					///< Number of Between subject (BS) node pairs.   // for this models pairs and BS pairs are eqivalent (for c less than number of neighbours) but this may change for joint approach
    boost::shared_ptr<SegmentationCostFunction> costfct;  // costfunction object
    Matrix SEGMENTATION;

    myparam params;

    int m_count;
    int m_num_subjects;
  public:
    /**
     * Constructor.
     */
  GroupSegmentationModel(): m_count(1),m_num_subjects(0){  costfct=boost::shared_ptr<SegmentationCostFunction>(new SegmentationCostFunction());};
    GroupSegmentationModel(myparam & P){ m_count=1; set_parameters(P); 
      costfct=boost::shared_ptr<SegmentationCostFunction>(new SegmentationCostFunction());
      costfct->set_parameters(P);
    }
    
    /**
     * Destructor.
     */
    ~GroupSegmentationModel(){}

    void Initialize();
   
    void setupCostFunction();

    void set_parameters(myparam PAR){
      params=PAR; myparam::iterator it;
      it=PAR.find("bscount");m_count=boost::get<int>(it->second); 
      it=PAR.find("verbosity");m_verbosity=boost::get<bool>(it->second);
      it=PAR.find("outdir");m_outdir=boost::get<string>(it->second);
    }

    void set_SURFACE(const newmesh & MESH){m_SOURCE=MESH;     
      if(costfct.get()){costfct->set_SURFACE(MESH);}
      else {throw  DISCRETEOPTHOCRException("Group segmentation Model:: You have not initialised the discrete costfunction before passing the SURFACE");}
    };

    void set_DATA(const vector<Matrix> &DATA){ 
      m_num_subjects=DATA.size(); 
      if(costfct.get()){ costfct->set_DATA(DATA); }
      else {throw  DISCRETEOPTHOCRException("Group segmentation Model:: You have not initialised the discrete costfunction before passing the data");}}

    boost::shared_ptr<DiscreteCostFunction> getCostFunction() { boost::shared_ptr<DiscreteCostFunction> dcostfct=costfct; return dcostfct; }// upcast
    void applyLabeling(){applyLabeling(labeling);} ;
    void applyLabeling(int *discreteLabeling);

    Matrix get_segmentation(){return SEGMENTATION;}
    /**
     * Invokes the unary costs look-up table computation.
     */
    void computeUnaryCosts(){ if(costfct) costfct->computeUnaryCosts(); }

    /**
     * Invokes the unary potential for a the given node
     */
    double computeUnaryCost(int node, int label) {return (costfct) ? costfct->computeUnaryCost(node,label): 0.0f; }

    /**
     * Invokes the pairwise costs look-up table computation.
     */
    void computePairwiseCosts() { if(costfct) costfct->computePairwiseCosts(pairs); }

    double computePairwiseCost(int pair, int labelA, int labelB){return (costfct) ? costfct->computePairwiseCost(pair,labelA,labelB): 0.0f; }
  
    /**
     * Invokes the pairwise costs look-up table computation.
     */
    //  void computeTripletCosts() { if(costfct) costfct->computeTripletCosts(); }

    double computeTripletCost(int triplet, int labelA, int labelB, int labelC){return (costfct) ?costfct->computeTripletCost(triplet,labelA,labelB,labelC): 0.0f; }

   
    /**
     * Evaluates the total cost w.r.t to the current labeling.
     */
    double evaluateTotalCostSum(){ return (costfct) ? costfct->evaluateTotalCostSum(labeling,pairs,triplets,quartets) : 0.0f; }
    

  };

}
