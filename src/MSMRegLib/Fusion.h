/*  Fusion.h

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

#include "QPBO/QPBO.hpp"
#include "ELC/ELC.h"
#if !defined(__DiscreteModel_h)
#define __DiscreteModel_h
#include "DiscreteModel.h"
#endif

#include <iostream>

#ifdef HAS_TBB

#include "tbb/parallel_for.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/tick_count.h"
using namespace tbb;
#endif

using namespace detail;

typedef	double REAL;
struct UnaryData	{ REAL buffer[2]; };
struct PairData		{ REAL buffer[4]; };
struct TripletData	{ REAL buffer[8]; };
struct QuartetData	{ REAL buffer[16]; };


namespace ELCReduce {
  template<typename T> class PBF;
};

template<typename T> class QPBO;


namespace DISCRETEOPT{

#ifdef HAS_TBB
  class computeTripletCosts{
    const boost::shared_ptr<DiscreteModel> my_energy;
    const int* my_trip;
    const int* my_labels;
    const int my_current_label;

    std::vector<TripletData> & my_DATA;

  public:
    // constructor copies the arguments into local storage
  computeTripletCosts(const boost::shared_ptr<DiscreteModel> & energy, const int *trip, const  int *labels, const int & label, std::vector<TripletData> & DATA) :
    my_energy(energy), my_trip(trip), my_labels(labels), my_current_label(label),my_DATA(DATA)
    {}
    // overload () so it does a vector multiply
    void operator() (const blocked_range<int> &r) const {
      for(size_t triplet=r.begin(); triplet!=r.end(); triplet++){
	const int nodeA = my_trip[triplet*3];
	const int nodeB = my_trip[triplet*3+1];
	const int nodeC = my_trip[triplet*3+2];

	my_DATA[triplet].buffer[0] = my_energy->computeTripletCost(triplet,my_labels[nodeA],my_labels[nodeB],my_labels[nodeC]);	//000
	my_DATA[triplet].buffer[1] = my_energy->computeTripletCost(triplet,my_labels[nodeA],my_labels[nodeB],my_current_label);			//001
	my_DATA[triplet].buffer[2] = my_energy->computeTripletCost(triplet,my_labels[nodeA],my_current_label,my_labels[nodeC]);			//010
	my_DATA[triplet].buffer[3] = my_energy->computeTripletCost(triplet,my_labels[nodeA],my_current_label,my_current_label);						//011
	my_DATA[triplet].buffer[4] = my_energy->computeTripletCost(triplet,my_current_label,my_labels[nodeB],my_labels[nodeC]);			//100
	my_DATA[triplet].buffer[5] = my_energy->computeTripletCost(triplet,my_current_label,my_labels[nodeB],my_current_label);						//101 
	my_DATA[triplet].buffer[6] = my_energy->computeTripletCost(triplet,my_current_label,my_current_label,my_labels[nodeC]);						//110 
	my_DATA[triplet].buffer[7] = my_energy->computeTripletCost(triplet,my_current_label,my_current_label,my_current_label);
      }
    }

  };

 class computeQuartetCosts{
    const boost::shared_ptr<DiscreteModel> my_energy;
    const int* my_quartet;
    const int* my_labels;
    const int my_current_label;

    std::vector<QuartetData> & my_DATA;

  public:
    // constructor copies the arguments into local storage
 computeQuartetCosts(const boost::shared_ptr<DiscreteModel> & energy, const int *quart, const  int *labels, const int & label, std::vector<QuartetData> & DATA) :
    my_energy(energy), my_quartet(quart), my_labels(labels), my_current_label(label),my_DATA(DATA)
    {}
    // overload () so it does a vector multiply
    void operator() (const blocked_range<int> &r) const {
      for(size_t quartet=r.begin(); quartet!=r.end(); quartet++){
	const int nodeA = my_quartet[quartet*4];
	const int nodeB = my_quartet[quartet*4+1];
	const int nodeC = my_quartet[quartet*4+2];
	const int nodeD = my_quartet[quartet*4+3];
	
	my_DATA[quartet].buffer[0]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_labels[nodeB],my_labels[nodeC],my_labels[nodeD]);
	my_DATA[quartet].buffer[1]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_labels[nodeB],my_labels[nodeC],my_current_label);
	my_DATA[quartet].buffer[2]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_labels[nodeB],my_current_label,my_labels[nodeD]);
	my_DATA[quartet].buffer[3]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_labels[nodeB],my_current_label,my_current_label);
	
	my_DATA[quartet].buffer[4]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_current_label,my_labels[nodeC],my_labels[nodeD]);
	my_DATA[quartet].buffer[5]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_current_label,my_labels[nodeC],my_current_label);
	my_DATA[quartet].buffer[6]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_current_label,my_current_label,my_labels[nodeD]);
	my_DATA[quartet].buffer[7]  = my_energy->computeQuartetCost(quartet,my_labels[nodeA],my_current_label,my_current_label,my_current_label);

	my_DATA[quartet].buffer[8]  = my_energy->computeQuartetCost(quartet,my_current_label,my_labels[nodeB],my_labels[nodeC],my_labels[nodeD]);
	my_DATA[quartet].buffer[9]  = my_energy->computeQuartetCost(quartet,my_current_label,my_labels[nodeB],my_labels[nodeC],my_current_label);
	my_DATA[quartet].buffer[10] = my_energy->computeQuartetCost(quartet,my_current_label,my_labels[nodeB],my_current_label,my_labels[nodeD]);
	my_DATA[quartet].buffer[11] = my_energy->computeQuartetCost(quartet,my_current_label,my_labels[nodeB],my_current_label,my_current_label);

	my_DATA[quartet].buffer[12] = my_energy->computeQuartetCost(quartet,my_current_label,my_current_label,my_labels[nodeC],my_labels[nodeD]);
	my_DATA[quartet].buffer[13] = my_energy->computeQuartetCost(quartet,my_current_label,my_current_label,my_labels[nodeC],my_current_label);
	my_DATA[quartet].buffer[14] = my_energy->computeQuartetCost(quartet,my_current_label,my_current_label,my_current_label,my_labels[nodeD]);
	my_DATA[quartet].buffer[15] = my_energy->computeQuartetCost(quartet,my_current_label,my_current_label,my_current_label,my_current_label);


      }
    }

  };
class computeUnaryCosts{
    const boost::shared_ptr<DiscreteModel> my_energy;
    const int* my_labels;
    const int my_current_label;

    std::vector<UnaryData> & my_DATA;
    double &my_label_diff;

  public:
    // constructor copies the arguments into local storage
 computeUnaryCosts(const boost::shared_ptr<DiscreteModel> & energy, const  int *labels, const int & label, std::vector<UnaryData> & DATA, double & diff) :
    my_energy(energy), my_labels(labels), my_current_label(label),my_DATA(DATA), my_label_diff(diff)
    {}
    // overload () so it does a vector multiply
    void operator() (const blocked_range<int> &r) const {
      for(size_t node=r.begin(); node!=r.end(); node++){

	my_DATA[node].buffer[0] = my_energy->computeUnaryCost(node,my_labels[node]);	//0
        my_DATA[node].buffer[1] = my_energy->computeUnaryCost(node,my_current_label);
	my_label_diff+=abs(my_current_label-my_labels[node]);
      }
    }

  };
#endif

 enum Reduction
  {
    ELC_HOCR,
    ELC_APPROX,
    HOCR
  };
  class Fusion
  {
  public:
    /**
     * Constructor.
     */
    Fusion() {}
    /**
     * Destructor.
     */
    ~Fusion() {}
  
    static void reduce_and_convert(ELCReduce::PBF<REAL>&, QPBO<REAL>&, Reduction);
    /**
     * Runs the optimization.
     */
    static double optimize(boost::shared_ptr<DiscreteModel> energy, Reduction reductionMode, bool debug=false);
  };

}
