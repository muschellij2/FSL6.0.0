/*  Fusion.cpp

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
#include "Fusion.h"
#include "QPBO/QPBO.hpp"
#include "ELC/ELC.h"
#include <iostream>
#include <time.h>

#define NUM_SWEEPS 2
#define MAX_IMPROVEMENTS 0

using namespace ELCReduce;

namespace DISCRETEOPT{


  void Fusion::reduce_and_convert(PBF<REAL>& pbf, QPBO<REAL>& qpbo, Reduction mode)
  {
      switch(mode)
      {
      case ELC_HOCR:
      {
        PBF<REAL> qpbf;
       
        pbf.reduceHigher();
	pbf.toQuadratic(qpbf, pbf.maxID()+1); // Reduce the remaining higher-order terms using HOCR adding auxiliary variables
	qpbf.convert(qpbo, qpbf.maxID()+1);
        pbf.clear();
        qpbf.clear();
        break;
      }
      case ELC_APPROX:
      { 
        //PBF<REAL> qpbf;
	//qpbf = pbf;
	pbf.reduceHigherApprox();

        pbf.convert(qpbo, pbf.maxID()+1);
        pbf.clear();

        break;
      }
      case HOCR:
      {
        PBF<REAL> qpbf;
        pbf.toQuadratic(qpbf, pbf.maxID()+1); // Reduce to Quadratic pseudo-Boolean function using HOCR.
        qpbf.convert(qpbo, qpbf.maxID()+1);
        pbf.clear();
        qpbf.clear();
        break;

      }
     
      }
    }

//====================================================================================================================//
  double Fusion::optimize( boost::shared_ptr<DiscreteModel> energy, Reduction reductionMode, bool verbose )
{
  const int numNodes    = energy->getNumNodes();
  const int numPairs    = energy->getNumPairs();
  const int numTriplets = energy->getNumTriplets();
  const int numQuartets = energy->getNumQuartets();
  const int numthreads  = energy->getNumthreads();
  const int *pairs    = energy->getPairs();
  const int *triplets = energy->getTriplets();
  const int *quartets = energy->getQuartets();

  const int numGraphNodes = numNodes + numTriplets + 5 * numQuartets;
  const int numGraphEdges = numPairs + 3 * numTriplets + 17 * numQuartets;
  
  
  QPBO<REAL> qpbo(numGraphNodes,numGraphEdges) ; //numGraphEdges);

  double max_unlabeled_ratio = 0;
  double avg_unlabeled_ratio = 0;

  const int numSweeps = NUM_SWEEPS;
  const int numLabels = energy->getNumLabels();
  int *labeling = energy->getLabeling();

 
  double initEnergy = energy->evaluateTotalCostSum();
  
  
  double lastEnergy = initEnergy;
  double unlabeledMax = 0, unlabeledAvg = 0;
  double sumlabeldiff=0;

 
  for(int sweep = 0; sweep < numSweeps; ++sweep)
  {
    for(int label = 0; label < numLabels; ++label)
    {
      
      sumlabeldiff=0;

      PBF<REAL> pbf;
      int unlabeledNodes = 0;
      int improveCounter = 0;

      double ratioUnlabeled = 0;
      int nodesChanged = 0;
      
     
      
      std::vector<UnaryData> unary_data(numNodes);
    
      for(int node = 0; node < numNodes; ++node)
	{
	  unary_data[node].buffer[0] = energy->computeUnaryCost(node,labeling[node]);	//0
	  unary_data[node].buffer[1] = energy->computeUnaryCost(node,label);				//1
	
	  sumlabeldiff+=abs(label-labeling[node]);
	}
     

      if(sumlabeldiff>0){
	for(int node = 0; node < numNodes; ++node)
	  {
	   
	    pbf.AddUnaryTerm(node, unary_data[node].buffer[0], unary_data[node].buffer[1]);
	  }

    
	std::vector<PairData> pair_data(numPairs);

	for(int pair = 0; pair < numPairs; ++pair)
	  {
	    const int nodeA = pairs[pair*2];
	    const int nodeB = pairs[pair*2+1];
	   
	      pair_data[pair].buffer[0] = energy->computePairwiseCost(pair,labeling[nodeA],labeling[nodeB]);	//00
	      pair_data[pair].buffer[1] = energy->computePairwiseCost(pair,labeling[nodeA],label);			//01
	      pair_data[pair].buffer[2] = energy->computePairwiseCost(pair,label,labeling[nodeB]);			//10
	      pair_data[pair].buffer[3] = energy->computePairwiseCost(pair,label,label);	
					//11
	  }

	for(int pair = 0; pair < numPairs; ++pair)
	  {
	    const int nodeA = pairs[pair*2];
	    const int nodeB = pairs[pair*2+1];
	    int node_ids[2] = { nodeA, nodeB};

	    pbf.AddPairwiseTerm(node_ids[0], node_ids[1], pair_data[pair].buffer[0], pair_data[pair].buffer[1], pair_data[pair].buffer[2], pair_data[pair].buffer[3]);
	  }
      
	std::vector<TripletData> triplet_data(numTriplets);
	//std::vector<TripletData> triplet_data_par(numTriplets);
#ifdef HAS_TBB
	if(label==1)	cout << " has tbb " << endl;

	tbb::task_scheduler_init init(numthreads);
	parallel_for( blocked_range<int>(0,numTriplets), computeTripletCosts(energy,triplets,labeling,label,triplet_data) );

#else
	if(label==1)cout << "does not has tbb " << endl;

	for (int triplet = 0; triplet < numTriplets; ++triplet)
	  {
	    const int nodeA = triplets[triplet*3];
	    const int nodeB = triplets[triplet*3+1];
	    const int nodeC = triplets[triplet*3+2];
	 	    
	    triplet_data[triplet].buffer[0] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],labeling[nodeC]);	//000
	    triplet_data[triplet].buffer[1] = energy->computeTripletCost(triplet,labeling[nodeA],labeling[nodeB],label);			//001
	    triplet_data[triplet].buffer[2] = energy->computeTripletCost(triplet,labeling[nodeA],label,labeling[nodeC]);			//010
	    triplet_data[triplet].buffer[3] = energy->computeTripletCost(triplet,labeling[nodeA],label,label);						//011
	    triplet_data[triplet].buffer[4] = energy->computeTripletCost(triplet,label,labeling[nodeB],labeling[nodeC]);			//100
	    triplet_data[triplet].buffer[5] = energy->computeTripletCost(triplet,label,labeling[nodeB],label);						//101 
	    triplet_data[triplet].buffer[6] = energy->computeTripletCost(triplet,label,label,labeling[nodeC]);						//110 
	    triplet_data[triplet].buffer[7] = energy->computeTripletCost(triplet,label,label,label);								//111 

	  }
	
#endif

	for (int triplet = 0; triplet < numTriplets; ++triplet)
	  {
	    const int nodeA = triplets[triplet*3];
	    const int nodeB = triplets[triplet*3+1];
	    const int nodeC = triplets[triplet*3+2];
	      int node_ids[3] = { nodeA, nodeB, nodeC };
	      pbf.AddHigherTerm(3, node_ids, triplet_data[triplet].buffer);
	     
	  }
	

	std::vector<QuartetData> quartet_data(numQuartets);
#ifdef HAS_TBB
	
	parallel_for( blocked_range<int>(0,numQuartets), computeQuartetCosts(energy,quartets,labeling,label,quartet_data) );

#else

	for (int quartet = 0; quartet < numQuartets; ++quartet)
	  {
	    const int nodeA = quartets[quartet*4];
	    const int nodeB = quartets[quartet*4+1];
	    const int nodeC = quartets[quartet*4+2];
	    const int nodeD = quartets[quartet*4+3];
	    if(! (label==labeling[nodeA] && label==labeling[nodeB] && label==labeling[nodeC] && label==labeling[nodeD])){

	    quartet_data[quartet].buffer[0]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],labeling[nodeC],labeling[nodeD]);
	    quartet_data[quartet].buffer[1]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],labeling[nodeC],label);
	    quartet_data[quartet].buffer[2]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],label,labeling[nodeD]);
	    quartet_data[quartet].buffer[3]  = energy->computeQuartetCost(quartet,labeling[nodeA],labeling[nodeB],label,label);

	    quartet_data[quartet].buffer[4]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,labeling[nodeC],labeling[nodeD]);
	    quartet_data[quartet].buffer[5]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,labeling[nodeC],label);
	    quartet_data[quartet].buffer[6]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,label,labeling[nodeD]);
	    quartet_data[quartet].buffer[7]  = energy->computeQuartetCost(quartet,labeling[nodeA],label,label,label);
	    quartet_data[quartet].buffer[8]  = energy->computeQuartetCost(quartet,label,labeling[nodeB],labeling[nodeC],labeling[nodeD]);
	    quartet_data[quartet].buffer[9]  = energy->computeQuartetCost(quartet,label,labeling[nodeB],labeling[nodeC],label);
	    quartet_data[quartet].buffer[10] = energy->computeQuartetCost(quartet,label,labeling[nodeB],label,labeling[nodeD]);
	    quartet_data[quartet].buffer[11] = energy->computeQuartetCost(quartet,label,labeling[nodeB],label,label);
	    quartet_data[quartet].buffer[12] = energy->computeQuartetCost(quartet,label,label,labeling[nodeC],labeling[nodeD]);
	    quartet_data[quartet].buffer[13] = energy->computeQuartetCost(quartet,label,label,labeling[nodeC],label);
	    quartet_data[quartet].buffer[14] = energy->computeQuartetCost(quartet,label,label,label,labeling[nodeD]);
	    quartet_data[quartet].buffer[15] = energy->computeQuartetCost(quartet,label,label,label,label);
	    }
	  }
#endif

	for (int quartet = 0; quartet < numQuartets; ++quartet)
	  {
	    const int nodeA = quartets[quartet*4];
	    const int nodeB = quartets[quartet*4+1];
	    const int nodeC = quartets[quartet*4+2];
	    const int nodeD = quartets[quartet*4+3];
	    if(! (label==labeling[nodeA] && label==labeling[nodeB] && label==labeling[nodeC] && label==labeling[nodeD])){

	      int node_ids[4] = { nodeA, nodeB, nodeC, nodeD };
	      pbf.AddHigherTerm(4, node_ids, quartet_data[quartet].buffer);
	    }
	  }

       	qpbo.Reset();
	reduce_and_convert(pbf, qpbo, reductionMode);
	qpbo.MergeParallelEdges();
	qpbo.Solve();

	for(int node = 0; node < numNodes; ++node) {
	  if(labeling[node] != label)
	    {
	      if(qpbo.GetLabel(node) < 0) unlabeledNodes++;
	    }
	}

	//TRY QPBO-I to improve the solution
	if(unlabeledNodes > 0)
	  {
	    srand ( static_cast<unsigned int>(time(NULL)) );

	    int numTrials = MAX_IMPROVEMENTS;
	    const double ratioThresh = 0.3;
	    ratioUnlabeled = static_cast<double>(unlabeledNodes) / static_cast<double>(numNodes);
	    if (ratioUnlabeled < ratioThresh) numTrials = static_cast<int>(0.5+ratioUnlabeled*numTrials/ratioThresh);

	    if(MAX_IMPROVEMENTS > 0)
	      {
		for(int i = 0; i < numTrials; ++i)
		  {
		    if(qpbo.Improve()) improveCounter++;
		  }
	      }

	    if(ratioUnlabeled > unlabeledMax) unlabeledMax = ratioUnlabeled;
	    unlabeledAvg += ratioUnlabeled;
	  }
    
      
	  
      

	for(int node = 0; node < numNodes; ++node)
	  {
	    if(labeling[node] != label)
	      {
		if(qpbo.GetLabel(node) == 1) { labeling[node] = label; nodesChanged++; }
	      }
	  }

	double newEnergy;
	
	if(verbose){

	  newEnergy = energy->evaluateTotalCostSum();
	  energy->report();
	  std::cout << "LAB " << label << ":\t" << lastEnergy << " -> " << newEnergy << " / " << ratioUnlabeled * 100 << "% UNL / " << nodesChanged / static_cast<double>(numNodes) * 100 << "% CHN / IMP: " << improveCounter << std::endl;
	  lastEnergy = newEnergy;
      }

      
  
      }
    
      }
    }
    
    unlabeledAvg /= static_cast<double>(numLabels*numSweeps);
#ifndef PRINT_ENERGY
    lastEnergy = energy->evaluateTotalCostSum();
#endif
    
    return lastEnergy;
  
  }
}
