/*  DiscreteModel.cpp

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
#include "DiscreteModel.h"
#include <algorithm>
#include <iostream>


namespace DISCRETEOPT{

  void sort_nodeids(int *nodes, int length){
    std::vector<int> myvector (nodes, nodes+length);
    std::sort(myvector.begin(), myvector.begin()+length);

    for(unsigned int i=0;i<myvector.size();i++)
      nodes[i] =myvector[i];

  }

 
  //===
  //====================================================================================================================//
  
  //====================================================================================================================//
  void DiscreteModel::resetLabeling()
  {
    if(labeling) std::fill(labeling,labeling+m_num_nodes,0);
  }

  //====================================================================================================================//
  void DiscreteModel::initLabeling()
  {
    if(m_num_nodes)
      {
	if(labeling) delete[] labeling;
	labeling = new int[m_num_nodes];
	resetLabeling();
      }
  }

  //====================================================================================================================//
  void DiscreteModel::printLabeling()
  {
    std::cout << "Labeling: ";
    for(int i = 0; i < m_num_nodes; ++i)
      std::cout << labeling[i] << " ";
    std::cout << std::endl;
  }
  
  
 //====================================================================================================================//
  void SRegDiscreteModel::set_parameters(myparam PAR){
      myparam::iterator it;
      
      it=PAR.find("lambda");m_lambda=boost::get<float>(it->second);
      it=PAR.find("CPres");m_CPres=boost::get<int>(it->second);
      it=PAR.find("SGres");m_SGres=boost::get<int>(it->second); 
      it=PAR.find("simmeasure");m_simmeasure=boost::get<int>(it->second); 
      it=PAR.find("regularisermode");m_regoption=boost::get<int>(it->second);
      it=PAR.find("multivariate");m_multivariate=boost::get<bool>(it->second);
      it=PAR.find("verbosity");m_verbosity=boost::get<bool>(it->second);
      it=PAR.find("outdir");m_outdir=boost::get<string>(it->second);
      it=PAR.find("TriLikelihood");m_triclique=boost::get<bool>(it->second);
      it=PAR.find("numthreads");_numthreads=boost::get<int>(it->second);
      it=PAR.find("quartet");_estquartet=boost::get<bool>(it->second);

      if(m_regoption==1) _pairwise=true;


      
     
  } 
            
  void SRegDiscreteModel::Initialize(const newmesh &CONTROLGRID){
    ColumnVector vMAXmvd;
    double MVDmax=0.0;
    MVD=0;
    double tot=0;
    double dist=0;
    Pt CP;
    m_CPgrid=CONTROLGRID;
   
    /////////////// SET LOW RES DEFORMATION GRID & INITIALISE ASSOCIATED MRF PARAMS /////////////////////////
    m_num_nodes=m_CPgrid.nvertices(); m_num_pairs=0; 
    initLabeling();
 
    //////////////////////// CALCULATE (MEAN) VERTEX SEPARATIONS FOR EACH VERTEX
    vMAXmvd.ReSize(m_CPgrid.nvertices()); vMAXmvd=0;
    for (int k=0;k<m_CPgrid.nvertices();k++){
   
       CP=m_CPgrid.get_coord(k); 
      for (vector<int>::const_iterator it=m_CPgrid.nbegin(k);it !=m_CPgrid.nend(k);it++){

	dist= 2*RAD*asin((CP-m_CPgrid.get_coord(*it)).norm()/(2*RAD));
	MVD+=(CP-m_CPgrid.get_coord(*it)).norm();
	tot++; 

	if(dist > vMAXmvd(k+1))
	  vMAXmvd(k+1)=dist;
	if(vMAXmvd(k+1)>MVDmax) MVDmax=vMAXmvd(k+1);
      }
    }
    MVD/=tot;
    
    m_maxs_dist=0.4*Calculate_MaxVD(m_CPgrid);
    ////////////////// INITIALIZE COSTFCT ///////////////////
    
    if(costfct.get()){
      costfct->set_meshes(m_TARGET,m_SOURCE,m_CPgrid);
      costfct->set_spacings(vMAXmvd,MVDmax);
    }
    else {throw  DISCRETEOPTHOCRException("Discrete Model:: You have not initialised the discrete costfunction before initialising the model");}


    m_iter=1;
   

   
  }
  
  void SRegDiscreteModel::initialize_cost_function(const bool &MV, const int & sim, myparam &P){
      //////////////////////// CREATE COST FUNCTION //////////////////////
      if(MV){

	if(m_triclique)
	  costfct=boost::shared_ptr<SRegDiscreteCostFunction>(new HOMultivariateNonLinearSRegDiscreteCostFunction());
	else
	  costfct=boost::shared_ptr<SRegDiscreteCostFunction>(new MultivariateNonLinearSRegDiscreteCostFunction()); // choose CF based on choice of method of calculation of the data term
	
      }
      else if (sim==4) {throw  DISCRETEOPTHOCRException("DiscreteModel ERROR:: alpha MI not in this version of MSM yet");}  //costfct=boost::shared_ptr<SRegDiscreteCostFunction>(new alphaMINonLinearSRegDiscreteCostFunction());
      else{
	if(m_triclique)
	  costfct=boost::shared_ptr<SRegDiscreteCostFunction>(new HOUnivariateNonLinearSRegDiscreteCostFunction());
	else
	  costfct=boost::shared_ptr<SRegDiscreteCostFunction>(new UnivariateNonLinearSRegDiscreteCostFunction());}

  
      costfct->set_parameters(P);
      
    }

  void SRegDiscreteModel::Initialize_sampling_grid(){
    //////////////////////// LABELS USING HIGHER RES GRID //////////////////////
    m_samplinggrid.make_mesh_from_icosa(m_SGres); true_rescale(m_samplinggrid,RAD); 
    ////// find the first centroid with 6 neighbours
    for(int i=0;i<m_samplinggrid.nvertices();i++){
      if(m_samplinggrid.get_total_neighbours(i)==6){
	m_centroid=i;
	  break;
      }
    }
   
    label_sampling_grid(m_centroid,m_maxs_dist,m_samplinggrid);
    

   
  }

  void SRegDiscreteModel::label_sampling_grid(const int & centroid, const double & dist, NEWMESH::newmesh &Grid)
  {
    m_samples.clear(); m_barycentres.clear(); 
    vector<int> getneighbours;
    vector<int> newneighbours;
    int label;
    Pt centre=Grid.get_coord(centroid);
    Pt sample;
    ColumnVector found(Grid.nvertices()),found_tr(Grid.ntriangles());
   
    found=0;found_tr=0;
    getneighbours.push_back(centroid);
    m_samples.push_back(centre);
    m_barycentres.push_back(centre);
    label=1;
    //// searches for neighbours of the cnetroid that are within the max sampling distance
    while(getneighbours.size()>0){
      for(unsigned int i=0;i< getneighbours.size();i++){
	for(vector<int>::const_iterator j=Grid.nbegin(getneighbours[i]); j!=Grid.nend(getneighbours[i]);j++){
	  sample=Grid.get_coord(*j);
	  if((sample-centre).norm()<=dist && (found(*j+1)==0) && *j!=centroid){
	    m_samples.push_back(Grid.get_coord(*j));  /// pt-centroid equals deformation vector
	    newneighbours.push_back(*j);
	    found(*j+1)=1;
	   
	  }
	}
	for(vector<int>::const_iterator j=Grid.tIDbegin(getneighbours[i]); j!=Grid.tIDend(getneighbours[i]);j++){
	    NEWMESH::Pt v1 = Grid.get_triangle_vertex(*j,0),  v2 = Grid.get_triangle_vertex(*j,1), v3 = Grid.get_triangle_vertex(*j,2);
	    Pt bary((v1.X+v2.X+v3.X)/3,(v1.Y+v2.Y+v3.Y)/3,(v1.Z+v2.Z+v3.Z)/3);
	    bary.normalize(); bary=bary*RAD; 
	    
	    if((bary-centre).norm()<=dist && (bary-centre).norm()> 0 && found_tr(*j+1)==0){
	      for(unsigned int k=0;k<m_barycentres.size();k++){
           	if(abs(1-(((bary-centre)|(m_barycentres[k]-centre))/((bary-centre).norm()*(m_barycentres[k]-centre).norm()))) < 1e-2){ found_tr(*j+1)=1;}
              }
              if(!found_tr(*j+1)){m_barycentres.push_back(bary);   Grid.set_pvalue(Grid.get_triangle_vertexID(*j,0),label);
		label++;}
	      found_tr(*j+1)=1;
	    }
	
	}
      }
    
      getneighbours=newneighbours;
      newneighbours.clear();
    }

  }

  void NonLinearSRegDiscreteModel::estimate_pairs(){
    int pair=0;

    for(int i=0;i<m_CPgrid.nvertices();i++){ /// estimate the total number of edge pairs ////
      m_num_pairs+=m_CPgrid.get_total_neighbours(i);
    }
    m_num_pairs/=2;
    pairs  = new int [m_num_pairs*2];
    for(int i=0;i<m_CPgrid.nvertices();i++){
      for(vector<int>::const_iterator j=m_CPgrid.nbegin(i); j!=m_CPgrid.nend(i);j++){
	if(*j>i){
	  int node_ids[2] = { i, *j};

	  sort_nodeids(node_ids,2);
	  pairs[2*pair]=node_ids[0];
	  pairs[2*pair+1]=node_ids[1];
	  pair++;
	}
      }
   
    }
  }
 
  void NonLinearSRegDiscreteModel::estimate_triplets(){
    m_num_triplets=m_CPgrid.ntriangles();
   
    triplets  = new int [m_num_triplets*3];

    for(int i=0;i<m_CPgrid.ntriangles();i++){
      int node_ids[3] = {m_CPgrid.get_triangle_vertexID(i,0),m_CPgrid.get_triangle_vertexID(i,1), m_CPgrid.get_triangle_vertexID(i,2)};
       sort_nodeids(node_ids,3);

      triplets[3*i]=node_ids[0];
      triplets[3*i+1]=node_ids[1];
      triplets[3*i+2]=node_ids[2];
    }
   
    
  }
    
  void NonLinearSRegDiscreteModel::Initialize(const newmesh &CONTROLGRID){
    
    SRegDiscreteModel::Initialize(CONTROLGRID);

    if(_pairwise){
      //////////////////////// GET PAIRS //////////////
      estimate_pairs(); 
    }
    else{
      /////////////////////////// INITIALIZE TRIPLETS /////////////////////
      estimate_triplets();
    }
   
    /////////////////// CALCULATE FIXED GRID SPACINGS////////
     
    double MVDinput=Calculate_MVD(m_TARGET);  /// estimate spacing of data grid from first vertex
  
    /// INITIALIAZE LABEL GRID/////////
    Initialize_sampling_grid();  
    get_rotations(m_ROT);  /// enables rotation of sampling grid onto every CP

    //////////// INITIALIZE NEIGHBOURHOODS ////////////////

    m_inputrel=boost::shared_ptr<RELATIONS >(new RELATIONS (m_SOURCE,m_TARGET,1.5*asin(MVDinput/RAD)));
    
  }


  
 
  void NonLinearSRegDiscreteModel::get_rotations(vector<Matrix> &ROT){ /// rotates sampling grid to each control point
    Matrix R,R2,R_n,R_n2,R_diff;
    Pt ci,CP,CP_n;
    ROT.clear();
    ci=m_samplinggrid.get_coord(m_centroid);

    for (int k=0;k<m_CPgrid.nvertices();k++){
      CP=m_CPgrid.get_coord(k);
      //// rotate sampling points to overlap with control point
      R=estimate_rotation_matrix(ci,CP);    
      ROT.push_back(R);
    
    }
    
  }

  void NonLinearSRegDiscreteModel::setupCostFunction()
  {
    resetLabeling(); // initialise label array to zero

    ////// use geodesic distances ////////
    costfct->reset_CPgrid(m_CPgrid);

    if(m_iter==1){ 
      costfct->initialize_regulariser();
      m_cp_neighbourhood=boost::shared_ptr<RELATIONS >(new RELATIONS (m_SOURCE,m_CPgrid,3*asin(MVD/RAD))); 
    

      m_cp_neighbourhood->update_RELATIONS(m_SOURCE);
      costfct->set_relations(m_cp_neighbourhood,m_inputrel);
    }

    costfct->reset_anatomical(m_outdir,m_iter);

    get_rotations(m_ROT);  /// instead of recalulating the source->CP neighbours, these are now constant (as source is moving with CPgrid) we we just need to recalculate the rotations of the label grid to the cp vertices
 
   
    if(m_debug){
      char filename[1000];
      sprintf(filename,"%sSOURCE-%d.surf",m_outdir.c_str(),m_iter); m_SOURCE.save(filename);
      sprintf(filename,"%sSOURCE-%d.func",m_outdir.c_str(),m_iter); m_SOURCE.save(filename);
     
      if(m_iter==1){sprintf(filename,"%sTARGET.surf",m_outdir.c_str());  m_TARGET.save(filename); }
      sprintf(filename,"%sCPgrid-%d.surf",m_outdir.c_str(),m_iter); m_CPgrid.save(filename);
      
    }
    //////////// set samples (labels vary from iter to iter ) ///////////////
    if(m_iter%2==0) m_labels=m_samples;
    else m_labels=m_barycentres;

    m_num_labels=m_labels.size();

    costfct->set_labels(m_labels,m_ROT);
    if(m_verbosity)cout << " initialize cost function " << m_iter <<  endl;
    costfct->initialize(m_num_nodes,m_num_labels,m_num_pairs,m_num_triplets);
   
    costfct->get_source_data();

   
    if(_pairwise) costfct->setPairs(pairs);
    else costfct->setTriplets(triplets);

      
   
    m_iter++;
  }


  void NonLinearSRegDiscreteModel::applyLabeling(int *dlabels){
   
    if(dlabels){
      for (int i=0;i<m_CPgrid.nvertices();i++){
	//// rotate sampling points to overlap with control point
	//// transform grid point to new position given by label
	Pt newpt=m_ROT[i]*m_labels[dlabels[i]];

	m_CPgrid.set_coord(i,newpt);
    }

    }
   
  } 

  void MetricDistortionDiscreteModel::Initialize(const newmesh &CONTROLGRID){
    m_CPgrid=CONTROLGRID;
 
    m_num_nodes=m_CPgrid.nvertices();
    m_maxs_dist=0.4*Calculate_MaxVD(m_CPgrid);
    initLabeling();
    /////////////////////////// INITIALIZE TRIPLETS /////////////////////
    estimate_triplets();
    Initialize_sampling_grid();  
    get_rotations(m_ROT);  /// enables rotation of sampling grid onto every CP
    m_iter=1;
   
  }

  void MetricDistortionDiscreteModel::setupCostFunction()
  {
    resetLabeling(); // initialise label array to zero

    ////// use geodesic distances ////////
    costfct->reset_CPgrid(m_CPgrid);

    get_rotations(m_ROT);  /// instead of recalulating the source->CP neighbours, these are now constant (as source is moving with CPgrid) we we just need to recalculate the rotations of the label grid to the cp vertices
 
    //////////// set samples (labels vary from iter to iter ) ///////////////
    if(m_iter%2==0) m_labels=m_samples;
    else m_labels=m_barycentres;

    m_num_labels=m_labels.size();
   


    costfct->set_labels(m_labels,m_ROT);
    if(m_verbosity)cout << " initialize cost function " << m_iter << endl;
    costfct->initialize(m_num_nodes,m_num_labels,m_num_pairs,m_num_triplets);
    
    
    m_iter++;
  }

  void GroupSegmentationModel::Initialize(){
   

    if (m_num_subjects ==0){throw  DISCRETEOPTHOCRException("Group segmentation Model:: DATA has not been set");}

    m_num_nodes=m_SOURCE.nvertices()*m_num_subjects;
    ///////// INIT LABELLING
    m_num_labels=2;
    initLabeling();

    // GET TRIPLETS: for intra subject edges ////////
    m_num_triplets=m_num_subjects*m_SOURCE.ntriangles();
    triplets  = new int [m_num_triplets*3];
 
    int triplet=0;
    for(int num=0;num<m_num_subjects;num++){
      for(int i=0;i<m_SOURCE.ntriangles();i++){
	if(triplet>=m_num_triplets){ cout << i << " too many triplets " << triplet <<  "num " << num << endl; exit(1);}
	triplets[3*triplet]=m_SOURCE.nvertices()*num + m_SOURCE.get_triangle_vertexID(i,0);
	triplets[3*triplet+1]=m_SOURCE.nvertices()*num + m_SOURCE.get_triangle_vertexID(i,1);
	triplets[3*triplet+2]=m_SOURCE.nvertices()*num + m_SOURCE.get_triangle_vertexID(i,2);
	
	triplet++;
      }
    }

    /////////////////// CALCULATE FIXED GRID SPACINGS////////
     
    double MVDmax= Calculate_MaxVD(m_SOURCE);  /// estimate spacing of data grid from first vertex
  
   
    //////////// INITIALIZE NEIGHBOURHOODS ////////////////

    m_neighbourhood=boost::shared_ptr<RELATIONS >(new RELATIONS (m_SOURCE,m_SOURCE,3*asin(MVDmax/RAD))); 
    m_neighbourhood->update_RELATIONS(m_SOURCE);
    
    //////////////////////// GET PAIRS //////////////
    int pair=0;
    m_num_pairs=(m_num_subjects-1)*m_count*m_SOURCE.nvertices();  // iterates through ALL intra subjects pairs (for all subjects) and then ALL inter subject pairs

    pairs  = new int [m_num_pairs*2];
   
    for (int i=1;i<=m_SOURCE.nvertices();i++){
      for (int j=1;j<=m_count;j++){
	if(m_neighbourhood->Nrows(i)<m_count)  {throw  DISCRETEOPTHOCRException("Group segmentation Model:: m_neighbourhood is not big enough");}
	  for(int num=1;num<m_num_subjects;num++){
	   
	    if(pair>=m_num_pairs){ cout << i << " too many pairs " << pair <<  " m_count " << m_count << endl; exit(1);}
	    pairs[2*pair]=(i-1);
	    pairs[2*pair+1]=m_SOURCE.nvertices()*num+((*m_neighbourhood)(j,i)-1);
	    pair++;
	  }
      }
    }

   
    
  }

  void GroupSegmentationModel::setupCostFunction()
  {
    resetLabeling(); // initialise label array to zero

    costfct->set_relations(m_neighbourhood);

    if(m_verbosity)cout << " initialize segmentation cost function " << endl;
    costfct->initialize(m_num_nodes,m_num_labels,m_num_pairs,m_num_triplets);

      
   if(m_verbosity) cout << " numpoints " << m_num_nodes << " m_num_labels " << m_num_labels << " m_num_pairs " << m_num_pairs <<  endl;
  }

  void GroupSegmentationModel::applyLabeling(int *dlabels){
    
    SEGMENTATION.ReSize(m_num_subjects,m_SOURCE.nvertices());
    if(dlabels){
      for(int num=0;num< m_num_subjects;num++){
	for (int i=0;i<m_SOURCE.nvertices();i++){
	  SEGMENTATION(num+1,i+1)=dlabels[num*m_SOURCE.nvertices()+i];
	}
      }
      
    } 
  }
}
