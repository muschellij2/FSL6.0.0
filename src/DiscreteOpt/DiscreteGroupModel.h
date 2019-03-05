/*  DiscreteGroupModel.h

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

#pragma once

#include "DiscreteGroupCostFunction.h"
#include "DiscreteModel.h"
#include <vector>

 

namespace DISCRETEOPTHOCR{
 
 

  class GroupDiscreteModel: public SRegDiscreteModel
  {     
  protected:
    vector<NEWMESH::newmesh> m_DATAMESHES; // TARGET MESH
    vector<NEWMESH::newmesh> m_CONTROLMESHES; // TARGET MESH
   
    NEWMESH::newmesh m_TEMPLATE;   /// template data grid
    NEWMESH::newmesh m_TEMPLATE_LR; /// template control grid

    /////////////////////// NEIGHBOURHOOD INFO ////////////////////
    boost::shared_ptr<RELATIONS>  m_TEMPLATE_LR_ALL_RELATIONS; /// hold neighbourhood between control grids and to LR target
    vector<vector<vector<int> > >  BETWEEN_SUBJECT_PAIRS;
    vector<int> subjects;   
  
    int m_num_subjects;
    int control_grid_size;
  public:
    /**
     * Constructor.
     */
    GroupDiscreteModel(){};  
    ~GroupDiscreteModel(){};

    GroupDiscreteModel(myparam & P){  
      m_CPres=2; m_SGres=4; m_simmeasure=2; m_multivariate=false; m_verbosity=false; m_outdir="";m_debug=false;  
      set_parameters(P);
      cout << " create cf group 1 " << endl;
      costfct=boost::shared_ptr<GroupDiscreteCostFunction>(new GroupDiscreteCostFunction());
      m_inputrel=boost::shared_ptr<RELATIONS >(new RELATIONS());
      m_cp_neighbourhood=boost::shared_ptr<RELATIONS >(new RELATIONS ()); 
      costfct->set_parameters(P);
     
    };

    void set_meshspace(const NEWMESH::newmesh & target,const NEWMESH::newmesh & source, const int num=1){  m_TEMPLATE=target;
      m_DATAMESHES.clear();
      m_num_subjects=num;
      for (int i=0;i<num;i++)
	m_DATAMESHES.push_back(source);
    }

    void reset_meshspace(const NEWMESH::newmesh & source, int num=0){  
      char filename[1000];
      cout << " reset meshspace" << endl;
      m_DATAMESHES[num]=source; 
      //sprintf(filename,"SOURCE_set_up_it%d_res%d_num%d.surf.gii",m_iter, m_CPres,num);  m_DATAMESHES[num].save(filename);
      if(costfct.get()){costfct->reset_source(source,num);}
      else{throw  DISCRETEOPTHOCRException("Discrete Model:: You cannot reset the source mesh withou initialising the discrete costfunction");} 
    }
    
    void reset_CPgrid(const NEWMESH::newmesh & grid, int num=0){  
       m_CONTROLMESHES[num]=grid;
    }
  
    void warp_CPgrid(NEWMESH::newmesh & START, NEWMESH::newmesh & END, int num=0){ barycentric_mesh_interpolation(m_CONTROLMESHES[num],START,END); unfold(m_CONTROLMESHES[num]);	}


    void applyLabeling(){applyLabeling(labeling);} ;
    void applyLabeling(int *discreteLabeling);

    void initialize_quartets();
    void initialize_pairs();

    void estimate_pairs();
    void estimate_triplets();
    void estimate_quartets();
    void estimate_combinations(int, int*);
 

    void Initialize(const newmesh &);
    void Initialize(){};
    void get_rotations(vector<Matrix>  &);
    void setupCostFunction();

    //// functions for keeping track of between mesh relationships
    void get_between_subject_pairs();
    void resample_to_template();

    NEWMESH::newmesh get_CPgrid(int num=0){return m_CONTROLMESHES[num];}

  };


}
