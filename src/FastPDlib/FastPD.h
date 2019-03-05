//#############################################################################
//#
//# FastPD Optimization (c) 2008.
//#
//# Original Author: Nikos Komodakis
//# Reimplementation: Ben Glocker
//#
//# This is an implementation of the FastPD optimizer as described in:
//# N. Komodakis, G. Tziritas, N. Paragios,
//# Fast, Approximately Optimal Solutions for Single and Dynamic MRFs,
//# IEEE Conference on Computer Vision and Pattern Recognition (CVPR), 2007.
//#
//# THE WORK IS ONLY FOR RESEARCH AND NON-COMMERCIAL PURPOSES. THE OPTIMIZATION
//# CODE IS PROTECTED FROM SEVERAL US/EU/CHINA PENDING PATENT APPLICATIONS.
//# IF YOU WOULD LIKE TO USE THIS SOFTWARE FOR COMMERCIAL PURPOSES OR LICENSING
//# THE TECHNOLOGY, PLEASE CONTACT:
//# ECOLE CENTRALE DE PARIS, PROF. NIKOS PARAGIOS (nikos.paragios@ecp.fr).
//#
//# If you intend to use this code or results obtained with it, the above
//# mentioned paper should be cited within your publication.
//#
//#############################################################################

#ifndef _FASTPDOPTIMIZATION_
#define _FASTPDOPTIMIZATION_

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <vector>
#include "graph.h"

#if !defined(__DiscreteModel_h)
#define __DiscreteModel_h
#include "DiscreteOpt/DiscreteModel.h"
#endif

using namespace DISCRETEOPT;


namespace FPD{

  class FastPD
  {
  public:
    typedef Graph::Real Real;
    typedef short TIME;
    
		// Auxiliary data structures and variables
    struct Node_info
    {
      Graph::Label label; 
      Real height; 
      TIME time;
      int next;    
      int prev;
      int *pairs; 
      int numpairs;
    };

    struct Pair_info
    {
			int i0, i1;
      TIME time;
    };
    
    struct Arc_info
    {
      int head, tail;
      Real balance;
    };
    
    FastPD( boost::shared_ptr<DiscreteModel> discreteModel, int max_iterations, bool copy_unary = false );
    
    ~FastPD();
    
    double run();
    
    double getInitialEnergy() { return m_initial_energy; }
    
    void getLabeling(int *labeling);

    void setDisplayFunc(void (*display)(void *), void *inst) { m_display = display;  m_displayinst = inst; }
    
    private:
    void init_duals_primals();
    
    void inner_iteration( Graph::Label label );
    
    void inner_iteration_adapted( Graph::Label label );
		
    void track_source_linked_nodes( Graph::Label label );        
    
    void fillGraph( Graph *_graph );
    
    void createNeighbors();
    
    boost::shared_ptr<DiscreteModel>	model;				///< Discrete model.
    boost::shared_ptr<DiscreteCostFunction>	costfct;			///< Cost function to be minimized.
    int				m_num_nodes;		///< Number of MRF nodes.
    int				m_num_labels;		///< Number of labels.
    int				m_num_pairs;		///< Number of node pairs (i.e. MRF edges).
    int				m_max_iterations;	///< Number of maximum iterations.
    int				*pairs;				///< Node pairs look up table.
    int				m_iterations;		///< Iterations counter.
    
    Real			*height;			///< Height variables.
    Real			*balance;			///< Balance variables.
    bool			m_delete_height;	///< Indicates whether the height variables have to be deleted in destructor.
    
    Graph			**graphs;			///< List of all graph instances.
    Graph::node		*graph_nodes;		///< List of graph node instances.
    Graph::arc		*graph_edges;		///< List of graph edge instances.
    Graph::node		**children;			///< List of all node children instances.
    Node_info		*node_info;			///< List of node information.
    Arc_info		*edge_info;			///< List of edge information.
    Pair_info		*pair_info;			///< List of pair information.
    int				*pairs_arr;
    
    double			 m_initial_energy;	///< Initial MRF energy.
    double			 m_energy;			///< MRF energy.
    int				 m_time;
    int				 m_active_list;
    int				 m_energy_change_time;
    
    int				*source_nodes_tmp1; 
    int				*source_nodes_tmp2; 			
    
    void			(*m_display)(void *);		///< Callback function.
    void			*m_displayinst;				///< Callback object instance.
    
    // Assignment or copying are not allowed
    FastPD(const FastPD &other);
    FastPD operator=( const FastPD &other );
    
    static void err_fun(char * msg)
    {
      printf("%s",msg);
    }
  };
} 
#endif
  
