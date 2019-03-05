/*  msm_metricmath.cc

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
/* this program is designed to downsample freesurfer label files to be used in combination with the SPH6.vtk or other downsampled meshes*/

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <sstream>
#include <limits>
#include "utils/options.h"
#include "newmesh/meshfns.h"
#include "MeshReg/meshreg.h"


using namespace std;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

enum sValue { evND, 
	      ev1, 
	      ev2, 
	      ev3, 
	      ev4, 
	      ev5,
	      ev6,
	      evfin };

  // Map to associate the strings with the enum values
static std::map<std::string, sValue> s_mapValues;
 
void set_method(const string &M){
    s_mapValues["-mean"] = ev1;
    s_mapValues["-variance"] = ev2;
    s_mapValues["-sum"] = ev3;
    s_mapValues["-max"] = ev4;
    s_mapValues["-min"] = ev5;  
    s_mapValues["-rank"] = ev6;      
    s_mapValues["-percentile"] = evfin;    

  }
  
void Usage()
{
  cout << "msm_metricmath  <inputmetric> -operations <operations> -options <options>   " << endl;
  cout << "options:  " << endl;
  cout << "-exclude lthr upthr use exclusion mesh for cut defined using upper (uthr) and lower (lthr) thresholds" << endl;
  cout << "-vertex-wise estimate statistics across columns (not implemented yet)" << endl;
  cout << " -exclusionmask use another input file to define exclusion " << endl;
  cout << " -abs estimate absolute average " << endl;

  cout << "operations are one or more of:  " << endl;
  cout << "-mean   " << endl;
  cout << "-variance  " << endl;  
  cout << "-sum  " << endl;  
  cout << "-max  " << endl;  
  cout << "-min  " << endl;  
  cout << "-percentile X " << endl;  
}


vector<float> estimate_average(vector<float> perc, bool makeabs,newmesh &IN, boost::shared_ptr<NEWMESH::newmesh> EXCL){
  
  vector<float> average(IN.get_dimension(),0);
  vector<int> totalcounts(IN.get_dimension(),0.0);
  cout << " EXCL.get() " << EXCL.get() <<  endl;
  
  for (int d=0;d<IN.get_dimension();d++){
    for(int i=0;i<IN.npvalues();i++){
      //cout << i << " (!EXCL.get() | EXCL->get_pvalue(i,d)) " << (!EXCL.get() || EXCL->get_pvalue(i,d)) << " IN.get_pvalue(i,d)) " << IN.get_pvalue(i,d) << " (abs(IN.get_pvalue(i,d)) > EPSILON) " << (abs(IN.get_pvalue(i,d)) > EPSILON) << endl;
      if((!EXCL.get() || EXCL->get_pvalue(i)) && (abs(IN.get_pvalue(i,d)) > EPSILON) && (abs(IN.get_pvalue(i,d)) < perc[d])){
	//EXCL->set_pvalue(i,10);
	if(makeabs) average[d]+= abs(IN.get_pvalue(i,d));
	else average[d]+= IN.get_pvalue(i,d);
	totalcounts[d]+=1;
      }
    }		
    average[d]/=totalcounts[d];
    if(abs(average[d]) < 0.001) average[d]=0;
    cout  << "average of dimension " << d << " is " << average[d] <<endl;
    //    if(abs(average[d]) < 0.001) average[d]=0;
 
  }	
  // EXCL->save("EXCL2.func.gii");
  return average;
}

void estimate_pca(newmesh &IN){
  
   
  Matrix DATA=IN.get_pvalues();
  DiagonalMatrix eigenvals;
  double sumvar=0;
  SVD(DATA,eigenvals);

  vector<float> vals(eigenvals.Nrows(),0);
 
  double tolerance = Max(DATA.Nrows(),DATA.Ncols()) * eigenvals.Maximum() * 1e-16;

  int therank=0;

  for(int i=0; i<eigenvals.Nrows(); i++){
    sumvar+=eigenvals(i+1);
    cout  << "eigval  " << i << " is " <<eigenvals(i+1) <<endl;
    if (eigenvals(i+1)>tolerance)
      therank++;
  
  }

  for(int i=0; i<eigenvals.Nrows(); i++)
    cout  << " percentage variance " << i << " is " <<eigenvals(i+1)/sumvar <<endl;
  cout << " rank " << therank << endl;	
  // EXCL->save("EXCL2.func.gii");
}

void estimate_variance(bool  makeabs,newmesh &IN, boost::shared_ptr<NEWMESH::newmesh> EXCL, vector<float> average = vector<float>()){
  
  vector<float> thresh(IN.get_dimension(),1);  
  if(average.size()!=IN.get_dimension()) average=estimate_average(thresh,makeabs,IN,EXCL);
  vector<float> variance(IN.get_dimension(),0);
  vector<int> totalcounts(IN.get_dimension(),0.0);
  
  for (int d=0;d<IN.get_dimension();d++){
    for(int i=0;i<IN.npvalues();i++){
      if((!EXCL.get() || EXCL->get_pvalue(i)) && (abs(IN.get_pvalue(i,d)) > EPSILON)){
	if(makeabs) variance[d]+=( abs(IN.get_pvalue(i,d))-average[d])*(abs(IN.get_pvalue(i,d))-average[d]);
	else variance[d]+=( IN.get_pvalue(i,d)-average[d])*( IN.get_pvalue(i,d)-average[d]);
	totalcounts[d]+=1;
      }
    }			
    cout  << "variance of dimension " << d << " is " << variance[d]/totalcounts[d] << " average " <<  average[d] << endl;
		
  }	
  
  
}


void estimate_sum(bool makeabs,  newmesh &IN, boost::shared_ptr<NEWMESH::newmesh> EXCL){
    

  vector<float> sum(IN.get_dimension(),0);
  
  for (int d=0;d<IN.get_dimension();d++){
    for(int i=0;i<IN.npvalues();i++){
      if((!EXCL.get() || EXCL->get_pvalue(i)) && abs(IN.get_pvalue(i,d)) > EPSILON){
	if( makeabs) sum[d]+= abs(IN.get_pvalue(i,d));
	else sum[d]+= IN.get_pvalue(i,d);
	
      }
    }		
    cout  << "sum of dimension " << d << " is " << sum[d] << endl;
    
  }	
  
 
}

void estimate_max(bool  makeabs,newmesh &IN,  boost::shared_ptr<NEWMESH::newmesh> EXCL){
  
  

  vector<float> max(IN.get_dimension(),0);
  double val;
  for (int d=0;d<IN.get_dimension();d++){
    for(int i=0;i<IN.npvalues();i++){
      if((!EXCL.get() || EXCL->get_pvalue(i)) && abs(IN.get_pvalue(i,d)) > EPSILON){
	if( makeabs) val=abs(IN.get_pvalue(i,d));
	else val=IN.get_pvalue(i,d);
	if(val > max[d])
	  max[d]=val;
      }
    }		
    cout  << "max of dimension " << d << " is " << max[d] <<endl;
    
  }	
  
 
}

vector<float> estimate_min(bool  makeabs, newmesh &IN,  boost::shared_ptr<NEWMESH::newmesh> EXCL){
  
  

  vector<float> min(IN.get_dimension(),std::numeric_limits<float>::max());
  double val;
  for (int d=0;d<IN.get_dimension();d++){
    for(int i=0;i<IN.npvalues();i++){
      if((!EXCL.get() || EXCL->get_pvalue(i)) && abs(IN.get_pvalue(i,d)) > EPSILON){
	if( makeabs) val=abs(IN.get_pvalue(i,d));
	else val=IN.get_pvalue(i,d);
	if(val < min[d])
	  min[d]=val;
      }
    }		
    cout  << "min of dimension " << d << " is " << min[d] <<endl;
    
  }	
  
  return min;
}
vector<float> estimate_percentile(float perc, int bins, bool  makeabs, newmesh &IN, boost::shared_ptr<NEWMESH::newmesh> EXCL){
  
  ColumnVector DATA(IN.npvalues());

  vector<float> percentile(IN.get_dimension(),0);
  cout << " makeabs" << makeabs << endl;
  double val;
  for (int d=0;d<IN.get_dimension();d++){
	DATA=0;
    for(int i=0;i<IN.npvalues();i++){
    
	if( makeabs) DATA(i+1)=abs(IN.get_pvalue(i,d));
	else DATA(i+1)=IN.get_pvalue(i,d);
	//if(DATA(i+1)<0) cout << i << " DATA " << DATA(i+1) << endl;
    }	
    Histogram hist=build_histogram(DATA,EXCL,bins);
  //  hist.generateCDF(); 
  //  ColumnVector CDF=hist.getCDF();
    percentile[d]=hist.getPercentile(perc);	
    cout  << "percentile " << perc << "of dimension " << d << " is " << percentile[d] <<endl;
    
  }	
  
  return percentile;
}

int main(int argc, char **argv){


  boost::shared_ptr<NEWMESH::newmesh> EXCL;
  vector<float> thresh;
  NEWMESH::newmesh IN,exclusionmask;
  int ok,_bins=256;
  bool _exclude=false,_vertexwise=false;
  bool _abs=false;
  boost::shared_ptr<BFMatrix> DATA;
  string output;
  vector<string> operations;
  double uthr=1,lthr=0;
  char filename[1000];
  float _perc;
  float limit=0.99999;
  if(argc < 2){

    Usage();
    exit(0);
  }

  cout << argv[1] << endl;
  set_data(argv[1],DATA,IN);
  argc--; 
  argv++;

  
  cout << " here " << endl;
//operations
  if(strcmp(argv[1], "-operations")==0) {
    argc--;
    argv++;
    while(argc > 1 && strcmp(argv[1], "-options") != 0){
      operations.push_back(argv[1]);
      if(strcmp(argv[1], "-percentile")== 0){
		  cout << " 1 " <<  argv[1] << endl;
	  argc--;
      argv++;  
      cout << argv[1] << endl;
	  _perc=atof(argv[1]);
      argc--;
      argv++;
      _bins=atoi(argv[1]);
       argc--;
      argv++;
	  }else{
		argc--;
		argv++;
	}
    }
  }
  
 //options
  while (argc > 1){
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-options") == 0)){
      argc--;
      argv++;    
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-vertex-wise") == 0)){
      argc--;
      argv++;
      _vertexwise=true;    
      ok = 1;
      cout << " not implemented yet " <<endl; exit(1);
    }
    else if((ok == 0) && (strcmp(argv[1], "-exclude") == 0)){
      argc--;
      argv++;
      _exclude=true;
      lthr=atof(argv[1]);
      argc--;
      argv++;
      uthr=atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-exclusionmask") == 0)){
      argc--;
      argv++;
      set_data(argv[1],DATA,exclusionmask);
       cout << exclusionmask.npvalues() << " " << argv[1] <<  endl;
        _exclude=true;
      argc--;
      argv++;
      ok = 1;
    }
     else if((ok == 0) && (strcmp(argv[1], "-abs") == 0)){
      argc--;
      argv++;
      cout << " abs " << endl;
      _abs=true;
    
      ok = 1;
    }
     else if((ok == 0) && (strcmp(argv[1], "-limit") == 0)){
       argc--;
       argv++;
       cout << " abs " << endl;
       limit=atof(argv[1]);
       argc--;
       argv++;

       ok = 1;
     }

     else{cout << argv[1] << " option doesn't exist " << endl; exit(1);}
  }

  if(_exclude){
    if(exclusionmask.npvalues()>0){
      EXCL= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(exclusionmask,exclusionmask.get_pvalues(),lthr,uthr)));
      cout << exclusionmask.npvalues() << endl;
      
      EXCL->save("exclusion.func"); exclusionmask.save("exclusion1.func"); 
      for (int i=0;i<EXCL->npvalues();i++){
	if(EXCL->get_pvalue(i)==0) EXCL->set_pvalue(i,1);
	else EXCL->set_pvalue(i,0);

      }
    }
    else
      EXCL= boost::shared_ptr<NEWMESH::newmesh>(new NEWMESH::newmesh(create_exclusion(IN,IN.get_pvalues(),lthr,uthr)));
  }

  cout << " limit " << limit << endl;
  thresh=estimate_percentile(limit,100, true,IN,EXCL);

  
 // EXCL->save("EXCL.func");
  vector<float> av;
  ////////////////////////////
  for (unsigned int op=0;op<operations.size();op++){
    set_method(operations[op]);

    
    switch(s_mapValues[operations[op]])
      {
      case ev1:
	av=estimate_average(thresh,_abs,IN,EXCL);
	break;
      case ev2:
	estimate_variance(_abs,IN,EXCL,av);
	break;
      case ev3:
	estimate_sum(_abs,IN,EXCL);
	break;
      case ev4:
	estimate_max(_abs,IN,EXCL);
	break;
      case ev5:
	estimate_min(_abs,IN,EXCL);
	break; 
      case ev6:
	estimate_pca(IN);
	break;
      case evfin:
	estimate_percentile(_perc,_bins,_abs,IN,EXCL);            
	break;
      default:
	cout << "'" <<operations[op] 
	     << "' is an invalid string or not yet implemented " << endl;
	break;
      }
  }
	 
}
