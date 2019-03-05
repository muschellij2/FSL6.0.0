/*  calculatevertexwisecorr.cc

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
#include <time.h>
#include <stdio.h>
#include "newmesh/meshfns.h"
#include "MeshReg/meshreg.h"

using namespace std;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "calculatevertexwisecorr  <surface1> <surface2> <outbase>  " << endl;
  cout << "Input are func files assumes surfaces have been resampled to same mesh   " << endl;
  cout << " options: " << endl;
  cout << "-compare X compare correlation to this file " << endl;
}

int main( int argc, char **argv )
{
  newmesh SURFACE1, SURFACE2, OUTSURF,COMP;
  vector<double> data1,data2;
  string outputname;
  sparsesimkernel<double> sim;
  double correlation,sum=0.0;
  bool _compare=false;
  int ok=0;

  if(argc < 3){

    Usage();
    exit(0);
  }

  SURFACE1.load(argv[1],false,false);
  argc--; 
  argv++;
  SURFACE2.load(argv[1],false,false);
  argc--; 
  argv++;
  outputname=argv[1];
  argc--; 
  argv++;

  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-compare") == 0)){
      argc--;
      argv++;
      _compare=true;
      COMP.load(argv[1],false,true);
      argc--; 
      argv++;
      ok = 1;
    }
  }
  OUTSURF=SURFACE1;

  OUTSURF.clear_data();
  ColumnVector allcorr(SURFACE1.npvalues()); 

  if(SURFACE1.npvalues()!=SURFACE2.npvalues() || SURFACE1.get_dimension()!=SURFACE2.get_dimension()) { cout << " surfaces have different number of vertices, or features abort! " << endl; exit(1);}

  if(SURFACE1.get_dimension()==1){
	for (int i=0;i<SURFACE1.npvalues();i++){	
		data1.push_back(SURFACE1.get_pvalue(i));
		data2.push_back(SURFACE2.get_pvalue(i));
      
	}
	correlation=sim.corr(data1,data2);
	cout << " mean=" << correlation << endl;
  }
  else{
	for (int i=0;i<SURFACE1.npvalues();i++){
		data1.clear();data2.clear();
    
		for (int j=0;j<SURFACE1.get_dimension();j++){
		data1.push_back(SURFACE1.get_pvalue(i,j));
		data2.push_back(SURFACE2.get_pvalue(i,j));
      
		}
		correlation=sim.corr(data1,data2);
		allcorr(i+1)= correlation;
    
		sum+=correlation;
	}
  
  OUTSURF.set_pvalues(allcorr);
  double mean,var=0.0;
  mean=sum/SURFACE1.npvalues();
  for (int i=0;i<SURFACE1.npvalues();i++)
    var+=(OUTSURF.get_pvalue(i) -mean)*(OUTSURF.get_pvalue(i) -mean);

  var/=SURFACE1.npvalues();

  cout << " mean=" << mean << endl;
  cout << " var=" << var << endl;

  OUTSURF.save(outputname+".func.gii");
 }
  if(_compare){
    for (int i=0;i<SURFACE1.npvalues();i++){
      double diff=OUTSURF.get_pvalue(i)-COMP.get_pvalue(i);
      OUTSURF.set_pvalue(i,diff);
    }
  OUTSURF.save(outputname+"DIFF.func.gii");

  }

}
