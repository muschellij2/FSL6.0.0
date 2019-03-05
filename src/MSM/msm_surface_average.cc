/*  msm_surface_average.cc

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
#include "utils/options.h"
#include "newmesh/meshfns.h"
#include "MeshReg/meshreg.h"


using namespace std;
using namespace MISCMATHS;
using namespace NEWMESH;
using namespace MESHREG;

void Usage()
{
  cout << "msm_surface_average  <asciii list of data> <output> <-options>   " << endl;
  cout << "options:  " << endl;
  cout << "-weighted X supply a list of distances from the mean" << endl;
  cout << "-sigma X  sigma values " << endl;
}


int main(int argc, char **argv){

 
  NEWMESH::newmesh IN, AVERAGE;
  
  string output;
  vector<string> INlist,distances;
  bool _weighted=false;
  double _sigma;
  double weight;
  double sumweighted=0;
  int ok=0;
 

  if(argc < 2){

    Usage();
    exit(0);
  }

  INlist=read_ascii_list(argv[1]);
  argc--; 
  argv++;
  output=argv[1];
  argc--; 
  argv++;
  
  

  while (argc > 1) {
    ok = 0;
    if((ok == 0) && (strcmp(argv[1], "-weighted") == 0)){
      argc--;
      argv++;
      _weighted=true;
      distances=read_ascii_list(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else if((ok == 0) && (strcmp(argv[1], "-sigma") == 0)){
      argc--;
      argv++;
      _sigma=atof(argv[1]);
      argc--;
      argv++;
      ok = 1;
    }
    else{cout << argv[1] << " option doesn't exist " << endl; exit(1);}
  }

  AVERAGE.load(INlist[0]);  

  for(int i=0;i<AVERAGE.nvertices();i++){
    Pt p;
    AVERAGE.set_coord(i,p);
  }
  double mean=0;
  //for (unsigned int i=0;i<INlist.size();i++){
	//mean+=atof(distances[i].c_str());
//}
///mean/=INlist.size();
  /////////////////////////////////////////////////////
  for (unsigned int i=0;i<INlist.size();i++){
    cout << i << " " << INlist[i] << endl;

    IN.load(INlist[i]);
  //  cout << " recentre " << endl;
   // recentre(IN);
    if(_weighted){
      if(_sigma==0){ cout << "must set sigma for weighted" <<endl; exit(1);}
      weight=1/(_sigma*(sqrt(2*PI)))*exp(-0.5*(((atof(distances[i].c_str())-mean))/_sigma)*(((atof(distances[i].c_str())-mean))/_sigma));
      sumweighted+=weight;
      cout << "distance " <<  distances[i] << " weight " << weight << endl;
    }else{ weight=1;sumweighted+=1;}

    
    for( int n=0;n<IN.nvertices();n++){
      Pt p=AVERAGE.get_coord(n);
      Pt p2=IN.get_coord(n);
    
      Pt newpt;
      newpt.X=p.X+weight*p2.X;      
      newpt.Y=p.Y+weight*p2.Y;
      newpt.Z=p.Z+weight*p2.Z;
      AVERAGE.set_coord(n,newpt);
    }

  }

  for( int n=0;n<AVERAGE.nvertices();n++){
      Pt p=AVERAGE.get_coord(n);
      p.X=p.X/sumweighted;
      p.Y=p.Y/sumweighted;
      p.Z=p.Z/sumweighted;
      AVERAGE.set_coord(n,p);
  }

  //recentre(AVERAGE);
  AVERAGE.save(output+"AVERAGE.surf");
 
}
