/*  msmGOptions.h

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2008 University of Oxford  

    Some sections of code inspired by A. Petrovic.
*/
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

#if !defined(msmGOptions_h)
#define msmGOptions_h

#include <string> 
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;
using namespace NEWMESH;

class msmGOptions {

public:
  static msmGOptions& getInstance();
  ~msmGOptions() { delete gopt; }
  
 
  Option<bool>   help;
  Option<bool>   verbose;
  Option<bool>   printoptions;
  Option<bool> version;
  Option<bool> debug;

  Option<string> meshes; 
  Option<string> templatemesh; 
  Option<string> data;
 
  Option<string> outbase; 
  Option<string> parameters;  /// slowly replace most of these with the parameter file??
  // Option<string> L1matlabpath;

  Option<int> multiresolutionlevels;
  Option<float>  smoothoutput;

  bool parse_command_line(int argc, char** argv,Log& logger);
  
  vector<int> return_datarange(NEWMESH::newmesh);
 private:
  msmGOptions();  
  const msmGOptions& operator=(msmGOptions&);
  msmGOptions(msmGOptions&);

  OptionParser options;
      
  static msmGOptions* gopt;
  
};

 inline msmGOptions& msmGOptions::getInstance(){
 
   if(gopt == NULL)
     gopt = new msmGOptions();
  
  return *gopt;
 }

inline msmGOptions::msmGOptions() :
		  help(string("-h,--help"), false,
		       string("display this message"),
		       false, no_argument),
		  verbose(string("-v,--verbose"), false,
		       string("switch on diagnostic messages"),
		       false, no_argument),
		  printoptions(string("-p,--printoptions"), false,
		       string("print configuration file options"),
		       false, no_argument),
		  version(string("--version"), false,
		       string("report version informations"),
		       false, no_argument),
		  debug(string("--debug"), false,
		       string("run debugging or optimising options"),
		       false, no_argument),
		  meshes(string("--meshes"), string(""),
			    string("list of paths to input meshes (available formats: VTK, ASCII, GIFTI). Needs to be a sphere"),
			    true , requires_argument),
		  templatemesh(string("--template"), string(""),
			    string("templates sphere for resampling (available formats: VTK, ASCII, GIFTI). Needs to be a sphere"),
			    true , requires_argument),
		  data(string("--data"), string(""),
				string("list of paths to the data"),
				false , requires_argument),
		  outbase(string("-o,--out"), string(""),
			  string("output basename"),
			  true, requires_argument),
		  parameters(string("--conf"), string(""),
			     string("\tconfiguration file "),
			     false, requires_argument),
		   //  L1matlabpath(string("--L1path"), string(""),
		   //	     string("\tPath to L1min matlab code"),
		   //		       false, requires_argument),
		  multiresolutionlevels(string("--levels"),0, 
					string("number of resolution levels (default = number of resolution levels specified by --opt in config file)"), 
					false, requires_argument),
		  smoothoutput(string("--smoothout"), 0, 
			       string("smooth tranformed output with this sigma (default=0)"), 
			       false, requires_argument),
		  
		  options("msm", "msm [options]\n")
{
  try {
   
        options.add(help); 
	options.add(verbose);
	options.add(printoptions);
	options.add(version);
	options.add(debug);
        options.add(meshes);
        options.add(templatemesh);
        options.add(data);
	options.add(outbase);
	options.add(parameters);
	//	options.add(L1matlabpath);
	options.add(multiresolutionlevels);
	options.add(smoothoutput);
     }
     catch(X_OptionError& e) {
       options.usage();
       cerr << endl << e.what() << endl;
     } 
     catch(std::exception &e) {
       cerr << e.what() << endl;
     }    
     

     
}

inline bool msmGOptions::parse_command_line(int argc, char** argv,Log& logger){
	 
	 
  for(int a = options.parse_command_line(argc, argv); a < argc; a++) ;
  if(!(printoptions.value() || version.value()) ){
  if(help.value() || ! options.check_compulsory_arguments())
    {
      options.usage();
      //throw NEWMAT::Exception("Not all of the compulsory arguments have been provided");
      exit(2);
    }      
  
  logger.makeDir(outbase.value()+"logdir","MSM.log");
	  
 
  
  // do again so that options are logged
  for(int a = 0; a < argc; a++)
    logger.str() << argv[a] << " ";

  logger.str() << endl << "---------------------------------------------" << endl << endl;
  }  
  return true;
}



/* inline vector<int> msmOptions::return_datarange(NEWMESH::newmesh in){ */

/*   vector<int> V; */
/*   int range2; */
/*   int i; */
/*   int indexval=index.value(); */
 
/*   if(range.value()<1E-7){ */
/*     range2=in.nvertices(); */
/*     cout << " range2 " << range2 << endl;; */
/*   } */
/*   else{ */
/*     range2=range.value(); */
/*     cout << " here2 " << endl; */
/*   } */
/*   if(patch.value()==""){ */
/*       //      D.ReSize(range); */
/*       cout << " in return datarange" << range2 << " range.value() " << range.value() <<  endl; */

/*       // cout << " range " << range <<  " patch " <<opts.patch.value() << " V.Nrows() " << V.Nrows() << endl; */
/*       //shared(V,indexval,range2) private(i) */


/* #pragma omp parallel */

/*       { */
/*   // cout << " here " << endl; */
 
/* //cout << "omp_get_thread_num" << omp_get_thread_num() << endl; */
/* #pragma omp for nowait */
/*       for( i=0;i<range2;i++) */
/* 	V.push_back(indexval+i); // labels index from 1  */
/*       }     */
/*   } */
/*   else{ */
      
/*       Matrix patchm=read_ascii_matrix(patch.value());  */
      
/*       if(in.get_pvalue(0)){ */
/* 	for (int i=0;i<in.nvertices();i++) */
/* 	  in.set_pvalue(i,0); */
/*       } */

/*       cout << " label load not available for newmesh yet " << endl; */
/*       exit(0); */
/*       // in.load_fs_label(patch.value()); */
/*       // cout << "here " << endl; */
/*       int ind=0; */
/*       for(vector<boost::shared_ptr<Mpoint> >::const_iterator i=in.vbegin();i!=in.vend();i++){ */
/* 	//  cout << " here 1 " << (*i)->get_value() << endl; */
/* 	if(in.get_pvalue((*i)->get_no()) > 0){ */
	  
/* 	  V.push_back((*i)->get_no()+1); /// labels indexing runs from 1 !!! */
/* 	  ind++; */
/* 	  /// note this was in error before inasmuch as if I was using a patch I would count indexes from 0 but the datarange assumed indexing was running from 1 */
/* 	// 	cout << ind << " V(ind) " << V(ind) << endl; */
/*        } */
/*     } */
/*   } */

/*   return V; */

/* } */

#endif

