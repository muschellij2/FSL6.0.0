/*  probtrackxOptions.cc

    Tim Behrens, Saad Jbabdi, FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

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

#define WANT_STREAM
#define WANT_MATH

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "probtrackxOptions.h"
#include "utils/options.h"
//#include "newmat.h"
using namespace Utilities;

namespace TRACT {

probtrackxOptions* probtrackxOptions::gopt = NULL;

  void probtrackxOptions::parse_command_line(int argc, char** argv, Log& logger)
  {
    //Do the parsing;
    try{
      for(int a = options.parse_command_line(argc, argv); a < argc; a++) ;
      // setup logger directory
      // if( mode.value()=="help"){
// 	modehelp();
// 	exit(2);
//       }
      if(help.value() || ! options.check_compulsory_arguments())
	{
	  options.usage();
	  exit(2);
	}   

      else{
	//modecheck(); // check all the correct options are set for this mode.	  
	if(forcedir.value())
	  logger.setthenmakeDir(logdir.value(),"probtrackx.log");
	else
	  logger.makeDir(logdir.value(),"probtrackx.log");
	
	cout << "Log directory is: " << logger.getDir() << endl;
	
	// do again so that options are logged
	for(int a = 0; a < argc; a++)
	  logger.str() << argv[a] << " ";
	logger.str() << endl << "---------------------------------------------" << endl << endl;	
	
      }
    
      
    }
    catch(X_OptionError& e){
      cerr<<e.what()<<endl;
      cerr<<"try: probtrackx2 --help"<<endl;
      exit(2);
    }
    
    
    
    
  }
  
  
  void probtrackxOptions::modecheck()
  {
//     bool check=true;
//     string mesg="";
//     if(mode.value()=="simple"){
//       if(outfile.value()==""){
// 	mesg+="You must set an output name in simple mode: -o\n";
// 	check=false;
//       }
//     }
    
    
//     cerr<<mesg;
//     exit(2);
  }
  
  
  
  void probtrackxOptions::status()
  {
    cout<<"basename   "<<basename.value()<<endl;
    cout<<"maskfile   "<<maskfile.value()<<endl;
    cout<<"seeds      "<<seedfile.value()<<endl;
    cout<<"output     "<<outfile.value()<<endl;
    cout<<"verbose    "<<verbose.value()<<endl;
    cout<<"nparticles "<<nparticles.value()<<endl;
    cout<<"nsteps     "<<nsteps.value()<<endl;
    cout<<"usef       "<<usef.value()<<endl;
    cout<<"rseed      "<<rseed.value()<<endl; 
    cout<<"randfib    "<<randfib.value()<<endl;
    cout<<"fibst      "<<fibst.value()<<endl;
}
  
}










