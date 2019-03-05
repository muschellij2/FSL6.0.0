/*  tractOptions.h

    Tim Behrens, FMRIB Image Analysis Group

    Copyright (C) 2004 University of Oxford  */

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

#if !defined(tractOptions_h)
#define tractOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "commonopts.h"

//#include "newmatall.h"
using namespace Utilities;

namespace TRACT {

class tractOptions {
 public:
  static tractOptions& getInstance();
  ~tractOptions() { delete gopt; }
  
  Option<int> verbose;
  Option<bool> help;
  Option<string> basename;
  Option<string> maskfile;
  Option<string> rubbish_file;
  Option<string> skipmask;
  Option<string> outfile;
  Option<string> seedfile; 
  Option<int> nparticles;
  Option<int> nsteps;
  Option<float> steplength;
  Option<float> c_thr;
  Option<bool> modeuler;
  Option<bool> noloopcheck;
  Option<bool> usef;
  Option<float> rseed;
  void parse_command_line(int argc, char** argv);
  void status();
 private:
  tractOptions();  
  const tractOptions& operator=(tractOptions&);
  tractOptions(tractOptions&);

  OptionParser options; 
      
  static tractOptions* gopt;
  
};

 inline tractOptions& tractOptions::getInstance(){
   if(gopt == NULL)
     gopt = new tractOptions();
   
   return *gopt;
 }

 inline tractOptions::tractOptions() :
  verbose(string("-V,--verbose"), 0, 
	  string("verbose level, [0-2]"), 
	  false, requires_argument),
   help(string("-h,--help"), false,
	string("display this message"),
	false, no_argument),
   basename(string("-s,--samples"), string("DTI"),
	       string("basename for samples files"),
	       true, requires_argument),  
   maskfile(string("-m,--mask"), string("mask"),
	    string("Bet binary mask file"),
	    true, requires_argument),
   rubbish_file(string("--rubbish"), string(""),
	    string("Rubbish file (tracking stops here)."),
	    false, requires_argument),
  skipmask(string("--no_integrity"), string(""),
	   string("no explanation needed"),
	   false, requires_argument),
   outfile(string("-o,--out"), string("out"),
	    string("Output file"),
	    true, requires_argument),
   seedfile(string("-x,--seed"), string("Seed"),
	    string("File with seed points in. e.g 68 68 32\n                                                      78 23 52"),
	    true, requires_argument),
   nparticles(string("-P,--nparticles"), 10000,
	 string("Number of particles"),
	 false, requires_argument),
   nsteps(string("-S,--nsteps"), 1000,
	    string("Number of steps per particle"),
	    false, requires_argument),
   steplength(string("-l,steplength"), 0.5, 
	      string("Steplength"), 
	      false, requires_argument),
  c_thr(string("-c,--cthr"), 0.2, 
	string("Curvature threshold"), 
	false, requires_argument),
  modeuler(string("--modeuler"), false, 
	      string("Do modified euler integration instead of simple euler"), 
	      false, no_argument),
  noloopcheck(string("--noloopcheck"), false, 
	 string("Don't perform loopchecking"), 
	 false, no_argument),
  usef(string("-f,--usef"), false, 
	 string("Use anisotropy to constrain tracking"), 
	 false, no_argument),
   rseed(string("--rseed"), 0.324571,
	 string("Random seed"),
	 false, requires_argument), 
   options("tract2","tract2 -s <basename> -m <maskname> -x <seedfile> -o <output>\n tract2 --help\n")
   {
     
    
     try {
       options.add(verbose);
       options.add(help);
       options.add(basename);
       options.add(maskfile);
       options.add(rubbish_file);
       options.add(skipmask);
       options.add(seedfile); 
       options.add(outfile);
       options.add(nparticles);
       options.add(nsteps);
       options.add(steplength);
       options.add(c_thr);
       options.add(modeuler);
       options.add(noloopcheck);
       options.add(usef);
       options.add(rseed);
     }
     catch(X_OptionError& e) {
       options.usage();
       cerr << endl << e.what() << endl;
     } 
     catch(std::exception &e) {
       cerr << e.what() << endl;
     }    
     
   }
}

#endif







