/*  readoptions.h

    Michael Chappell - FMRIB Image Analysis Group

    Copyright (C) 2009 University of Oxford  */

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

#if !defined(ReadOptions_h)
#define ReadOptions_h

#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include "utils/options.h"
#include "utils/log.h"

using namespace Utilities;

namespace OXASL {

class ReadOptions {
 public:
  static ReadOptions& getInstance();
  ~ReadOptions() { delete ropt; }
  
  Option<bool> help;

  Option<string> datafile;
  Option<string> maskfile;
  Option<string> outname;
  Option<string> aif;
  Option<float> dt;

  Option<string> metric;
  Option<float> mthresh;

  // timing correction
  Option<bool> tcorrect;
  Option<string> bata;
  Option<string> batt;
  Option<bool> batout;
  Option<float> T1;
  Option<float> fa;

  // std deviation estimation
  Option<bool> std;
  Option<int> nwb;

  void parse_command_line(int argc, char** argv);
  
 private:
  ReadOptions();  
  const ReadOptions& operator=(ReadOptions&);
  ReadOptions(ReadOptions&);

  OptionParser options; 
      
  static ReadOptions* ropt;
  
};

 inline ReadOptions& ReadOptions::getInstance(){
   if(ropt == NULL)
     ropt = new ReadOptions();
   
   return *ropt;
 }

 inline ReadOptions::ReadOptions() :

help(string("-h,--help"), false,
		    string("display this message"),
		    false, no_argument),
   //input files
   datafile(string("--data,--datafile"), string("data"),
			  string("ASL data file"),
		     true, requires_argument),  
   maskfile(string("--mask"), string("maskfile"),
	    string("mask"),
	    true, requires_argument),
   outname(string("--out"), string("asl_mfree"),
	   string("Output directory name"),
	   true, requires_argument),
   aif(string("--aif"),string(""),
       string("Arterial input functions for deconvolution (4D volume, one aif for each voxel within mask)"),
       true,requires_argument),
   dt(string("--dt"),1.0,
      string("Temporal spacing in data (s)\n"),
      true,requires_argument),

   metric(string("--metric"),string(""),
	  string("Metric image"),
	  false,requires_argument),
   mthresh(string("--mthresh"),0.0,
	   string("Metric threshold\n"),
	   false,requires_argument),

   tcorrect(string("--tcorrect"),false,
	    string("Apply correction for timing difference between AIF and tissue curve"),
	    false,no_argument),
   bata(string("--bata"),string(""),
	string("arterial BAT image"),
	false,requires_argument),
   batt(string("--batt"),string(""),
	string("tissue BAT image"),
	false,requires_argument),
   batout(string("--bat"),false,
	  string("Estimate tissue BAT from data (and save this image)"),
	  false,no_argument),
   T1(string("--t1"),1.6,
      string("T1 (of blood) value"),
      false,requires_argument),
   fa(string("--fa"),0.0,
      string("Flip anlge (is using LL readout)"),
      false,requires_argument),
   
  
   std(string("--std"),false,
       string("Calculate standard deviations on perfusion values using wild bootstrapping"),
       false,no_argument),
   nwb(string("--nwb"),1000,
       string("Number of permuatations for wild bootstrapping"),
       false,requires_argument),
   

   options("asl_mfree","asl_mfree --verbose\n")
   {
     try {
       options.add(help);

       options.add(datafile);
       options.add(maskfile);
       options.add(outname);
       options.add(aif);
       options.add(dt);

       options.add(metric);
       options.add(mthresh);

       options.add(tcorrect);
       options.add(bata);
       options.add(batt);
       options.add(batout);
       options.add(T1);
       options.add(fa);
       
       options.add(std);
       options.add(nwb);

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



