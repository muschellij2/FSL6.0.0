/*  Copyright (C) 2009 University of Oxford  */

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

#include <iostream>
#include <cmath>
#include "utils/options.h"
#include "miscmaths/miscmaths.h"
#include "newmat.h"
#include "newimage/newimageall.h"

using namespace std;
using namespace NEWMAT;
using namespace MISCMATHS;
using namespace NEWIMAGE;
using namespace Utilities;



string title="dtigen - generate diffusion data using tensor model";
string examples="dtigen -t <input4Dtensor> -o <output4Ddata> -b <bvals> -r <bvecs> -m <brainmask> --s0=<s0file>";

Option<bool> help(string("-h,--help"),false,
		       string("display this message"),
		       false,no_argument);
Option<string> itensor(string("-t,--tensor"),string(""),
		       string("input tensor"),
		       true,requires_argument);
Option<string> s0file(string("--s0"),string(""),
		       string("input S0"),
		       true,requires_argument);
Option<string> odata(string("-o,--output"),string(""),
		       string("output data"),
		       true,requires_argument);
Option<string> bvecsfile(string("-r,--bvecs"),string(""),
		       string("bvecs ASCII text file"),
		       true,requires_argument);
Option<string> bvalsfile(string("-b,--bvals"),string(""),
		       string("bvals ASCII text file"),
		       true,requires_argument);
Option<string> maskfile(string("-m,--mask"),string(""),
		       string("brain mask"),
		       true,requires_argument);
Option<string> kurtfile(string("--kurt"),string(""),
		       string("mean kurtosis map"),
		       false,requires_argument);



Matrix form_Amat(const Matrix& r,const Matrix& b)
{
  Matrix A(r.Ncols(),7);
  Matrix tmpvec(3,1), tmpmat;
  
  for( int i = 1; i <= r.Ncols(); i++){
    tmpvec << r(1,i) << r(2,i) << r(3,i);
    tmpmat = tmpvec*tmpvec.t()*b(1,i);
    A(i,1) = tmpmat(1,1);
    A(i,2) = 2*tmpmat(1,2);
    A(i,3) = 2*tmpmat(1,3);
    A(i,4) = tmpmat(2,2);
    A(i,5) = 2*tmpmat(2,3);
    A(i,6) = tmpmat(3,3);
    A(i,7) = 1;
  }
  return A;
}


Matrix form_Amat_kurt(const Matrix& r,const Matrix& b)
{
  Matrix A(r.Ncols(),8);
  Matrix tmpvec(3,1), tmpmat;
  
  ColumnVector v(r.Ncols());
  for( int i = 1; i <= r.Ncols(); i++){
    v(i) = -b(1,i)*b(1,i)/6;
  }
  Matrix M(r.Ncols(),7);
  M = form_Amat(r,b);

  for( int i = 1; i <= r.Ncols(); i++){
    tmpvec << r(1,i) << r(2,i) << r(3,i);
    tmpmat = tmpvec*tmpvec.t()*b(1,i);
    A(i,1) = tmpmat(1,1);
    A(i,2) = 2*tmpmat(1,2);
    A(i,3) = 2*tmpmat(1,3);
    A(i,4) = tmpmat(2,2);
    A(i,5) = 2*tmpmat(2,3);
    A(i,6) = tmpmat(3,3);
    A(i,7) = 1;
    A(i,8) = v(i); 
  }
  return A;
}



int do_dtigen(){
  volume<float> mask,S0,kurt;
  volume4D<float> data,tensor;

  read_volume(mask,maskfile.value());
  read_volume(S0,s0file.value());
  read_volume4D(tensor,itensor.value());
  
  if(kurtfile.set()){
    read_volume(kurt,kurtfile.value());
  }

  Matrix r = read_ascii_matrix(bvecsfile.value());
  if(r.Nrows()>3) r=r.t();
  for(int i=1;i<=r.Ncols();i++){
    float tmpsum=sqrt(r(1,i)*r(1,i)+r(2,i)*r(2,i)+r(3,i)*r(3,i));
    if(tmpsum!=0){
      r(1,i)=r(1,i)/tmpsum;
      r(2,i)=r(2,i)/tmpsum;
      r(3,i)=r(3,i)/tmpsum;
    }  
  }
  Matrix b = read_ascii_matrix(bvalsfile.value());
  if(b.Nrows()>1) b=b.t();
  if( b.Ncols() != r.Ncols() ){ cerr << "Error: bvecs and bvals don't have the same number of entries" << endl; return(-1);}
  if( r.Nrows() !=3 ){cerr << "Error: bvecs must be either 3xN or Nx3" << endl; return(-1);}
 
  data.reinitialize(mask.xsize(),mask.ysize(),mask.zsize(),b.Ncols());
  copybasicproperties(tensor[0],data);

  Matrix Amat = form_Amat(r,b);
  ColumnVector Dvec(7);
  if(kurtfile.set()){
    Amat = form_Amat_kurt(r,b);
    Dvec.ReSize(8);
  }
  ColumnVector logpred(data.tsize());
  float md;
  cout << "generate data" << endl << endl;;
  for(int z=mask.minz();z<=mask.maxz();z++){
    cout << "processing slice" << z << endl;
    for(int y=mask.miny();y<=mask.maxy();y++)
      for(int x=mask.minx();x<=mask.maxx();x++){
	if(mask(x,y,z)==0)continue;		

	
	  Dvec(1) = tensor(x,y,z,0);
	  Dvec(2) = tensor(x,y,z,1);
	  Dvec(3) = tensor(x,y,z,2);
	  Dvec(4) = tensor(x,y,z,3);
	  Dvec(5) = tensor(x,y,z,4);
	  Dvec(6) = tensor(x,y,z,5);

	  Dvec(7) = -log(S0(x,y,z));
	  
	  if(kurtfile.set()){
	    md = (Dvec(1)+Dvec(4)+Dvec(6))/3.0;
	    Dvec(8) = kurt(x,y,z)*md*md;
	  }
	  logpred = -Amat*Dvec;

	  for(int t=1;t<=data.tsize();t++){
	    data(x,y,z,t-1) = exp(logpred(t));	  
	}

      }
  }
  cout<<"saving results" << endl;
  data.setDisplayMaximumMinimum(1000,0);
  save_volume4D(data,odata.value());


  return 0;
}

int main(int argc,char *argv[]){

  Tracer tr("main");
  OptionParser options(title,examples);

  try{
    options.add(help);
    options.add(itensor);
    options.add(s0file);
    options.add(odata);
    options.add(bvecsfile);
    options.add(bvalsfile);
    options.add(maskfile);
    options.add(kurtfile);

    options.parse_command_line(argc,argv);

    
    if ( (help.value()) || (!options.check_compulsory_arguments(true)) ){
      options.usage();
      exit(EXIT_FAILURE);
    }
  }
  catch(X_OptionError& e) {
    options.usage();
    cerr << endl << e.what() << endl;
    exit(EXIT_FAILURE);
  } 
  catch(std::exception &e) {
    cerr << e.what() << endl;
  } 
  
  return do_dtigen();
  
  
}
