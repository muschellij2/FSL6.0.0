/*  Copyright (C) 2010 University of Oxford  */
/*  Stam Sotiropoulos, Saad Jbabdi    */
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
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "newimage/newimageall.h"

using namespace std;
using namespace NEWIMAGE;
using namespace NEWMAT;


//For a voxel (x,y,z) at Low_res, returns the overlapping HR voxels coordinates
ReturnMatrix get_HRindices(const int x,const int y,const int z, const int xratio, const int yratio, const int zratio){
  Matrix HRindices(xratio*yratio*zratio,3);  //num_HR_voxels x 3 coordinates per voxel

  int count=1;
  for (int Hx=1; Hx<=xratio; Hx++)
    for (int Hy=1; Hy<=yratio; Hy++)
      for (int Hz=1; Hz<=zratio; Hz++){
	HRindices(count,1)=(x+1)*xratio-Hx;
	HRindices(count,2)=(y+1)*yratio-Hy;
	HRindices(count,3)=(z+1)*zratio-Hz;
	count++;
      }

  HRindices.Release();
  return HRindices;
}



//Computes the model predicted signal from the mode of the rubix samples
//An f0 compartment could also be included in the model   
int main ( int argc, char *argv[]){
  if(argc<=2){
    cerr<<" "<<endl;
    cerr<<"usage: rubix_pred rubix_dir modelnum"<<endl<<endl;
    cerr<<"       rubix_dir is the RubiX output directory "<<endl; 
    cerr<<"       modelnum chooses between the spatial model for the LR predictions "<<endl; 
    cerr<<"       1: for sum of attenuations, 2:attenuation of sums "<<endl;
    cerr<<" "<<endl;
    exit(1);
  }
  
  Matrix bvecsLR,bvalsLR, bvecsHR, bvalsHR;   //Read Input 
  volume<float> maskLR, maskHR;
  volume<float> d,d_std,S0,S0LR,f0,R, temp;
  volume4D<float> temp4D;
  vector< volume<float> > f;  
  vector< volume4D<float> > dyads;  
  
  string dir_name=argv[1], temp_name;
  int PVmodelnum=atoi(argv[2]);
  temp_name=dir_name+"/bvalsHR";
  bvalsHR=read_ascii_matrix(temp_name);
  temp_name=dir_name+"/bvecsHR";
  bvecsHR=read_ascii_matrix(temp_name);

  temp_name=dir_name+"/bvalsLR";
  bvalsLR=read_ascii_matrix(temp_name);
  temp_name=dir_name+"/bvecsLR";
  bvecsLR=read_ascii_matrix(temp_name);

  
  if(bvecsLR.Nrows()>3) bvecsLR=bvecsLR.t();   //Make sure the bvecs entries are normalized unit vectors
  if(bvalsLR.Nrows()>1) bvalsLR=bvalsLR.t();
  for(int i=1;i<=bvecsLR.Ncols();i++){
    float tmpsum=sqrt(bvecsLR(1,i)*bvecsLR(1,i)+bvecsLR(2,i)*bvecsLR(2,i)+bvecsLR(3,i)*bvecsLR(3,i));
    if(tmpsum!=0){
      bvecsLR(1,i)=bvecsLR(1,i)/tmpsum;
      bvecsLR(2,i)=bvecsLR(2,i)/tmpsum;
      bvecsLR(3,i)=bvecsLR(3,i)/tmpsum;
    }  
  }
  if(bvecsHR.Nrows()>3) bvecsHR=bvecsHR.t();   
  if(bvalsHR.Nrows()>1) bvalsHR=bvalsHR.t();
  for(int i=1;i<=bvecsHR.Ncols();i++){
    float tmpsum=sqrt(bvecsHR(1,i)*bvecsHR(1,i)+bvecsHR(2,i)*bvecsHR(2,i)+bvecsHR(3,i)*bvecsHR(3,i));
    if(tmpsum!=0){
      bvecsHR(1,i)=bvecsHR(1,i)/tmpsum;
      bvecsHR(2,i)=bvecsHR(2,i)/tmpsum;
      bvecsHR(3,i)=bvecsHR(3,i)/tmpsum;
    }  
  }

  int num_fibres=0;  //Check how many fibres exist
  while (fsl_imageexists(dir_name+"/dyads"+num2str(num_fibres+1)))
    num_fibres++;
 
  temp_name=dir_name+"/mean_dsamples";   //Read rubix local model results 
  if (!fsl_imageexists(temp_name)){
    cout<<"No mean_dsamples file exists!"<<endl;
    exit(1); 
  }
  else read_volume(d,temp_name);
 
  temp_name=dir_name+"/mean_S0samples";   
  if (!fsl_imageexists(temp_name)){
    cerr<<"No mean_S0samples or data file exists!"<<endl;
    exit(1);
  } 
  else read_volume(S0,temp_name);
  
  temp_name=dir_name+"/mean_S0_LRsamples";   
  if (!fsl_imageexists(temp_name)){
    cerr<<"No mean_S0_LRsamples or data file exists!"<<endl;
    exit(1);
  } 
  else read_volume(S0LR,temp_name);

  for (int n=0; n<num_fibres; n++){   //Read dyads
    temp_name=dir_name+"/dyads"+num2str(n+1);
    if (!fsl_imageexists(temp_name)){
      cerr<<"No dyads"<<n+1<<" file exists!"<<endl;
      exit(1); }
    else{
      read_volume4D(temp4D,temp_name);
      dyads.push_back(temp4D);
    }
    temp_name=dir_name+"/mean_f"+num2str(n+1)+"samples";
    if (!fsl_imageexists(temp_name)){
      cerr<<"No mean_f"<<n+1<<"samples file exists!"<<endl;
      exit(1); }
    else{
      read_volume(temp,temp_name);
      f.push_back(temp);
    }
  }

  int modelnum=1; 

  temp_name=dir_name+"/mean_d_stdsamples";    
  if (fsl_imageexists(temp_name)){   //Read d_std if model2
    modelnum=2;
    read_volume(d_std,temp_name);
  }
  temp_name=dir_name+"/mean_dstd_samples";    //some older versions had a different dstd filename
  if (fsl_imageexists(temp_name)){   //Read d_std if model2
    modelnum=2;
    read_volume(d_std,temp_name);
  }

  temp_name=dir_name+"/mean_Rsamples";    
  if (fsl_imageexists(temp_name)){   //Read R if model3
    modelnum=3;
    read_volume(R,temp_name);
  }


  int f0_incl=0;
  temp_name=dir_name+"/mean_f0samples";    
  if (fsl_imageexists(temp_name)){
    f0_incl=1;
    read_volume(f0,temp_name); 
  }

  cout<<"Files for local model"<<modelnum<<" with "<<num_fibres<<" fibres found"<<endl;
  if (f0_incl==1)
    cout<<"Also an f0 noise floor is assumed"<<endl<<endl;

  temp_name=dir_name+"/nodif_brain_mask";    
  if (fsl_imageexists(temp_name))
    read_volume(maskHR,temp_name);
  else{ 
    maskHR.reinitialize(d.xsize(),d.ysize(),d.zsize());
    maskHR=1;
  } 

  temp_name=dir_name+"/nodif_brain_maskLR";    
  if (fsl_imageexists(temp_name))
    read_volume(maskLR,temp_name);
  else{ 
    maskLR.reinitialize(S0LR.xsize(),S0LR.ysize(),S0LR.zsize());
    maskLR=1;
  }


  volume4D<float> outputHR, outputLR;
  outputHR.reinitialize(S0.xsize(),S0.ysize(),S0.zsize(),bvalsHR.Ncols());
  copybasicproperties(S0,outputHR);
  outputHR.setdims(S0.xdim(),S0.ydim(),S0.zdim(),1.0);
  outputHR=0;
  
  outputLR.reinitialize(S0LR.xsize(),S0LR.ysize(),S0LR.zsize(),bvalsLR.Ncols());
  copybasicproperties(S0LR,outputLR);
  outputLR.setdims(S0LR.xdim(),S0LR.ydim(),S0LR.zdim(),1.0);
  outputLR=0;

  float xratio=round(S0LR.xdim()/S0.xdim());
  float yratio=round(S0LR.ydim()/S0.ydim());
  float zratio=round(S0LR.zdim()/S0.zdim());

  for(int k=S0LR.minz();k<=S0LR.maxz();k++){     //Compute predicted signal for each LR and HR voxel
    for(int j=S0LR.miny();j<=S0LR.maxy();j++){
      for(int i=S0LR.minx();i<=S0LR.maxx();i++){ //For each LR voxel
	if (maskLR(i,j,k)!=0){
	  Matrix HRindices;
	  HRindices=get_HRindices(i,j,k,(int)xratio,(int)yratio,(int)zratio);  //find overlapping HR voxels
	  Matrix LRpart; LRpart.ReSize(HRindices.Nrows(),bvalsLR.Ncols()); LRpart=0;
	  ColumnVector S0part; S0part.ReSize(HRindices.Nrows()); S0part=0;

	  for (int v=1; v<=HRindices.Nrows(); v++){ //For each overlapping HR voxel
	    int x=HRindices(v,1); int y=HRindices(v,2); int z=HRindices(v,3);
	    //Predict HR signal
	    for (int l=0; l<bvalsHR.Ncols(); l++){ //for each HR datapoint
	      //Aniso Signal first 
	      float sumf=0; float sig2=0; float dalpha=0;
	      if (modelnum==1 || (modelnum==2 && d_std(x,y,z)<=1e-5)){ //model1 or model2 with small dstd
		sumf=0;
		for (int n=0; n<num_fibres; n++){
		  sumf+=f[n](x,y,z);
		  float angp=dyads[n](x,y,z,0)*bvecsHR(1,l+1)+dyads[n](x,y,z,1)*bvecsHR(2,l+1)+dyads[n](x,y,z,2)*bvecsHR(3,l+1);
		  outputHR(x,y,z,l)+=f[n](x,y,z)*std::exp(-bvalsHR(1,l+1)*d(x,y,z)*angp*angp);
		}
	      }
	      else if (modelnum>=2) { //model2 or model3
		sig2=d_std(x,y,z)*d_std(x,y,z);
		dalpha=d(x,y,z)*d(x,y,z)/sig2;    
		sumf=0;
		for (int n=0; n<num_fibres; n++){
		  sumf+=f[n](x,y,z);
		  float angp=dyads[n](x,y,z,0)*bvecsHR(1,l+1)+dyads[n](x,y,z,1)*bvecsHR(2,l+1)+dyads[n](x,y,z,2)*bvecsHR(3,l+1);
		  if (modelnum==2)
		    outputHR(x,y,z,l)+=f[n](x,y,z)*std::exp(dalpha*log(d(x,y,z)/(d(x,y,z) + bvalsHR(1,l+1)*angp*angp*sig2)));
		  if (modelnum==3)
		    outputHR(x,y,z,l)+=f[n](x,y,z)*std::exp(-bvalsHR(1,l+1)*3*d(x,y,z)/(2*R(x,y,z)+1.0)*((1-R(x,y,z))*angp*angp+R(x,y,z)));
		}
	      }
	      //Iso signal Now
	      if (modelnum==1 || d_std(x,y,z)<=1e-5){
		if (f0_incl==1)
		  outputHR(x,y,z,l)+=f0(x,y,z)+(1-sumf-f0(x,y,z))*std::exp(-bvalsHR(1,l+1)*d(x,y,z));
		else  
		  outputHR(x,y,z,l)+=(1-sumf)*std::exp(-bvalsHR(1,l+1)*d(x,y,z));
	      }
	      else{
		if (f0_incl==1)
		  outputHR(x,y,z,l)+=f0(x,y,z)+(1-sumf-f0(x,y,z))*std::exp(dalpha*log(d(x,y,z)/(d(x,y,z)+bvalsHR(1,l+1)*sig2)));
		else  
		  outputHR(x,y,z,l)+=(1-sumf)*std::exp(dalpha*log(d(x,y,z)/(d(x,y,z)+bvalsHR(1,l+1)*sig2)));
	      }
	      outputHR(x,y,z,l)*=S0(x,y,z); 
	    }

	    //Predict contribution to LR signal
	    for (int l=0; l<bvalsLR.Ncols(); l++){ //for each LR datapoint
	      float sumf=0; float sig2=0; float dalpha=0;
	      //Aniso Signal first
	      if (modelnum==1 || (modelnum==2 && d_std(x,y,z)<=1e-5)){ //model1 or model2 with small dstd
		sumf=0;
		for (int n=0; n<num_fibres; n++){
		  sumf+=f[n](x,y,z);
		  float angp=dyads[n](x,y,z,0)*bvecsLR(1,l+1)+dyads[n](x,y,z,1)*bvecsLR(2,l+1)+dyads[n](x,y,z,2)*bvecsLR(3,l+1);
		  LRpart(v,l+1)+=f[n](x,y,z)*std::exp(-bvalsLR(1,l+1)*d(x,y,z)*angp*angp);
		}
	      }
	      else if (modelnum>=2) { //model2 or model3
		sig2=d_std(x,y,z)*d_std(x,y,z);
		dalpha=d(x,y,z)*d(x,y,z)/sig2; 
		sumf=0;
		for (int n=0; n<num_fibres; n++){
		  sumf+=f[n](x,y,z);
		  float angp=dyads[n](x,y,z,0)*bvecsLR(1,l+1)+dyads[n](x,y,z,1)*bvecsLR(2,l+1)+dyads[n](x,y,z,2)*bvecsLR(3,l+1);
		  if (modelnum==2)
		     LRpart(v,l+1)+=f[n](x,y,z)*std::exp(dalpha*log(d(x,y,z)/(d(x,y,z) + bvalsLR(1,l+1)*angp*angp*sig2)));
		  if (modelnum==3)
		    LRpart(v,l+1)+=f[n](x,y,z)*std::exp(-bvalsLR(1,l+1)*3*d(x,y,z)/(2*R(x,y,z)+1.0)*((1-R(x,y,z))*angp*angp+R(x,y,z)));
		}
	      }
	      //Iso signal Now
	      if (modelnum==1 || d_std(x,y,z)<=1e-5){
		if (f0_incl==1)
		  LRpart(v,l+1)+=f0(x,y,z)+(1-sumf-f0(x,y,z))*std::exp(-bvalsLR(1,l+1)*d(x,y,z));
		else  
		  LRpart(v,l+1)+=(1-sumf)*std::exp(-bvalsLR(1,l+1)*d(x,y,z));
	      }
	      else{
		if (f0_incl==1)
		  LRpart(v,l+1)+=f0(x,y,z)+(1-sumf-f0(x,y,z))*std::exp(dalpha*log(d(x,y,z)/(d(x,y,z)+bvalsLR(1,l+1)*sig2)));
		else  
		  LRpart(v,l+1)+=(1-sumf)*std::exp(dalpha*log(d(x,y,z)/(d(x,y,z)+bvalsLR(1,l+1)*sig2)));
	      }
	      LRpart(v,l+1)*=S0(x,y,z);
	      S0part(v)=S0(x,y,z);
	    } 
	  }
	  //Predict LR signal out of the HR voxels
	  if (PVmodelnum==1){ //old partial volume model (sum of attenuations)
	    for (int l=0; l<bvalsLR.Ncols(); l++){ //for each LR datapoint
	      float sumf=0;
	      for (int v=1; v<=HRindices.Nrows(); v++) //For each HR voxel
		sumf+=1.0/HRindices.Nrows()*LRpart(v,l+1)/S0part(v); 
	      outputLR(i,j,k,l)=S0LR(i,j,k)*sumf;
	    }
	  }
	  else{    //attenuation of sums
	    for (int l=0; l<bvalsLR.Ncols(); l++){ //for each LR datapoint
	      float sumS0=0; float sumf=0;
	      for (int v=1; v<=HRindices.Nrows(); v++){ //For each HR voxel
		sumf+=LRpart(v,l+1);
		sumS0+=S0part(v);
	      }
	      outputLR(i,j,k,l)=S0LR(i,j,k)*sumf/sumS0;
	    }
	  } 
	}
      }
    }
    cout<<k+1<<" LR slices processed"<<endl;
  }
  cout<<"saving results"<<endl;
  temp_name=dir_name+"/dataHR_pred";    
  save_volume4D(outputHR,temp_name); 
  temp_name=dir_name+"/dataLR_pred";    
  save_volume4D(outputLR,temp_name); 
  //  temp_name=dir_name+"/dataHRLR_pred";    
  //save_volume4D(outputHRLR,temp_name); 
  return 0;
}









