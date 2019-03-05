/*  split_parts_gpu.cc

    Tim Behrens, Saad Jbabdi, Stam Sotiropoulos, Moises Hernandez  - FMRIB Image Analysis Group

    Copyright (C) 2005 University of Oxford  */

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

#include "newimage/newimageall.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

void save_part(Matrix data, string name, int idpart){

	int nvox = data.Ncols();
	int ndir = data.Nrows();

	string file_name;
	file_name = name+"_"+num2str(idpart);

	ofstream out;
	out.open(file_name.data(), ios::out | ios::binary);
	out.write((char*)&data(1,1),nvox*ndir*sizeof(Real));
	out.close();
}


// parameters:
// 1. data.nii.gz
// 2. mask.nii.gz
// 3. bvals
// 4. bvecs
// 5. grad_dev.nii.gz
// 6. gflag
// 7. nparts
// 8. output directory


int main(int argc, char *argv[]){

	///////////////////////////////////////////
	///////////// Check Arguments /////////////
	///////////////////////////////////////////
	if (argc!=9){
		cerr << "\nsplit_parts_gpu. Wrong number of parameters\nUsage:\nsplit_parts_gpu  Datafile  Maskfile   bvals_file   bvecs_file  Grad_file(can be null)  Use_grad_file(0 or 1)  TotalNumParts  OutputDirectory\n" << endl;
    		exit (EXIT_FAILURE);
	}

	std::string data_str = argv[1];	
	std::string mask_str = argv[2];

	Matrix bvals,bvecs;
    	bvals=read_ascii_matrix(argv[3]);
    	bvecs=read_ascii_matrix(argv[4]);
    	if(bvecs.Nrows()>3) bvecs=bvecs.t();
    	if(bvals.Nrows()>1) bvals=bvals.t();
	if(bvecs.Nrows()!=3){
		cerr << "split_parts_gpu error: bvecs is not 3xN or Nx3 format\n" << endl;
    		exit (EXIT_FAILURE);
    	}
    	for(int i=1;i<=bvecs.Ncols();i++){
      		float tmpsum=sqrt(bvecs(1,i)*bvecs(1,i)+bvecs(2,i)*bvecs(2,i)+bvecs(3,i)*bvecs(3,i));
      		if(tmpsum!=0){
			bvecs(1,i)=bvecs(1,i)/tmpsum;
			bvecs(2,i)=bvecs(2,i)/tmpsum;
			bvecs(3,i)=bvecs(3,i)/tmpsum;
      		}  
    	}
	
	std::string grad_str = argv[5];

	istringstream ss_gflag(argv[6]);
	int gflag;
	if (!(ss_gflag >> gflag)){
		cerr << "\nsplit_parts_gpu. Wrong flag Use_grad_file: " << argv[6] <<  "\nUsage:\nsplit_parts_gpu  Datafile  Maskfile   bvals_file   bvecs_file  Grad_file(can be null)  Use_grad_file(0 or 1)  TotalNumParts  OutputDirectory\n" << endl;
    		exit (EXIT_FAILURE);
	}

	istringstream ss_nparts(argv[7]);
	int nparts;
	if (!(ss_nparts >> nparts)){
		cerr << "\nsplit_parts_gpu. Wrong number of parts: " << argv[7] <<  "\nUsage:\nsplit_parts_gpu Datafile  Maskfile   bvals_file   bvecs_file  Grad_file(can be null)  Use_grad_file(0 or 1)  TotalNumParts  OutputDirectory\n" << endl;
    		exit (EXIT_FAILURE);
	}

	std::string out_dir = argv[8];
	struct stat sb;
	if (stat(out_dir.data(), &sb) != 0 || !S_ISDIR(sb.st_mode)){
		cerr << "\nsplit_parts_gpu. Wrong output directory: " << argv[8] <<  "\nUsage:\nsplit_parts_gpu  Datafile  Maskfile  Grad_file(can be null)  Use_grad_file(0 or 1)  TotalNumParts  OutputDirectory\n" << endl;
    		exit (EXIT_FAILURE);
	}
	///////////////////////////////////////////

	NEWIMAGE::volume4D<float> data;
	NEWIMAGE::volume<float> mask;
    	read_volume4D(data,data_str);
    	read_volume(mask,mask_str);
	Matrix datam;
    	datam=data.matrix(mask); 

	int ndirections=bvals.Ncols();
    	if(bvecs.Ncols()!=ndirections || datam.Nrows()!=ndirections){
		cerr << "split_parts_gpu error: The number of elements in bvals, number of vectors in bvecs and number of vols in data must be the same\n" << endl;
    		exit (EXIT_FAILURE);
    	}
	int nvoxels=datam.Ncols();
	if(nvoxels<=0 || ndirections<=0){
		cerr << "The number of voxels and gradient directions must be greater than 0" << endl;
    		exit (EXIT_FAILURE);
	}
	
	NEWIMAGE::volume4D<float> grad; 
	Matrix gradm;
	int dirs_grad=0;
	if(gflag){
		read_volume4D(grad,grad_str);
      		gradm=grad.matrix(mask);
		dirs_grad = gradm.Nrows();
	}
	
	int size_part=nvoxels/nparts;

	Matrix data_part;
	Matrix grad_part;
	string out_data;
	string out_grad;
	out_data.append(out_dir);
	out_grad.append(out_dir);
	out_data.append("/data");
	out_grad.append("/grad_dev");
	for(int i=0;i<(nparts-1);i++){
		data_part = datam.SubMatrix(1,ndirections,i*size_part+1,(i+1)*size_part);
		save_part(data_part,out_data,i);
		if(gflag){
			grad_part = gradm.SubMatrix(1,dirs_grad,i*size_part+1,(i+1)*size_part);
			save_part(grad_part,out_grad,i);
		}
	}
	// last part
	data_part = datam.SubMatrix(1,ndirections,(nparts-1)*size_part+1,nvoxels);
	save_part(data_part,out_data,(nparts-1));
	if(gflag){
		grad_part = gradm.SubMatrix(1,dirs_grad,(nparts-1)*size_part+1,nvoxels);
		save_part(grad_part,out_grad,(nparts-1));
	}
}
