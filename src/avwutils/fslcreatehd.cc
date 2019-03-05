//     fslcreatehd.cc - Copy certain parts of an AVW header
//     Mark Jenkinson, Steve Smith and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2001-2018 University of Oxford
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
#include <fstream>
#include <iostream>

using namespace NEWIMAGE;

void print_usage(const string& progname)
{
  cout << endl;
  cout << "Usage: fslcreatehd <xsize> <ysize> <zsize> <tsize> <xvoxsize> <yvoxsize> <zvoxsize> <tr> <xorigin> <yorigin> <zorigin> <datatype> <headername>" << endl;
  cout << "       fslcreatehd <nifti_xml_file> <headername>" << endl;
  cout << "  Datatype values: " << DT_UNSIGNED_CHAR << "=char, " << DT_SIGNED_SHORT  << "=short, " <<  DT_SIGNED_INT <<"=int, " << DT_FLOAT << "=float, " <<  DT_DOUBLE<< "=double" << endl;
  cout << "  In the first form, a radiological image will be created, with the origin input as a voxel co-ordinate." << endl;
  cout << "  If the output file already exists, its data ( but not geometric information ) will be copied if it has" << endl;
  cout << "  a matching number of elements." << endl;
  cout << "  In the second form, an XML-ish form of nifti header is read (as output by fslhd -x)" << endl;
  cout << "  Note that stdin is used if '-' is used in place of a filename" << endl;

}


int fslcreatehd_main(int argc, char *argv[])
{
  NiftiIO reader;
  vector <NiftiExtension> extensions;
  NiftiHeader header,originalHeader;
  char *buffer(NULL);
  string filename;
  int fileread(1);
  bool existingImage(false);

  if (argc==3) /* use the XML form of header specification */
    filename = string(argv[2]);
  else
    filename = string(argv[13]);

  /* check if file already exists and if so, read the image contents */
  if (FslFileExists(filename)) {
    existingImage = true;
    filename=return_validimagefilename(filename);
    header=reader.loadImage(filename,buffer,extensions);
    originalHeader=header;
  }


  if (argc>3) {
    /* set uninteresting defaults */
    if (!existingImage)
      header.datatype=atoi(argv[12]);
    header.dim[0]=4;
    for (int i=1;i<=4;i++)
      header.dim[i]=atoi(argv[i]);
    for (int i=1;i<=4;i++)
      header.pixdim[i]=atof(argv[i+4]);
      header.sformCode=2;
      header.sX[0]=-header.pixdim[1];
      header.sY[1]=header.pixdim[2];
      header.sZ[2]=header.pixdim[3];
      header.sX[3]=header.pixdim[1]*(header.dim[1]-1.0-atoi(argv[9]));
      header.sY[3]=-header.pixdim[2]*atoi(argv[10]);
      header.sZ[3]=-header.pixdim[3]*atoi(argv[11]);
      header.qformCode=2;
      header.setQForm(header.getSForm());

  } else {
      /* read XML form */
      char *newstr, *oldfname, *oldiname;
      ifstream inputfile;
      vector<string> settings;

      if (strcmp(argv[1],"-")==0) {fileread=0;}
      newstr = (char *)calloc(10000,1);
      oldfname = (char *)calloc(10000,1);
      oldiname = (char *)calloc(10000,1);
      if (fileread)
      {
	      inputfile.open (argv[1], ifstream::in | ifstream::binary);
        if (!inputfile.is_open())
        {
	      cerr << "Cannot open file " << argv[1] << endl;
	      return EXIT_FAILURE;
	}
      }

      do
      {
	  if (fileread)
          {
	    inputfile.getline(newstr,9999);  // maybe use > for delimiting character remove while increase size
	  }
          else
          {
	      if (fgets(newstr,9999,stdin)==NULL) break;
	  }
	  settings.push_back(string(newstr));
      }  while (settings[settings.size()-1]!="/>");

      header=NiftiHeader(settings);

      if (fileread)
	inputfile.close();
  }

  /* reset datatype in case it has been overwritten */
  if (existingImage) {
    header.datatype=originalHeader.datatype;
    header.bitsPerVoxel=header.bpvOfDatatype();
  }

 /* if previously read buffer is wrong size then make a zero image here */
  if ( !existingImage || header.nElements() != originalHeader.nElements() ) {
    if(buffer!=NULL)
      delete buffer;
    buffer = new char[header.nElements()*(header.bpvOfDatatype()/8)];
    fill(buffer,buffer+header.nElements(),0);
  }

    int filetype=FslGetEnvOutputType();
    if(!existingImage) {
      filename=make_basename(filename)+outputExtension(filetype);
    }

    header.bitsPerVoxel=header.bpvOfDatatype();
    header.setNiftiVersion(FslNiftiVersionFileType(filetype),FslIsSingleFileType(filetype));
    reader.saveImage(filename,buffer,extensions,header, FslIsCompressedFileType(filetype));

  return 0;
  }


int main(int argc,char *argv[])
{
  if (argc != 14 && argc != 3)
  {
    print_usage(string(argv[0]));
    return 1;
  }
  return fslcreatehd_main(argc,argv);
}
