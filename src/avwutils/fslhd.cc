//     fslhd.cc - show image header
//     Steve Smith, Mark Jenkinson and Matthew Webster, FMRIB Image Analysis Group
//     Copyright (C) 2000-2012 University of Oxford  
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
#include <iomanip>
#include <iostream>
using namespace NEWIMAGE;

void reportXMLHeader(const NiftiHeader& header)
{
  cout << "<nifti_image" << endl;
  cout << "  image_offset = '" << header.vox_offset  << "'" << endl;
  cout << "  ndim = '" << header.dim[0] << "'" << endl;
  cout << "  nx = '" << header.dim[1] << "'" << endl;
  cout << "  ny = '" << header.dim[2] << "'" << endl;
  cout << "  nz = '" << header.dim[3] << "'" << endl;
  cout << "  nt = '" << header.dim[4] << "'" << endl;
  if ( header.dim[0] > 4 ) {
    cout << "  nu = '" << header.dim[5] << "'" << endl;
    cout << "  nv = '" << header.dim[6] << "'" << endl;
    cout << "  nw = '" << header.dim[7] << "'" << endl;
  }
  cout << "  dx = '" << header.pixdim[1] << "'" << endl;
  cout << "  dy = '" << header.pixdim[2] << "'" << endl;
  cout << "  dz = '" << header.pixdim[3] << "'" << endl;
  cout << "  dt = '" << header.pixdim[4] << "'" << endl;
  if ( header.dim[0] > 4 ) {
    cout << "  du = '" << header.pixdim[5] << "'" << endl;
    cout << "  dv = '" << header.pixdim[6] << "'" << endl;
    cout << "  dw = '" << header.pixdim[7] << "'" << endl;
  }
  cout << "  datatype = '" << header.datatype << "'" << endl;
  cout << "  nvox = '" << header.nElements() << "'" << endl;
  cout << "  nbyper = '" << header.bitsPerVoxel/8 << "'" << endl;
  cout << "  scl_slope = '" << header.sclSlope << "'" << endl;
  cout << "  scl_inter = '" << header.sclInter << "'" << endl;
  cout << "  intent_code = '" << header.intentCode << "'" << endl;
  cout << "  intent_p1 = '" << header.intent_p1 << "'" << endl;
  cout << "  intent_p2 = '" << header.intent_p2 << "'" << endl;
  cout << "  intent_p3 = '" << header.intent_p3 << "'" << endl;
  cout << "  intent_name = '" << header.intentName << "'" << endl;
  cout << "  toffset = '" << header.toffset << "'" << endl;
  cout << "  xyz_units = '" << XYZT_TO_SPACE(header.units) << "'" << endl;
  cout << "  time_units = '" << XYZT_TO_TIME(header.units) << "'" << endl;
  cout << "  freq_dim = '" << (int)header.freqDim() << "'" << endl;
  cout << "  phase_dim = '" << (int)header.phaseDim() << "'" << endl;
  cout << "  slice_dim = '" << (int)header.sliceDim() << "'" << endl;
  cout << "  descrip = '" << header.description << "'" << endl;
  cout << "  aux_file = '" << header.auxillaryFile << "'" << endl;
  cout << "  qform_code = '" << header.qformCode << "'" << endl;
  cout << "  qfac = '" << header.leftHanded() << "'" << endl;
  cout << "  quatern_b = '" << header.qB << "'" << endl;
  cout << "  quatern_c = '" << header.qC << "'" << endl;
  cout << "  quatern_d = '" << header.qD << "'" << endl;
  cout << "  qoffset_x = '" << header.qX << "'" << endl;
  cout << "  qoffset_y = '" << header.qY << "'" << endl;
  cout << "  qoffset_z = '" << header.qZ << "'" << endl;
  cout << "  sform_code = '" << header.sformCode << "'" << endl;
  cout << "  sto_xyz_matrix = '"; 
  mat44 output(header.getSForm());
  for (int i=0;i<4;i++)
    for (int j=0;j<4;j++) {
     cout << output.m[i][j]; 
     if ( j==3 && i==3 )
       cout << "'" << endl;
     else
       cout << " ";
    }
  cout << "  slice_code = '" << header.sliceCode << "'" << endl;
  cout << "  slice_start = '" << header.sliceStart << "'" << endl;
  cout << "  scl_end = '" << header.sliceEnd << "'" << endl;
  cout << "  scl_duration = '" << header.sliceDuration << "'" << endl;
  cout << "/>" << endl;
}

int print_usage(const string& progname) 
{
  cout << endl;
  cout << "Usage: fslhd [-x] <input>" << endl;
  cout << "       -x : instead print an XML-style NIFTI header" << endl;
  return 1;
}

int main(int argc,char *argv[])
{
  if (argc < 2) 
    return print_usage(string(argv[0]));
  NiftiIO reader;
  try {
    string filename(return_validimagefilename(argv[argc-1]));
    NiftiHeader header(reader.loadHeader(filename));
    if ( argc==3 && strcmp(argv[1],"-x")==0 )
      reportXMLHeader(header);
    else {
      cout << "filename\t" << filename << endl;
      header.report();
    }
  }  catch ( exception& e ) { cerr << e.what() << endl; return 1;}
  return 0;
}


