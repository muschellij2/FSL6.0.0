/*  Multimodal Image Segmentation Tool (MIST)  */
/*  Eelke Visser  */

/*  Copyright (c) 2016 University of Oxford  */

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

#include "transformation.h"
#include "miscmaths/miscmaths.h"
#include "warpfns/fnirt_file_reader.h"

using namespace NEWIMAGE;
using namespace NEWMAT;

Transformation::Transformation()
{
}

ColumnVector IdentityTransformation::InverseTransformPoint(const ColumnVector &p) const
{
    return p;
}

AffineTransformation::AffineTransformation(const string &filename)
{
    BOOST_LOG_TRIVIAL(debug) << "AffineTransformation: Loading matrix " << filename;

    Matrix matrix = MISCMATHS::read_ascii_matrix(4, 4, filename);

    if (matrix.IsZero())
        throw TransformationException("Cannot load affine matrix");

    m_inverse = matrix.i();
}

AffineTransformation::AffineTransformation(const Matrix &matrix)
    : m_inverse(matrix.i())
{

}

AffineTransformation AffineTransformation::Concatenate(const AffineTransformation &second) const
{
    BOOST_LOG_TRIVIAL(warning) << "AffineTransformation::Concatenate(AffineTransformation&) is untested!";

    return AffineTransformation(m_inverse * second.m_inverse);
}

NonlinearTransformation AffineTransformation::Concatenate(const NonlinearTransformation &second) const
{
    volume4D<float> firstfield;
    affine2warp(m_inverse.i(), firstfield, second.AccessField()[0]);

    volume4D<float> combinedfield;
    concat_warps(firstfield, second.AccessField(), combinedfield);

    return NonlinearTransformation(combinedfield);
}

ColumnVector AffineTransformation::InverseTransformPoint(const ColumnVector &p) const
{
    ColumnVector ap(4);
    ap.Rows(1, 3) = p;
    ap(4) = 1.0;

    return (m_inverse * ap).Rows(1, 3);
}

NonlinearTransformation::NonlinearTransformation(const string &filename, bool includeaffine)
    : m_ready(false),
      m_backingFile(filename),
      m_includeAffineFromBackingFile(includeaffine)
{

}

NonlinearTransformation::NonlinearTransformation(const NEWIMAGE::volume4D<float> &warpfield)
    : m_ready(true),
      m_includeAffineFromBackingFile(false),
      m_field(warpfield)
{

}

void NonlinearTransformation::Load() const
{
    BOOST_LOG_TRIVIAL(debug) << "NonlinearTransformation: Loading warp field " << m_backingFile
                             << " (includeaffine = " << m_includeAffineFromBackingFile << ")";

    FnirtFileReader ffr(m_backingFile);
    m_field = ffr.FieldAsNewimageVolume4D(m_includeAffineFromBackingFile);
    convertwarp_rel2abs(m_field);

    m_ready = true;
}

volume4D<float> &NonlinearTransformation::AccessField() const
{
    if (!m_ready)
        Load();

    return m_field;
}

void NonlinearTransformation::ReleaseMemory() const
{
    if (!m_backingFile.empty())
    {
        BOOST_LOG_TRIVIAL(debug) << "NonlinearTransformation: Releasing warp field " << m_backingFile;

        m_ready = false;
        m_field = volume4D<float>();
    }
}

NonlinearTransformation NonlinearTransformation::Concatenate(const NonlinearTransformation &second) const
{
    BOOST_LOG_TRIVIAL(warning) << "NonlinearTransformation::Concatenate(NonlinearTransformation&) is untested!";

    volume4D<float> combinedfield;
    concat_warps(AccessField(), second.AccessField(), combinedfield);

    return NonlinearTransformation(combinedfield);
}

ColumnVector NonlinearTransformation::InverseTransformPoint(const ColumnVector &p) const
{
    if (!m_ready)
        Load();

    ColumnVector d(3);
    // Trilinear by default
    for (int i = 0; i < 3; i++)
        d(i + 1) = m_field[i].interpolate(p(1) / m_field.xdim(), p(2) / m_field.ydim(), p(3) / m_field.zdim());

    return d;
}
