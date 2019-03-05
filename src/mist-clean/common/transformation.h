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

#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include <string>
#include <vector>
#include <boost/log/trivial.hpp>
#include "newmat.h"
// removed because warpfns.h needs to define EXPOSE_TREACHEROUS before including this
// #include "newimage.h"
#include "warpfns/warpfns.h"

// NB: All classes in here have default operator=, copy ctor, etc!


class IdentityTransformation;
class AffineTransformation;
class NonlinearTransformation;

class Transformation
{
public:
    class TransformationException : public std::runtime_error
    {
    public:
        TransformationException(const std::string &what) : runtime_error(what) { };
    };

    Transformation();
    virtual ~Transformation() { }

    virtual void ReleaseMemory() const { }

    virtual NEWMAT::ColumnVector InverseTransformPoint(const NEWMAT::ColumnVector &p) const = 0;

    // Can't think of a better way :(
    template<class T, class U>
    NEWIMAGE::volume<T> TransformVolume(const NEWIMAGE::volume<T> &invol, const NEWIMAGE::volume<U> &refvol) const
    {
        return apply_derived_transformation(this, invol, refvol);
    }
};


class IdentityTransformation : public Transformation
{
public:
    virtual NEWMAT::ColumnVector InverseTransformPoint(const NEWMAT::ColumnVector &p) const;

    template<class T, class U>
    NEWIMAGE::volume<T> TransformVolumeImplementation(const NEWIMAGE::volume<T> &invol, const NEWIMAGE::volume<U> &refvol) const
    {
        BOOST_LOG_TRIVIAL(warning) << "IdentityTransformation::TransformVolumeImplementation() is untested!";

        Matrix identity(4, 4);
        identity << 1.0 << 0.0 << 0.0 << 0.0
                 << 0.0 << 1.0 << 0.0 << 0.0
                 << 0.0 << 0.0 << 1.0 << 0.0
                 << 0.0 << 0.0 << 0.0 << 1.0;

        NEWIMAGE::volume<T> outvol(refvol.xsize(), refvol.ysize(), refvol.zsize());
        NEWIMAGE::copybasicproperties(refvol, outvol);

        // No warp
        NEWIMAGE::volume4D<float> warp;
        vector<int> defdir;

        // No derivatives
        vector<int> derivdir;
        NEWIMAGE::volume4D<T> deriv;

        NEWIMAGE::raw_general_transform(invol, identity, warp, defdir, derivdir, outvol, deriv);

        return outvol;
    }
};


class AffineTransformation : public Transformation
{
public:
    AffineTransformation(const std::string &filename);
    AffineTransformation(const NEWMAT::Matrix &matrix);

    AffineTransformation Concatenate(const AffineTransformation &second) const;
    NonlinearTransformation Concatenate(const NonlinearTransformation &second) const;

    virtual NEWMAT::ColumnVector InverseTransformPoint(const NEWMAT::ColumnVector &p) const;

    template<class T, class U>
    NEWIMAGE::volume<T> TransformVolumeImplementation(const NEWIMAGE::volume<T> &invol, const NEWIMAGE::volume<U> &refvol) const
    {
        NEWIMAGE::volume<T> outvol(refvol.xsize(), refvol.ysize(), refvol.zsize());
        NEWIMAGE::copybasicproperties(refvol, outvol);

        // No warp
        NEWIMAGE::volume4D<float> warp;
        vector<int> defdir;

        // No derivatives
        vector<int> derivdir;
        NEWIMAGE::volume4D<T> deriv;

        NEWIMAGE::raw_general_transform(invol, m_inverse.i(), warp, defdir, derivdir, outvol, deriv);

        return outvol;
    }

private:
    friend class NonlinearTransformation;

    Matrix m_inverse;
};


class NonlinearTransformation : public Transformation
{
public:
    NonlinearTransformation(const std::string &filename, bool includeaffine);
    NonlinearTransformation(const NEWIMAGE::volume4D<float> &warpfield);

    NonlinearTransformation Concatenate(const NonlinearTransformation &second) const;

    template<class T>
    NonlinearTransformation Concatenate(const AffineTransformation &second, const NEWIMAGE::volume<T> &reference) const
    {
        BOOST_LOG_TRIVIAL(warning) << "NonlinearTransformation::Concatenate(AffineTransformation&) is untested!";

        NEWIMAGE::volume4D<float> secondfield;
        affine2warp(second.m_inverse.i(), secondfield, reference);

        NEWIMAGE::volume4D<float> combinedfield;
        concat_warps(AccessField(), secondfield, combinedfield);

        return NonlinearTransformation(combinedfield);
    }

    virtual NEWMAT::ColumnVector InverseTransformPoint(const NEWMAT::ColumnVector &p) const;

    template<class T, class U>
    NEWIMAGE::volume<T> TransformVolumeImplementation(const NEWIMAGE::volume<T> &invol, const NEWIMAGE::volume<U> &refvol) const
    {
        NEWIMAGE::volume<T> outvol(refvol.xsize(), refvol.ysize(), refvol.zsize());
        NEWIMAGE::copybasicproperties(refvol, outvol);

        NEWIMAGE::volume4D<float> relwarp(AccessField());
        NEWIMAGE::convertwarp_abs2rel(relwarp);
        NEWIMAGE::apply_warp(invol, outvol, relwarp);

        return outvol;
    }

    virtual void ReleaseMemory() const;

private:
    friend class AffineTransformation;

    mutable bool m_ready;
    std::string m_backingFile;
    bool m_includeAffineFromBackingFile;

    mutable NEWIMAGE::volume4D<float> m_field;

    NEWIMAGE::volume4D<float> &AccessField() const;

    void Load() const;
};


template<class T, class U>
NEWIMAGE::volume<T> apply_derived_transformation(const Transformation *xfm, const NEWIMAGE::volume<T> &invol,
                                                 const NEWIMAGE::volume<U> &refvol)
{
    if (auto ix = dynamic_cast<const IdentityTransformation *>(xfm))
        return ix->TransformVolumeImplementation(invol, refvol);
    else if (auto ax = dynamic_cast<const AffineTransformation *>(xfm))
        return ax->TransformVolumeImplementation(invol, refvol);
    else if (auto nx = dynamic_cast<const NonlinearTransformation *>(xfm))
        return nx->TransformVolumeImplementation(invol, refvol);

    throw std::logic_error("Cannot downcast Transformation object");
}


#endif // TRANSFORMATION_H
