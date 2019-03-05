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

#ifndef SERIALISATION_H
#define SERIALISATION_H

#include <boost/log/trivial.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/collections_save_imp.hpp>
#include <boost/serialization/collections_load_imp.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/split_free.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include "newmat.h"
#include <unordered_map>

// NOTE: According to string.hpp, strings are treated as primitive types
//       (important when creating temporary strings in load/save etc)

// NB: Removed code for serialising std::unordered_map, as the latest version of boost includes support.
//     Serialisation format may have changed...

BOOST_CLASS_EXPORT_KEY2(NEWMAT::Matrix, "Matrix")
BOOST_CLASS_EXPORT_KEY2(NEWMAT::SymmetricMatrix, "SymmetricMatrix")
BOOST_CLASS_EXPORT_KEY2(NEWMAT::DiagonalMatrix, "DiagonalMatrix")
BOOST_CLASS_EXPORT_KEY2(NEWMAT::ColumnVector, "ColumnVector")
BOOST_CLASS_EXPORT_KEY2(NEWMAT::RowVector, "RowVector")

BOOST_CLASS_TRACKING(NEWMAT::Matrix, boost::serialization::track_never)
BOOST_CLASS_TRACKING(NEWMAT::SymmetricMatrix, boost::serialization::track_never)
BOOST_CLASS_TRACKING(NEWMAT::DiagonalMatrix, boost::serialization::track_never)
BOOST_CLASS_TRACKING(NEWMAT::ColumnVector, boost::serialization::track_never)
BOOST_CLASS_TRACKING(NEWMAT::RowVector, boost::serialization::track_never)

class SerialisationException : public std::logic_error
{
public:
    SerialisationException(const std::string &what) : logic_error(what) { }
};

namespace boost {
namespace serialization {
    template<class Archive>
    void save(Archive &ar, const arma::Mat<NEWMAT::Real> &m, const unsigned int version)
    {
        std::size_t storage = m.n_elem;
        ar << storage;
        for (std::size_t i = 0; i < storage; i++) {
            ar << m(i);
        }
    }

    template<class Archive>
    void load(Archive &ar, arma::Mat<NEWMAT::Real> &m, const unsigned int version)
    {
        std::size_t storage;

        ar >> storage;

        if (storage != m.n_elem)
            throw SerialisationException("The matrix that is being constructed has a different storage size from the one that was serialised");

        for (std::size_t i = 0; i < storage; i++)
            ar >> m(i);
    }

    template<class Archive>
    void save(Archive &ar, const NEWMAT::Matrix &m, const unsigned int version)
    {
        std::size_t rows = m.Nrows();
        ar << rows;
        std::size_t cols = m.Ncols();
        ar << cols;
        ar << boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }

    template<class Archive>
    void load(Archive &ar, NEWMAT::Matrix &m, const unsigned int version)
    {
        std::size_t rows, cols;
        ar >> rows;
        ar >> cols;
        m.ReSize(rows, cols);
        ar >> boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }

    template<class Archive>
    void save(Archive &ar, const NEWMAT::SymmetricMatrix &m, const unsigned int version)
    {
        std::size_t rows = m.Nrows();
        ar << rows;
        ar << boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }

    template<class Archive>
    void load(Archive &ar, NEWMAT::SymmetricMatrix &m, const unsigned int version)
    {
        std::size_t rows;
        ar >> rows;

        m.ReSize(rows);
        ar >> boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }

    template<class Archive>
    void save(Archive &ar, const NEWMAT::DiagonalMatrix &m, const unsigned int version)
    {
        std::size_t rows = m.Nrows();
        ar << rows;
        ar << boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }

    template<class Archive>
    void load(Archive &ar, NEWMAT::DiagonalMatrix &m, const unsigned int version)
    {
        std::size_t rows;
        ar >> rows;
        m.ReSize(rows);
        ar >> boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }

    template<class Archive>
    void save(Archive &ar, const NEWMAT::ColumnVector &m, const unsigned int version)
    {
        std::size_t rows = m.Nrows();
        ar << rows;
        ar << boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }

    template<class Archive>
    void load(Archive &ar, NEWMAT::ColumnVector &m, const unsigned int version)
    {
        std::size_t rows;
        ar >> rows;
        m.ReSize(rows);
        ar >> boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }

    template<class Archive>
    void save(Archive &ar, const NEWMAT::RowVector &m, const unsigned int version)
    {
        std::size_t cols = m.Ncols();
        ar << cols;
        ar << boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }

    template<class Archive>
    void load(Archive &ar, NEWMAT::RowVector &m, const unsigned int version)
    {
        std::size_t cols;
        ar >> cols;
        m.ReSize(cols);
        ar >> boost::serialization::base_object<arma::Mat<NEWMAT::Real>>(m);
    }
}
}

BOOST_SERIALIZATION_SPLIT_FREE(arma::Mat<NEWMAT::Real>)
BOOST_SERIALIZATION_SPLIT_FREE(NEWMAT::Matrix)
BOOST_SERIALIZATION_SPLIT_FREE(NEWMAT::SymmetricMatrix)
BOOST_SERIALIZATION_SPLIT_FREE(NEWMAT::DiagonalMatrix)
BOOST_SERIALIZATION_SPLIT_FREE(NEWMAT::ColumnVector)
BOOST_SERIALIZATION_SPLIT_FREE(NEWMAT::RowVector)

#endif // SERIALISATION_H
