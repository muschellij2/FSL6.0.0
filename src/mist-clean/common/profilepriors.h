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

#ifndef PROFILEPRIORS_H
#define PROFILEPRIORS_H

#include "newmat.h"
#include <boost/log/trivial.hpp>
#include "serialisation.h"

class FlatPrior;
class SimpleEdgePrior;
class BlockPrior;
class RightAngledPrior;
class ExponentialPrior;
class Exponential2Prior;
class DoubleExponentialPrior;
class DoubleExponential2Prior;
class ParabolicPrior;
BOOST_CLASS_EXPORT_KEY(FlatPrior)
BOOST_CLASS_EXPORT_KEY(SimpleEdgePrior)
BOOST_CLASS_EXPORT_KEY(BlockPrior)
BOOST_CLASS_EXPORT_KEY(RightAngledPrior)
BOOST_CLASS_EXPORT_KEY(ExponentialPrior)
BOOST_CLASS_EXPORT_KEY(Exponential2Prior)
BOOST_CLASS_EXPORT_KEY(DoubleExponentialPrior)
BOOST_CLASS_EXPORT_KEY(DoubleExponential2Prior)
BOOST_CLASS_EXPORT_KEY(ParabolicPrior)

namespace boost
{
namespace serialization
{
    template<class Archive>
    inline void save_construct_data(Archive &ar, const FlatPrior *pp, const unsigned int version);

    template<class Archive>
    inline void save_construct_data(Archive &ar, const SimpleEdgePrior *pp, const unsigned int version);

    template<class Archive>
    inline void save_construct_data(Archive &ar, const BlockPrior *pp, const unsigned int version);

    template<class Archive>
    inline void save_construct_data(Archive &ar, const RightAngledPrior *pp, const unsigned int version);

    template<class Archive>
    inline void save_construct_data(Archive &ar, const ExponentialPrior *pp, const unsigned int version);

    template<class Archive>
    inline void save_construct_data(Archive &ar, const Exponential2Prior *pp, const unsigned int version);

    template<class Archive>
    inline void save_construct_data(Archive &ar, const DoubleExponentialPrior *pp, const unsigned int version);
    
    template<class Archive>
    inline void save_construct_data(Archive &ar, const DoubleExponential2Prior *pp, const unsigned int version);

    template<class Archive>
    inline void save_construct_data(Archive &ar, const ParabolicPrior *pp, const unsigned int version);
}
}


class ProfilePrior
{
public:
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const = 0;
    virtual ~ProfilePrior();

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising base ProfilePrior";
    }
};

class FlatPrior : public ProfilePrior
{
public:
    FlatPrior(double intensity);
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const;

private:
    double m_intensity;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising FlatPrior";

        ar & boost::serialization::base_object<ProfilePrior>(*this);
    }

    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const FlatPrior *pp, const unsigned int version);
};

class SimpleEdgePrior : public ProfilePrior
{
public:
    SimpleEdgePrior(double intensitya, double intensityb);
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const;

private:
    double m_intensityA;
    double m_intensityB;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising SimpleEdgePrior";

        ar & boost::serialization::base_object<ProfilePrior>(*this);
    }

    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const SimpleEdgePrior *pp, const unsigned int version);
};

class BlockPrior : public ProfilePrior
{
public:
    BlockPrior(double width, double intensitya, double intensityb);
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const;

private:
    double m_width;
    double m_intensityA;
    double m_intensityB;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising BlockPrior";

        ar & boost::serialization::base_object<ProfilePrior>(*this);
    }

    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const BlockPrior *pp, const unsigned int version);
};

class RightAngledPrior : public ProfilePrior
{
public:
    RightAngledPrior(double width, double intensitya, double intensityb);
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const;

private:
    double m_width;
    double m_intensityA;
    double m_intensityB;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising RightAngledPrior";

        ar & boost::serialization::base_object<ProfilePrior>(*this);
    }

    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const RightAngledPrior *pp, const unsigned int version);
};

class ExponentialPrior : public ProfilePrior
{
public:
    ExponentialPrior(double timeconst, double intensitya, double intensityb);
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const;

private:
    double m_timeConst;
    double m_intensityA;
    double m_intensityB;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising ExponentialPrior";

        ar & boost::serialization::base_object<ProfilePrior>(*this);
    }

    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const ExponentialPrior *pp, const unsigned int version);
};

class Exponential2Prior : public ProfilePrior
{
public:
    Exponential2Prior(double timeconst, double intensitya, double intensityb, double intensityc);
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const;

private:
    double m_timeConst;
    double m_intensityA;
    double m_intensityB;
    double m_intensityC;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising Exponential2Prior";

        ar & boost::serialization::base_object<ProfilePrior>(*this);
    }

    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const Exponential2Prior *pp, const unsigned int version);
};

class DoubleExponentialPrior : public ProfilePrior
{
public:
    DoubleExponentialPrior(double timeconst1, double timeconst2, double constant, double intensitya, double intensityb);

    static double Objective(const std::vector<double> &x, std::vector<double> &grad, void *f_data);
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const;

private:
    double m_timeConst1;
    double m_timeConst2;
    double m_intensityA;
    double m_intensityB;
    double m_constant;

    // Temporary that needs to be passed to Objective() - no need to serialise!
    mutable double m_length;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising DoubleExponentialPrior";

        ar & boost::serialization::base_object<ProfilePrior>(*this);
    }

    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const DoubleExponentialPrior *pp, const unsigned int version);
};

class DoubleExponential2Prior : public ProfilePrior
{
public:
    DoubleExponential2Prior(double timeconst1, double timeconst2, double constant, double intensitya, double intensityb, double intensityc);

    static double Objective(const std::vector<double> &x, std::vector<double> &grad, void *f_data);
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const;

private:
    // NB: Names are confusing but consistent with DoubleExponentialPrior; m_constant is 'inside intensity'
    double m_timeConst1;
    double m_timeConst2;
    double m_intensityA;
    double m_intensityB;
    double m_intensityC;
    double m_constant;

    // Temporary that needs to be passed to Objective() - no need to serialise!
    mutable double m_length;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising DoubleExponential2Prior";

        ar & boost::serialization::base_object<ProfilePrior>(*this);
    }

    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const DoubleExponential2Prior *pp, const unsigned int version);
};

class ParabolicPrior : public ProfilePrior
{
public:
    ParabolicPrior(double minintensity, double maxintensity);
    virtual NEWMAT::ColumnVector Create(int points, double spacing) const;

private:
    double m_minIntensity;
    double m_maxIntensity;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "(De)serialising ParabolicPrior";

        ar & boost::serialization::base_object<ProfilePrior>(*this);
    }

    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const ParabolicPrior *pp, const unsigned int version);
};


namespace boost
{
namespace serialization
{
    template<class Archive>
    inline void save_construct_data(Archive &ar, const FlatPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for FlatPrior";

        ar << pp->m_intensity;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, FlatPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for FlatPrior";

        double intensity;
        ar >> intensity;

        ::new(pp)FlatPrior(intensity);
    }

    template<class Archive>
    inline void save_construct_data(Archive &ar, const SimpleEdgePrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for SimpleEdgePrior";

        ar << pp->m_intensityA;
        ar << pp->m_intensityB;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, SimpleEdgePrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for SimpleEdgePrior";

        double intensitya, intensityb;
        ar >> intensitya;
        ar >> intensityb;

        ::new(pp)SimpleEdgePrior(intensitya, intensityb);
    }

    template<class Archive>
    inline void save_construct_data(Archive &ar, const BlockPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for BlockPrior";

        ar << pp->m_width;
        ar << pp->m_intensityA;
        ar << pp->m_intensityB;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, BlockPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for BlockPrior";

        double width, intensitya, intensityb;
        ar >> width;
        ar >> intensitya;
        ar >> intensityb;

        ::new(pp)BlockPrior(width, intensitya, intensityb);
    }

    template<class Archive>
    inline void save_construct_data(Archive &ar, const RightAngledPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for RightAngledPrior";

        ar << pp->m_width;
        ar << pp->m_intensityA;
        ar << pp->m_intensityB;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, RightAngledPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for RightAngledPrior";

        double width, intensitya, intensityb;
        ar >> width;
        ar >> intensitya;
        ar >> intensityb;

        ::new(pp)RightAngledPrior(width, intensitya, intensityb);
    }

    template<class Archive>
    inline void save_construct_data(Archive &ar, const ExponentialPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for ExponentialPrior";

        ar << pp->m_timeConst;
        ar << pp->m_intensityA;
        ar << pp->m_intensityB;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, ExponentialPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for ExponentialPrior";

        double timeconst, intensitya, intensityb;
        ar >> timeconst;
        ar >> intensitya;
        ar >> intensityb;

        ::new(pp)ExponentialPrior(timeconst, intensitya, intensityb);
    }

    template<class Archive>
    inline void save_construct_data(Archive &ar, const Exponential2Prior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for Exponential2Prior";

        ar << pp->m_timeConst;
        ar << pp->m_intensityA;
        ar << pp->m_intensityB;
        ar << pp->m_intensityC;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, Exponential2Prior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for Exponential2Prior";

        double timeconst, intensitya, intensityb, intensityc;
        ar >> timeconst;
        ar >> intensitya;
        ar >> intensityb;
        ar >> intensityc;

        ::new(pp)Exponential2Prior(timeconst, intensitya, intensityb, intensityc);
    }

    template<class Archive>
    inline void save_construct_data(Archive &ar, const DoubleExponentialPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for DoubleExponentialPrior";

        ar << pp->m_timeConst1;
        ar << pp->m_timeConst2;
        ar << pp->m_intensityA;
        ar << pp->m_intensityB;
        ar << pp->m_constant;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, DoubleExponentialPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for DoubleExponentialPrior";

        double timeconst1, timeconst2, constant, intensitya, intensityb;
        ar >> timeconst1;
        ar >> timeconst2;
        ar >> intensitya;
        ar >> intensityb;
        ar >> constant;

        ::new(pp)DoubleExponentialPrior(timeconst1, timeconst2, intensitya, intensityb, constant);
    }

    template<class Archive>
    inline void save_construct_data(Archive &ar, const DoubleExponential2Prior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for DoubleExponential2Prior";

        ar << pp->m_timeConst1;
        ar << pp->m_timeConst2;
        ar << pp->m_intensityA;
        ar << pp->m_intensityB;
        ar << pp->m_intensityC;
        ar << pp->m_constant;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, DoubleExponential2Prior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for DoubleExponential2Prior";

        double timeconst1, timeconst2, constant, intensitya, intensityb, intensityc;
        ar >> timeconst1;
        ar >> timeconst2;
        ar >> intensitya;
        ar >> intensityb;
        ar >> intensityc;
        ar >> constant;

        ::new(pp)DoubleExponential2Prior(timeconst1, timeconst2, intensitya, intensityb, intensityc, constant);
    }

    template<class Archive>
    inline void save_construct_data(Archive &ar, const ParabolicPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for ParabolicPrior";

        ar << pp->m_minIntensity;
        ar << pp->m_maxIntensity;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, ParabolicPrior *pp, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for ParabolicPrior";

        double minintensity, maxintensity;
        ar >> minintensity;
        ar >> maxintensity;

        ::new(pp)ParabolicPrior(minintensity, maxintensity);
    }
}
}

#endif
