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

#include "profilepriors.h"
#include "nlopt.hpp"
#include <limits>

BOOST_CLASS_EXPORT_IMPLEMENT(FlatPrior)
BOOST_CLASS_EXPORT_IMPLEMENT(SimpleEdgePrior)
BOOST_CLASS_EXPORT_IMPLEMENT(BlockPrior)
BOOST_CLASS_EXPORT_IMPLEMENT(RightAngledPrior)
BOOST_CLASS_EXPORT_IMPLEMENT(ExponentialPrior)
BOOST_CLASS_EXPORT_IMPLEMENT(Exponential2Prior)
BOOST_CLASS_EXPORT_IMPLEMENT(DoubleExponentialPrior)
BOOST_CLASS_EXPORT_IMPLEMENT(DoubleExponential2Prior)
BOOST_CLASS_EXPORT_IMPLEMENT(ParabolicPrior)

using namespace NEWMAT;

ProfilePrior::~ProfilePrior()
{
}

FlatPrior::FlatPrior(double intensity)
    : m_intensity(intensity)
{
}

ColumnVector FlatPrior::Create(int points, double spacing) const
{
    BOOST_LOG_TRIVIAL(info) << "FlatPrior: Creating profile prior (intensity = " << m_intensity << ")";

    ColumnVector result(points);
    result = m_intensity;

    return result;
}

SimpleEdgePrior::SimpleEdgePrior(double intensitya, double intensityb)
    : m_intensityA(intensitya),
      m_intensityB(intensityb)
{
}

ColumnVector SimpleEdgePrior::Create(int points, double spacing) const
{
    BOOST_LOG_TRIVIAL(info) << "SimpleEdgePrior: Creating profile prior (intensity A = " << m_intensityA
                             << ", intensity B = " << m_intensityB << ")";

    ColumnVector result(points);
    std::size_t half = points / 2;
    result.Rows(1, half) = m_intensityA;
    result.Rows(half + 1, points) = m_intensityB;

    return result;
}

BlockPrior::BlockPrior(double width, double intensitya, double intensityb)
    : m_width(width),
      m_intensityA(intensitya),
      m_intensityB(intensityb)
{
}

ColumnVector BlockPrior::Create(int points, double spacing) const
{
    BOOST_LOG_TRIVIAL(info) << "BlockPrior: Creating profile prior (intensity A = " << m_intensityA
                             << ", intensity B = " << m_intensityB << ", width = " << m_width << ")";

    ColumnVector result(points);
    std::size_t half = points / 2;
    result = m_intensityA;
    result.Rows(half + 1, half + static_cast<int>(m_width / spacing + 0.5)) = m_intensityB;

    return result;
}

RightAngledPrior::RightAngledPrior(double width, double intensitya, double intensityb)
    : m_width(width),
      m_intensityA(intensitya),
      m_intensityB(intensityb)
{
}

ColumnVector RightAngledPrior::Create(int points, double spacing) const
{
    BOOST_LOG_TRIVIAL(info) << "RightAngledPrior: Creating profile prior (intensity A = " << m_intensityA
                             << ", intensity B = " << m_intensityB << ")";

    ColumnVector result(points);
    std::size_t half = points / 2;
    std::size_t width = static_cast<int>(m_width / spacing + 0.5);

    result = m_intensityA;

    for (std::size_t i = 0; i < width; i++)
        result(half + 1 + i) = static_cast<float>(width - i) / width * (m_intensityB - m_intensityA) + m_intensityA;

    return result;
}

ExponentialPrior::ExponentialPrior(double timeconst, double intensitya, double intensityb)
    : m_timeConst(timeconst),
      m_intensityA(intensitya),
      m_intensityB(intensityb)
{
}

ColumnVector ExponentialPrior::Create(int points, double spacing) const
{
    BOOST_LOG_TRIVIAL(info) << "ExponentialPrior: Creating profile prior (intensity A = " << m_intensityA
                             << ", intensity B = " << m_intensityB << ", time constant = " << m_timeConst << ")";

    ColumnVector result(points);
    std::size_t half = points / 2;

    result = m_intensityA;

    for (std::size_t i = 0; i < points - half; i++)
        result(half + 1 + i) = m_intensityA + (m_intensityB - m_intensityA) * std::exp(-static_cast<float>(i) * spacing / m_timeConst);

    return result;
}

Exponential2Prior::Exponential2Prior(double timeconst, double intensitya, double intensityb, double intensityc)
    : m_timeConst(timeconst),
      m_intensityA(intensitya),
      m_intensityB(intensityb),
      m_intensityC(intensityc)
{
}

ColumnVector Exponential2Prior::Create(int points, double spacing) const
{
    BOOST_LOG_TRIVIAL(info) << "Exponential2Prior: Creating profile prior (intensity A = " << m_intensityA
                             << ", intensity B = " << m_intensityB << ", intensity C = " << m_intensityC
                             << ", time constant = " << m_timeConst << ")";

    ColumnVector result(points);
    std::size_t half = points / 2;

    result = m_intensityA;

    for (std::size_t i = 0; i < points - half; i++)
        result(half + 1 + i) = m_intensityC + (m_intensityB - m_intensityC) * std::exp(-static_cast<float>(i) * spacing / m_timeConst);

    return result;
}

DoubleExponentialPrior::DoubleExponentialPrior(double timeconst1, double timeconst2,
                                               double constant, double intensitya, double intensityb)
    : m_timeConst1(timeconst1),
      m_timeConst2(timeconst2),
      m_intensityA(intensitya),
      m_intensityB(intensityb),
      m_constant(constant)
{
    if (timeconst1 >= timeconst2)
        throw std::logic_error("DoubleExponentialPrior: Time constant 1 should be smaller than time constant 2");
}

double DoubleExponentialPrior::Objective(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
    if (!grad.empty())
        throw std::logic_error("DoubleExponentialPrior::Objective is intended to be used with derivative-free optimisers");

    const DoubleExponentialPrior *obj = static_cast<const DoubleExponentialPrior *>(f_data);

    double coef1 = x[0];
    double coef2 = x[1];

    double xminmax = obj->m_timeConst1 * obj->m_timeConst2 / (obj->m_timeConst1 - obj->m_timeConst2)
            * std::log(- obj->m_timeConst1 * coef2 / obj->m_timeConst2 / coef1);

    if (xminmax < 0.0 || xminmax > obj->m_length / 2)
        return std::numeric_limits<double>::infinity();

    double ssq = std::pow(coef1 + coef2 + obj->m_constant - obj->m_intensityA, 2);
    ssq += std::pow(coef1 * std::exp(-xminmax / obj->m_timeConst1) + coef2 * std::exp(-xminmax / obj->m_timeConst2) + obj->m_constant - obj->m_intensityB, 2);

    BOOST_LOG_TRIVIAL(trace) << "xminmax = " << xminmax << ", coef1 = " << coef1 << ", coef2 = " << coef2 << ", ssq = " << ssq;

    return ssq;
}

ColumnVector DoubleExponentialPrior::Create(int points, double spacing) const
{
    BOOST_LOG_TRIVIAL(trace) << "DoubleExponentialPrior: Fitting parameters to requested intensities";

    std::vector<double> lb(2);
    std::vector<double> ub(2);
    std::vector<double> x(2);

    double eps = 1e-10;

    if (m_intensityA < m_intensityB)
    {
        lb[0] = -HUGE_VAL;
        ub[0] = -eps;
        x[0] = -1.0;
        lb[1] = eps;
        ub[1] = HUGE_VAL;
        x[1] = 1.0;
    }
    else
    {
        lb[0] = eps;
        ub[0] = HUGE_VAL;
        x[0] = 1.0;
        lb[1] = -HUGE_VAL;
        ub[1] = -eps;
        x[1] = -1.0;
    }

    m_length = (points - 1) * spacing;

    nlopt::opt optimiser(nlopt::LN_NELDERMEAD, 2);
    optimiser.set_min_objective(&Objective, const_cast<DoubleExponentialPrior *>(this));
    optimiser.set_lower_bounds(lb);
    optimiser.set_upper_bounds(ub);

    double ftol = 1e-7;
    optimiser.set_ftol_rel(ftol);

    double f;
    nlopt::result nlr = optimiser.optimize(x, f);
    double coef1 = x[0];
    double coef2 = x[1];

    if (f > ftol * std::pow(std::abs(m_intensityA) + std::abs(m_intensityB) + std::abs(m_constant), 2))
        throw std::runtime_error("DoubleExponentialPrior: Parameter fit did not converge");

    BOOST_LOG_TRIVIAL(info) << "DoubleExponentialPrior parameter fit successful. NLopt result code was " << nlr << ", parameters:"
                               << "\n\tTime constant 1 = " << m_timeConst1
                               << "\n\tTime constant 2 = " << m_timeConst2
                               << "\n\tCoefficient 1 = " << coef1
                               << "\n\tCoefficient 2 = " << coef2
                               << "\n\tConstant = " << m_constant;

    ColumnVector result(points);
    std::size_t half = points / 2;

    result = m_constant;

    for (std::size_t i = 0; i < points - half; i++)
        result(half + 1 + i) = coef1 * std::exp(-static_cast<float>(i) * spacing / m_timeConst1)
                            + coef2 * std::exp(-static_cast<float>(i) * spacing / m_timeConst2) + m_constant;

    return result;
}

DoubleExponential2Prior::DoubleExponential2Prior(double timeconst1, double timeconst2, double constant,
                                                 double intensitya, double intensityb, double intensityc)
    : m_timeConst1(timeconst1),
      m_timeConst2(timeconst2),
      m_intensityA(intensitya),
      m_intensityB(intensityb),
      m_intensityC(intensityc),
      m_constant(constant)
{
    if (timeconst1 >= timeconst2)
        throw std::logic_error("DoubleExponential2Prior: Time constant 1 should be smaller than time constant 2");
}

double DoubleExponential2Prior::Objective(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
    if (!grad.empty())
        throw std::logic_error("DoubleExponential2Prior::Objective is intended to be used with derivative-free optimisers");

    const DoubleExponential2Prior *obj = static_cast<const DoubleExponential2Prior *>(f_data);

    double coef1 = x[0];
    double coef2 = x[1];

    double xminmax = obj->m_timeConst1 * obj->m_timeConst2 / (obj->m_timeConst1 - obj->m_timeConst2)
            * std::log(- obj->m_timeConst1 * coef2 / obj->m_timeConst2 / coef1);

    if (xminmax < 0.0 || xminmax > obj->m_length / 2)
        return std::numeric_limits<double>::infinity();

    double ssq = std::pow(coef1 + coef2 + obj->m_intensityC - obj->m_intensityA, 2);
    ssq += std::pow(coef1 * std::exp(-xminmax / obj->m_timeConst1) + coef2 * std::exp(-xminmax / obj->m_timeConst2) + obj->m_intensityC - obj->m_intensityB, 2);

    BOOST_LOG_TRIVIAL(trace) << "xminmax = " << xminmax << ", coef1 = " << coef1 << ", coef2 = " << coef2 << ", ssq = " << ssq;

    return ssq;
}

ColumnVector DoubleExponential2Prior::Create(int points, double spacing) const
{
    BOOST_LOG_TRIVIAL(trace) << "DoubleExponential2Prior: Fitting parameters to requested intensities";

    std::vector<double> lb(2);
    std::vector<double> ub(2);
    std::vector<double> x(2);

    double eps = 1e-10;

    if (m_intensityA < m_intensityB)
    {
        lb[0] = -HUGE_VAL;
        ub[0] = -eps;
        x[0] = -1.0;
        lb[1] = eps;
        ub[1] = HUGE_VAL;
        x[1] = 1.0;
    }
    else
    {
        lb[0] = eps;
        ub[0] = HUGE_VAL;
        x[0] = 1.0;
        lb[1] = -HUGE_VAL;
        ub[1] = -eps;
        x[1] = -1.0;
    }

    m_length = (points - 1) * spacing;

    nlopt::opt optimiser(nlopt::LN_NELDERMEAD, 2);
    optimiser.set_min_objective(&Objective, const_cast<DoubleExponential2Prior *>(this));
    optimiser.set_lower_bounds(lb);
    optimiser.set_upper_bounds(ub);

    double ftol = 1e-7;
    optimiser.set_ftol_rel(ftol);

    double f;
    nlopt::result nlr = optimiser.optimize(x, f);
    double coef1 = x[0];
    double coef2 = x[1];

    if (f > ftol * std::pow(std::abs(m_intensityA) + std::abs(m_intensityB) + std::abs(m_intensityC), 2))
        throw std::runtime_error("DoubleExponential2Prior: Parameter fit did not converge");

    BOOST_LOG_TRIVIAL(info) << "DoubleExponential2Prior parameter fit successful. NLopt result code was " << nlr << ", parameters:"
                               << "\n\tTime constant 1 = " << m_timeConst1
                               << "\n\tTime constant 2 = " << m_timeConst2
                               << "\n\tCoefficient 1 = " << coef1
                               << "\n\tCoefficient 2 = " << coef2
                               << "\n\tInner intensity = " << m_constant
                               << "\n\tOuter constant = " << m_intensityC;

    ColumnVector result(points);
    std::size_t half = points / 2;

    result = m_constant;

    for (std::size_t i = 0; i < points - half; i++)
        result(half + 1 + i) = coef1 * std::exp(-static_cast<float>(i) * spacing / m_timeConst1)
                            + coef2 * std::exp(-static_cast<float>(i) * spacing / m_timeConst2) + m_intensityC;

    return result;
}

ParabolicPrior::ParabolicPrior(double minintensity, double maxintensity)
    : m_minIntensity(minintensity),
      m_maxIntensity(maxintensity)
{
}

ColumnVector ParabolicPrior::Create(int points, double spacing) const
{
    BOOST_LOG_TRIVIAL(info) << "ParabolicPrior: Creating profile prior (min intensity = " << m_minIntensity
                             << ", max intensity = " << m_maxIntensity << ")";

    ColumnVector result(points);
    int half = points / 2;
    double coef = (m_maxIntensity - m_minIntensity) / half / half;
    for (int i = 1; i <= points; i++)
        result(i) = m_minIntensity + coef * (i - half) * (i - half);

    return result;
}
