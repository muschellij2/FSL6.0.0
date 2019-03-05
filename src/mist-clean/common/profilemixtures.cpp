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

#include "profilemixtures.h"
#include "stats.h"
#include "nlopt.hpp"
//#include "optimisers.h"
#include <fstream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <random>
#include <numeric>
#include <boost/make_shared.hpp>
#include <boost/log/trivial.hpp>

#define WANT_STREAM
#include "newmatio.h"

BOOST_CLASS_EXPORT_IMPLEMENT(ProfileMixtures)

using namespace NEWMAT;

ProfileMixtures::ProfileMixtures(const std::vector<string> &modalitynames, int reflength,
                                 int datalength, int components)
    : ProfileModel(modalitynames, reflength, datalength),
    m_components(components)
{
    BOOST_LOG_TRIVIAL(debug) << "Initialising mixture model";

    // NB: These are not the actual default values, just some values to make sure the object has valid state

    int steps = GetNumberOfSteps();
    m_deltaMean = steps / 2.0;
    m_deltaPrecision = 1.0;
    m_mixingAlpha = std::vector<double>(m_components, 1.0);
    // This is the actual minimum! (Or .5 lower..?)
    m_alpha = (m_refLength - 1) / 2.0 + 1.0;
    m_n0 = 1.0;

    m_ftol = 1.0;
    m_ftolIgnore = 0;
    m_maxEval = 1;
    m_maxRetries = 0;

    m_useMCMC = false;
    m_mcmcIterations = 500000;
    m_mcmcBurnIn = 100000;
    m_mcmcThin = 4000;
    m_mcmcScaleMixing = 0.01;
    m_mcmcScaleMeans = 0.1;
    m_mcmcScaleCovarianceCoefs = 0.1;

    for (auto mn : modalitynames)
    {
        ModalityParameters p;

        for (int i = 0; i < components; i++)
        {
            ColumnVector mean(m_refLength);
            mean = 0.0;
            p.componentPriorMeans.push_back(mean);
            ColumnVector covcoefs(m_refLength);
            covcoefs = 1.0;
            p.componentPriorCovarianceCoefs.push_back(covcoefs);
        }

        m_modalityParameters[mn] = p;
        SetSmoothness(mn, 1.0);
    }
}

boost::shared_ptr<ProfileModel> ProfileMixtures::Clone() const
{
    return boost::make_shared<ProfileMixtures>(*this);
}

void ProfileMixtures::SetComponentPriorMean(const std::string &modality, int component, boost::shared_ptr<const ProfilePrior> mean, double spacing)
{
    if (m_modalityParameters.find(modality) == m_modalityParameters.end())
        throw ModelException("Specified modality does not exist in model");

    ColumnVector unfilteredmean = mean->Create(m_refLength, spacing);

    // TODO: Clean up
    ColumnVector filteredmean = unfilteredmean;

    double normfact = 0.0;
    if (unfilteredmean.Sum() > 0.0)
        normfact = unfilteredmean.Sum() / filteredmean.Sum();
    else
        normfact = 1.0;

    m_modalityParameters[modality].componentPriorMeans[component] = normfact * filteredmean;
}

void ProfileMixtures::SetComponentPriorCovarianceCoefs(const std::string &modality, int component, boost::shared_ptr<const ProfilePrior> coefs, double spacing)
{
    if (m_modalityParameters.find(modality) == m_modalityParameters.end())
        throw ModelException("Specified modality does not exist in model");

    m_modalityParameters[modality].componentPriorCovarianceCoefs[component] = coefs->Create(m_refLength, spacing);
}

void ProfileMixtures::SetSmoothness(const string &modality, double smoothness)
{
    if (m_modalityParameters.find(modality) == m_modalityParameters.end())
        throw ModelException("Specified modality does not exist in model");

    m_modalityParameters[modality].smoothness = smoothness;

    m_modalityParameters[modality].verySmallGpinv.clear();
    int steps = GetNumberOfSteps();
    int fitdatalength = m_dataLength;

    BOOST_LOG_TRIVIAL(debug) << "SetSmoothness(): steps = " << steps << " for " << modality;

    Matrix G = ProfileMixtureProbability::GetSmoothingMatrix(smoothness * smoothness, m_refLength);
    m_modalityParameters[modality].G = G;

    int rank = 0;
    for (std::size_t delta = 0; delta < steps; delta++)
        m_modalityParameters[modality].verySmallGpinv.push_back(
                Stats::pinv(G.Columns(delta + 1, delta + fitdatalength), 1e-5, rank).t());
}

void ProfileMixtures::SetN0(int n0)
{
    m_n0 = n0;
}

void ProfileMixtures::SetAlpha(double alpha)
{
    double alphamin = (m_refLength - 1) / 2.0 + 1.0;

    if (alpha < alphamin)
        throw ModelException(std::string("Alpha should be at least ") + std::to_string(alphamin));

    m_alpha = alpha;
}

void ProfileMixtures::SetDeltaStdev(double stdev)
{
    m_deltaPrecision = 1.0 / std::pow(stdev, 2);
}

void ProfileMixtures::SetMixingAlpha(const std::vector<double> &mixingalpha)
{
    m_mixingAlpha = mixingalpha;
}

void ProfileMixtures::SetFTolerance(double ftol)
{
    m_ftol = ftol;
}

void ProfileMixtures::SetFToleranceIgnoreCount(double count)
{
    m_ftolIgnore = count;
}

void ProfileMixtures::SetMaxEvaluations(int maxeval)
{
    m_maxEval = maxeval;
}

void ProfileMixtures::SetMaxRetries(int maxretries)
{
    m_maxRetries = maxretries;
}

void ProfileMixtures::SetUseMCMC(bool usemcmc)
{
    m_useMCMC = usemcmc;
}

int ProfileMixtures::GetNumberOfComponents() const
{
    return m_components;
}

std::unordered_map<std::string, std::vector<ColumnVector> > ProfileMixtures::GetComponentPriorMeans() const
{
    std::unordered_map<std::string, std::vector<ColumnVector> > result;
    for (auto &p : m_modalityParameters)
        result[p.first] = p.second.componentPriorMeans;

    return result;
}

std::unordered_map<std::string, std::vector<ColumnVector> > ProfileMixtures::GetComponentPriorCovarianceCoefs() const
{
    std::unordered_map<std::string, std::vector<ColumnVector> > result;
    for (auto &p : m_modalityParameters)
        result[p.first] = p.second.componentPriorCovarianceCoefs;

    return result;
}

std::vector<double> ProfileMixtures::FullMixingCoefs(const std::vector<double> &coefs)
{
    std::vector<double> v;
    double s = 0.0;

    for (auto &t : coefs)
    {
        s += t;
        v.push_back(t);
    }

    v.push_back(1.0 - s);

    return v;
}

std::unordered_map<std::string, std::vector<double> > ProfileMixtures::GetFullMixingCoefs() const
{
    if (m_modalityParameters.cbegin()->second.componentMeans.size() == 0)
        throw ModelException("Cannot get parameters from untrained model");

    std::unordered_map<std::string, std::vector<double> > result;
    for (auto &p : m_modalityParameters)
        result[p.first] = FullMixingCoefs(p.second.mixingCoefs);

    return result;
}

std::unordered_map<std::string, std::vector<ColumnVector> > ProfileMixtures::GetComponentMeans() const
{
    if (m_modalityParameters.cbegin()->second.componentMeans.size() == 0)
        throw ModelException("Cannot get parameters from untrained model");

    std::unordered_map<std::string, std::vector<ColumnVector> > result;
    for (auto &p : m_modalityParameters)
        result[p.first] = p.second.componentMeans;

    return result;
}

std::unordered_map<std::string, std::vector<ColumnVector> > ProfileMixtures::GetComponentCovarianceCoefs() const
{
    if (m_modalityParameters.cbegin()->second.componentMeans.size() == 0)
        throw ModelException("Cannot get parameters from untrained model");

    std::unordered_map<std::string, std::vector<ColumnVector> > result;
    for (auto &p : m_modalityParameters)
        result[p.first] = p.second.componentCovarianceCoefs;

    return result;
}

std::unordered_map<std::string, std::vector<std::vector<double> > > ProfileMixtures::GetFullMixingCoefsSamples() const
{
    if (m_modalityParameters.cbegin()->second.componentMeans.size() == 0)
        throw ModelException("Cannot get samples from untrained model");

    std::unordered_map<std::string, std::vector<std::vector<double> > > result;
    for (auto &p : m_modalityParameters)
        for (auto &s : p.second.mixingCoefsSamples)
            result[p.first].push_back(FullMixingCoefs(s));

    return result;
}

std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > ProfileMixtures::GetComponentMeansSamples() const
{
    if (m_modalityParameters.cbegin()->second.componentMeans.size() == 0)
        throw ModelException("Cannot get samples from untrained model");

    std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > result;
    for (const auto &p : m_modalityParameters)
        for (int r = 0; r < m_components; r++)
            result[p.first].push_back(p.second.componentMeansSamples[r]);

    return result;
}

std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > ProfileMixtures::GetComponentCovarianceCoefsSamples() const
{
    if (m_modalityParameters.cbegin()->second.componentMeans.size() == 0)
        throw ModelException("Cannot get samples from untrained model");

    std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > result;
    for (auto &p : m_modalityParameters)
        for (int r = 0; r < m_components; r++)
            result[p.first].push_back(p.second.componentCovarianceCoefsSamples[r]);

    return result;
}

std::unordered_map<std::string, double> ProfileMixtures::GetSmoothness() const
{
    if (m_modalityParameters.cbegin()->second.componentMeans.size() == 0)
        throw ModelException("Cannot get parameters from untrained model");

    std::unordered_map<std::string, double> result;
    for (auto &p : m_modalityParameters)
        result[p.first] = p.second.smoothness;

    return result;
}

std::vector<double> ProfileMixtures::Optimise(ProfileMixtureProbability &pmp)
{

    BOOST_LOG_TRIVIAL(info) << "Creating NLopt object:\n"
                                << "\tftol = " << m_ftol << "\n"
                                << "\tftol ignore count = " << m_ftolIgnore << "\n"
                                << "\tmaxeval = " << m_maxEval << "\n";

    nlopt::opt optimiser(nlopt::LD_MMA, pmp.LowerBounds().size());
    optimiser.set_max_objective(&ProfileMixtureProbability::Objective, &pmp);
    optimiser.set_ftol_rel(m_ftol);
    optimiser.set_maxeval(m_maxEval);
    optimiser.set_lower_bounds(pmp.LowerBounds());
    optimiser.set_upper_bounds(pmp.UpperBounds());
    // Is zero tolerance is ok?
    optimiser.add_inequality_mconstraint(&ProfileMixtureProbability::MixingInequality,
                                         &pmp, std::vector<double>(m_modalityParameters.size(), 0.0));

    BOOST_LOG_TRIVIAL(info) << "Starting optimisation";

    std::vector<double> x;
    nlopt::result result = nlopt::FAILURE;

    for (int attempt = 0; ; attempt++)
    {
        x = pmp.GetInitialGuess();

        try
        {
            double f;

            for (int ignore = 0; ignore <= m_ftolIgnore; ignore++)
            {
                result = optimiser.optimize(x, f);

                if (result == nlopt::FTOL_REACHED)
                {
                    if (ignore != m_ftolIgnore)
                        BOOST_LOG_TRIVIAL(info) << "FTol reached; restarting optimisation (will restart "
                                                << m_ftolIgnore - ignore - 1 << " more times)";
                }
                else
                    break;
            }
        }
        catch (std::runtime_error &e)
        {
            BOOST_LOG_TRIVIAL(warning) << "NLopt failed. The exception text was: " << e.what();

            if (attempt == m_maxRetries)
                throw ModelException("Maximum number of optimisation attempts exceeded");
            else
            {
                BOOST_LOG_TRIVIAL(info) << "Retrying ...";
                continue;
            }
        }

        break;
    }

    BOOST_LOG_TRIVIAL(info) << "Optimisation finished successfully with result code " << result << "; unpacking result";

    return x;
}

std::vector<std::vector<double> > ProfileMixtures::GenerateSamples(ProfileMixtureProbability &pmp,
                                    std::vector<double> current)
{
    // NOTE: Local maxima are probably not much of an issue for MCMC as the initial minimisation seems to always converge to the
    // desired optimum

    BOOST_LOG_TRIVIAL(info) << "Starting MCMC";

    // NB: Using fixed seed - might want to make this configurable
    std::mt19937 mt(0);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // NB: Need to get these from cost function object because of normalisation
    auto priorcovcoefs = pmp.GetProfilePriorCovarianceCoefs();

    std::vector<double> lowerbounds = pmp.LowerBounds();
    std::vector<double> upperbounds = pmp.UpperBounds();

    pmp.UnpackPoint(current);
    double lpcurr = pmp.LogProbability();

    std::vector<std::vector<double> > chain;
    int accepted = 0;

    for (int iter = 0; iter < m_mcmcIterations; iter++)
    {
        std::vector<double> proposal(current.size());

        for (auto p = proposal.begin(); p != proposal.end(); p++)
            *p = normal(mt);

        std::unordered_map<std::string, std::vector<double> > mixing;
        std::unordered_map<std::string, std::vector<ColumnVector> > means;
        std::unordered_map<std::string, std::vector<ColumnVector> > covcoefs;

        pmp.UnpackX(proposal, mixing, means, covcoefs);

        for (auto &t : mixing)
        {
            for (int r = 0; r < m_components - 1; r++)
                t.second[r] *= m_mcmcScaleMixing;

            for (int r = 0; r < m_components; r++)
            {
                // We don't need the actual cholesky decomposition to generate MVN samples as sigma = G * coefs * G
                double priorcovcoefsmax = priorcovcoefs[t.first][r].MaximumAbsoluteValue();
                means[t.first][r] = m_mcmcScaleMeans * std::sqrt(priorcovcoefsmax)
                                            * m_modalityParameters[t.first].G * means[t.first][r];
                covcoefs[t.first][r] *= m_mcmcScaleCovarianceCoefs * priorcovcoefsmax;
            }
        }

        proposal = pmp.PackX(mixing, means, covcoefs);

        for (auto p = proposal.begin(), c = current.begin(); p != proposal.end(); )
            *p++ += *c++;

        bool boundscheck = true;
        for (auto p = proposal.cbegin(), lb = lowerbounds.cbegin(), ub = upperbounds.cbegin(); p != proposal.cend(); p++, lb++, ub++)
            if (*p < *lb || *p > *ub)
                boundscheck = false;

        std::vector<double> constraint(m_modalityNames.size());
        pmp.MixingInequality(constraint.size(), constraint.data(), current.size(), proposal.data(), NULL, &pmp);

        if (boundscheck && std::accumulate(constraint.cbegin(), constraint.cend(), 0.0) <= 0.0)
        {
            pmp.UnpackPoint(proposal);
            double lpprop = pmp.LogProbability();
            double lpdiff = lpprop - lpcurr;

            BOOST_LOG_TRIVIAL(debug) << "lpcurr = " << lpcurr << ", lpdiff = " << lpdiff << ", exp(lpdiff) = " << std::exp(lpdiff);

            if (lpdiff > 0.0 || uniform(mt) <= std::exp(lpdiff))
            {
                BOOST_LOG_TRIVIAL(trace) << " -> accept";
                accepted++;
                current = proposal;
                lpcurr = lpprop;
            }
            else
                BOOST_LOG_TRIVIAL(trace) << "-> reject";
        }
        else
            BOOST_LOG_TRIVIAL(trace) << " -> fails constraints";

        if (iter % m_mcmcThin == 0)
        {
            if (iter >= m_mcmcBurnIn)
            {
                BOOST_LOG_TRIVIAL(info) << "MH Iteration " << iter << " (log p = " << lpcurr << ")";
                chain.push_back(current);
            }
            else
                BOOST_LOG_TRIVIAL(info) << "MH Iteration " << iter << " (burn-in, log p = " << lpcurr << ")";
        }
    }

    BOOST_LOG_TRIVIAL(info) << "MCMC acceptance fraction was " << static_cast<float>(accepted) / m_mcmcIterations;

    return chain;
}

void ProfileMixtures::Train(const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &trainingdata)
{
    BOOST_LOG_TRIVIAL(debug) << "Training mixture model";

    ModelException maperror("Training data should be provided for all modalities that were included when constructing the model");

    if (trainingdata.size() != m_modalityParameters.size())
        throw maperror;

    for (auto t : trainingdata)
        if (m_modalityParameters.find(t.first) == m_modalityParameters.end())
            throw maperror;

    std::unordered_map<std::string, std::vector<ColumnVector> > priormeans;
    std::unordered_map<std::string, std::vector<ColumnVector> > priorcovcoefs;
    std::unordered_map<std::string, double> stdevsmooth;
    std::unordered_map<std::string, std::vector<ColumnVector> > scaledtd;

    BOOST_LOG_TRIVIAL(debug) << "Normalising profiles";

    for (const auto &t : m_modalityParameters)
    {
        for (int r = 0; r < m_components; r++)
        {
            priormeans[t.first].push_back((t.second.componentPriorMeans[r] - m_globalMean[t.first]) / m_globalStdev[t.first]);
            priorcovcoefs[t.first].push_back(t.second.componentPriorCovarianceCoefs[r] / std::pow(m_globalStdev[t.first], 2));
        }

        for (const auto &v : trainingdata.find(t.first)->second)
            scaledtd[t.first].push_back((v - m_globalMean[t.first]) / m_globalStdev[t.first]);

        stdevsmooth[t.first] = t.second.smoothness;
    }

    BOOST_LOG_TRIVIAL(debug) << "Creating cost function object";
    ProfileMixtureProbability pmp(scaledtd, priormeans, priorcovcoefs, m_n0, m_alpha, m_deltaMean,
                                  m_deltaPrecision, m_mixingAlpha, stdevsmooth);

    std::vector<double> x = Optimise(pmp);

    std::unordered_map<std::string, std::vector<double> > mixing;
    std::unordered_map<std::string, std::vector<ColumnVector> > means;
    std::unordered_map<std::string, std::vector<ColumnVector> > covcoefs;
    pmp.UnpackX(x, mixing, means, covcoefs);

    for (auto &t : m_modalityParameters)
    {
        t.second.mixingCoefs = mixing[t.first];

        t.second.componentMeans.clear();
        t.second.componentCovarianceCoefs.clear();

        for (int r = 0; r < m_components; r++)
        {
            t.second.componentMeans.push_back(means[t.first][r] * m_globalStdev[t.first] + m_globalMean[t.first]);
            t.second.componentCovarianceCoefs.push_back(covcoefs[t.first][r] * std::pow(m_globalStdev[t.first], 2));
        }
    }

    for (auto &t : m_modalityParameters)
    {
        t.second.mixingCoefsSamples.clear();
        t.second.componentMeansSamples = std::vector<std::vector<ColumnVector> >(m_components);
        t.second.componentCovarianceCoefsSamples = std::vector<std::vector<ColumnVector> >(m_components);
    }

    if (m_useMCMC)
    {
        std::vector<std::vector<double> > samples = GenerateSamples(pmp, x);

        for (const auto &sample : samples)
        {
            std::unordered_map<std::string, std::vector<double> > mixing;
            std::unordered_map<std::string, std::vector<ColumnVector> > means;
            std::unordered_map<std::string, std::vector<ColumnVector> > covcoefs;

            pmp.UnpackX(sample, mixing, means, covcoefs);

            for (auto &t : m_modalityParameters)
            {
                t.second.mixingCoefsSamples.push_back(mixing[t.first]);

                for (int r = 0; r < m_components; r++)
                {
                    t.second.componentMeansSamples[r].push_back(means[t.first][r] * m_globalStdev[t.first] + m_globalMean[t.first]);
                    t.second.componentCovarianceCoefsSamples[r].push_back(covcoefs[t.first][r] * std::pow(m_globalStdev[t.first], 2));
                }
            }
        }
    }
}

std::vector<double> ProfileMixtures::GetDeltaLikelihoods(const std::unordered_map<std::string, ColumnVector> &data,
                                                         bool usedeltaprior) const
{
    // TODO: Log final result to debug

    int steps = GetNumberOfSteps();
    int fitdatalength = m_dataLength;

    if (m_modalityParameters.cbegin()->second.componentMeans.size() == 0)
        throw ModelException("Cannot fit new data to untrained profile model");

    for (const auto &t : data)
    {
        if (t.second.Nrows() != fitdatalength)
            throw ModelException("One or more profiles have incorrect length");

        if (m_modalityParameters.find(t.first) == m_modalityParameters.end())
            throw ModelException("Cannot fit a modality that was not in the training data");
    }

    std::unordered_map<std::string, std::vector<std::vector<double> > > mixingsamples;
    std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > meanssamples;
    std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > covcoefssamples;

    if (m_useMCMC)
    {
        if (m_modalityParameters.cbegin()->second.mixingCoefsSamples.empty())
            throw ModelException("Use of MCMC samples was requested, but the model contains no saved samples");

        BOOST_LOG_TRIVIAL(debug) << "GetDeltaLikelihoods(): Using MCMC samples";

        // Use MCMC samples
        for (const auto &t : m_modalityParameters)
        {
            mixingsamples[t.first] = t.second.mixingCoefsSamples;
            meanssamples[t.first] = t.second.componentMeansSamples;
            covcoefssamples[t.first] = t.second.componentCovarianceCoefsSamples;
        }
    }
    else
    {
        BOOST_LOG_TRIVIAL(debug) << "GetDeltaLikelihoods(): Using MAP estimates";

        // Use MAP estimates
        for (const auto &t : m_modalityParameters)
        {
            mixingsamples[t.first].push_back(t.second.mixingCoefs);

            for (int r = 0; r < m_components; r++)
            {
                meanssamples[t.first].push_back(
                            std::vector<ColumnVector>(1, t.second.componentMeans[r]));
                covcoefssamples[t.first].push_back(
                            std::vector<ColumnVector>(1, t.second.componentCovarianceCoefs[r]));
            }
        }
    }

    std::vector<double> logprobs(steps, -std::numeric_limits<double>::infinity());

    for (int i = 0; i < mixingsamples.cbegin()->second.size(); i++)
    {
        for (int delta = 0; delta < steps; delta++)
        {
            std::string tracestr = std::string("Sample ") + std::to_string(i)
                                               + ": logprobs[" + std::to_string(delta) + "] = ";
            bool tracefirstmodality = true;

            double lp = 0.0;

            // This is a replacement for the shape model during fitting
            if (usedeltaprior)
                lp += Stats::lognormal(delta, m_deltaMean, m_deltaPrecision);

            tracestr += std::to_string(lp) + " + ";

            for (const auto &t : data)
            {
                tracestr += tracefirstmodality ? "log(" : "+ log(";
                tracefirstmodality = false;

                std::vector<double> fullmixcoefs = FullMixingCoefs(mixingsamples[t.first][i]);
                double lse_r = -std::numeric_limits<double>::infinity();

                for (int r = 0; r < m_components; r++)
                {
                    ColumnVector mu = meanssamples[t.first][r][i].Rows(delta + 1, delta + fitdatalength);

                    const Matrix &Gsubinv = m_modalityParameters.find(t.first)->second.verySmallGpinv[delta];
                    DiagonalMatrix D(covcoefssamples[t.first][r][i].Nrows());
                    for (std::size_t j = 1; j <= covcoefssamples[t.first][r][i].Nrows(); j++)
                        D(j) = 1.0 / covcoefssamples[t.first][r][i](j);

                    Matrix lambda = Gsubinv.t() * D * Gsubinv;

                    double lmvn = Stats::MVN(mu, lambda).LogP(t.second);
                    lse_r = Stats::logsumexp(lse_r, lmvn + std::log(fullmixcoefs[r]));

                    tracestr += std::to_string(fullmixcoefs[r]) + " * exp(" + std::to_string(lmvn) + ")";
                    tracestr += r < m_components - 1 ? " + " : ")[" + t.first + "] ";
                }

                lp += lse_r;
            }

            tracestr += "= " + std::to_string(lp);
            BOOST_LOG_TRIVIAL(debug) << tracestr;

            logprobs[delta] = Stats::logsumexp(logprobs[delta], lp);
        }
    }

    return logprobs;
}

void ProfileMixtures::WritePlots(const Plotting::plotfunc &pfunc, int vertex,
                                 std::unordered_map<std::string, std::vector<ColumnVector> > *profiles) const
{
    for (auto t : GetComponentMeans())
        pfunc(t.second, std::vector<std::size_t>(), 0, vertex, t.first + "_means");

    for (auto t : GetComponentCovarianceCoefs())
    {
        pfunc(t.second, std::vector<std::size_t>(), 0, vertex, t.first + "_covcoefs");

        std::vector<ColumnVector> stdev;
        for (int r = 0; r < t.second.size(); r++)
        {
            const Matrix &G = m_modalityParameters.find(t.first)->second.G;
            Matrix C = G * t.second[r].AsDiagonal() * G;
            ColumnVector s(m_refLength);
            for (std::size_t i = 1; i <= m_refLength; i++)
                s(i) = std::sqrt(C(i, i));

            stdev.push_back(s);
        }

        pfunc(stdev, std::vector<std::size_t>(), 0, vertex, t.first + "_stdev");
    }

    if (profiles)
    {
        int steps = GetNumberOfSteps();
        std::vector<std::size_t> deltas;

        for (std::size_t i = 0; i < profiles->cbegin()->second.size(); i++)
        {
            std::unordered_map<std::string, ColumnVector> fitdata;

            for (const auto &t : *profiles)
            {
                fitdata[t.first] = t.second[i];
            }

            std::vector<double> logprobs = GetDeltaLikelihoods(fitdata, true);
            std::size_t delta = std::max_element(logprobs.cbegin(), logprobs.cend()) - logprobs.cbegin();

            deltas.push_back(delta);
        }

        for (const auto &t : *profiles)
        {
            pfunc(t.second, std::vector<std::size_t>(), 0, vertex, t.first + "_unaligned");
            pfunc(t.second, deltas, steps, vertex, t.first + "_aligned");
        }
    }
}

// -----------------------------------------------------------------------------------------------------------------------------

ProfileMixtureProbability::ProfileMixtureProbability(const std::unordered_map<std::string, std::vector<ColumnVector> > &data,
                            const std::unordered_map<std::string, std::vector<ColumnVector> > &profilepriormeans,
                            const std::unordered_map<std::string, std::vector<ColumnVector> > &profilepriorcovcoefs,
                            int profilen0, int profilealpha, double deltamean, double deltaprecision,
                            const std::vector<double> &mixingalpha, const std::unordered_map<std::string, double> &stdevsmooth) :
    m_data(data),
    m_deltaMean(deltamean),
    m_deltaPrecision(deltaprecision),
    m_mixingAlpha(mixingalpha),
    m_profileN0(profilen0),
    m_profileAlpha(profilealpha),
    m_profilePriorMeans(profilepriormeans),
    m_profilePriorCovarianceCoefs(profilepriorcovcoefs)
{
    BOOST_LOG_TRIVIAL(debug) << "Initialising cost function object";
    m_subjects = m_data.size();
    m_dataLength = m_data.begin()->second[0].Nrows();
    m_components = profilepriormeans.begin()->second.size();
    m_refLength = profilepriormeans.begin()->second[0].Nrows();

    m_steps = m_refLength - m_dataLength;

    for (const auto &t : stdevsmooth)
    {
        BOOST_LOG_TRIVIAL(debug) << " -> Smoothness for modality " << t.first << " is " << t.second;
        // NB: G * diagonal * G.t() is positive definite as required (as sqrt(diagonal) is still diagonal)
        double sigmasq = t.second * t.second;

        Matrix G = GetSmoothingMatrix(sigmasq, m_refLength);

        m_fullG[t.first] = G;
        m_fullGinv[t.first] = G.i();

        // We need all of these because G is not cyclical
        for (int delta = 0; delta < (m_refLength - m_dataLength); delta++)
        {
            Matrix sg = G.Rows(1 + delta, m_dataLength + delta);
            m_smallG[t.first].push_back(sg);
            int rank = 0;
            m_smallGpinv[t.first].push_back(Stats::pinv(sg.t(), 1e-5, rank).t());
        }

        for (int r = 0; r < m_components; r++)
            m_profileBetas[t.first].push_back(m_profileAlpha * ExpandCovariance(t.first, profilepriorcovcoefs.find(t.first)->second[r]));
    }

    std::string info = std::string("ProfileMixtures cost object initialised. Parameters:\n")
            + "\tReference length in points = " + std::to_string(m_refLength) + "\n"
            + "\tSampled profile length in points = " + std::to_string(m_dataLength) + "\n"
            + "\tN0 = " + std::to_string(m_profileN0) + "\n"
            + "\tAlpha = " + std::to_string(m_profileAlpha) + "\n"
            + "\tDelta prior mean = " + std::to_string(m_deltaMean) + "\n"
            + "\tDelta prior precision = " + std::to_string(m_deltaPrecision) + "\n"
            + "\tMixing Dirichlet parameters =";

    for (auto v : m_mixingAlpha)
        info += " " + std::to_string(v);

    info += "\n";

    for (const auto &t : data)
    {
        for (int i = 0; i < m_components; i++)
        {
            info += "\tModality " + t.first + ", component " + std::to_string(i) + " (normalised):\n"
                    + "\t\tPrior mean: " + std::to_string(m_profilePriorMeans[t.first][i](1))
                    + " ... " + std::to_string(m_profilePriorMeans[t.first][i](m_refLength / 2))
                    + " " + std::to_string(m_profilePriorMeans[t.first][i](m_refLength / 2 + 1))
                    + " ... " + std::to_string(m_profilePriorMeans[t.first][i](m_refLength)) + "\n"
                    + "\t\tPrior beta diagonal: " + std::to_string(m_profileBetas[t.first][i](1, 1))
                    + " ... " + std::to_string(m_profileBetas[t.first][i](m_refLength / 2, m_refLength / 2))
                    + " " + std::to_string(m_profileBetas[t.first][i](m_refLength / 2 + 1, m_refLength / 2 + 1))
                    + " ... " + std::to_string(m_profileBetas[t.first][i](m_refLength, m_refLength)) + "\n";
        }
    }

    BOOST_LOG_TRIVIAL(info) << info;
}

Matrix ProfileMixtureProbability::GetSmoothingMatrix(double sigmasq, std::size_t size)
{
    Matrix G(size, size);

    for (int i = 1; i <= size; i++)
    {
        for (int j = 1; j <= size; j++)
        {
            int d = i - j;
            G(i, j) = std::exp(- d * d / 2.0 / sigmasq);
        }
    }

    return G;
}

int ProfileMixtureProbability::GetProfileN0() const
{
    return m_profileN0;
}

double ProfileMixtureProbability::GetProfileAlpha() const
{
    return m_profileAlpha;
}

std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > ProfileMixtureProbability::GetProfilePriorMeans() const
{
    return m_profilePriorMeans;
}

std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > ProfileMixtureProbability::GetProfilePriorCovarianceCoefs() const
{
    return m_profilePriorCovarianceCoefs;
}

Matrix ProfileMixtureProbability::ExpandPrecision(const std::string &modality, const ColumnVector &coefs)
{
    Matrix &Ginv = m_fullGinv[modality];
    DiagonalMatrix D(Ginv.Nrows());
    D = 0.0;

    for (int i = 1; i <= Ginv.Nrows(); i++)
        D(i) = 1.0 / coefs(i);

    return Ginv * D * Ginv;
}

Matrix ProfileMixtureProbability::ExpandCovariance(const std::string &modality, const ColumnVector &coefs)
{
    Matrix &G = m_fullG[modality];

    return G * coefs.AsDiagonal() * G;
}

Matrix ProfileMixtureProbability::ExpandSmallPrecision(const std::string &modality, const ColumnVector &coefs, int delta)
{
    Matrix &Gpinv = m_smallGpinv[modality][delta];
    DiagonalMatrix D(coefs.Nrows());
    D = 0.0;

    for (int i = 1; i <= coefs.Nrows(); i++)
        D(i) = 1.0 / coefs(i);

    return Gpinv.t() * D * Gpinv;
}

Matrix ProfileMixtureProbability::ExpandSmallCovariance(const std::string &modality, const ColumnVector &coefs, int delta)
{
    Matrix &G = m_smallG[modality][delta];

    return G * coefs.AsDiagonal() * G.t();
}

void ProfileMixtureProbability::NewMixingCoefs(const std::unordered_map<string, std::vector<double> > &coefs)
{
    BOOST_LOG_TRIVIAL(debug) << "Setting new mixing coefficients";

    m_currentMixingCoefsFull.clear();

    for (const auto &t : coefs)
    {
        m_currentMixingCoefsFull[t.first] = std::vector<double>(t.second.size() + 1);

        double s(0.0);
        for (std::size_t i = 0; i < t.second.size(); i++)
        {
            m_currentMixingCoefsFull[t.first][i] = t.second[i];
            s += t.second[i];
        }

        m_currentMixingCoefsFull[t.first][t.second.size()] = 1.0 - s;
    }

    m_precalcValid = false;
}

void ProfileMixtureProbability::NewComponentMeans(const std::unordered_map<std::string, std::vector<ColumnVector> > &means)
{
    BOOST_LOG_TRIVIAL(debug) << "Setting new means";

    m_currentComponentMeans = means;
    m_precalcValid = false;
}

void ProfileMixtureProbability::NewCovarianceCoefs(const std::unordered_map<std::string, std::vector<ColumnVector> > &coefs)
{
    BOOST_LOG_TRIVIAL(debug) << "Setting new covariance coefficients";

    m_currentCovarianceCoefs = coefs;
    m_precalcValid = false;
}

void ProfileMixtureProbability::Precalc()
{
    BOOST_LOG_TRIVIAL(debug) << "Precalc (creating MVN objects)";
    for (auto t : m_currentComponentMeans)
    {
        std::vector<ColumnVector> &means = m_currentComponentMeans[t.first];
        std::vector<ColumnVector> &covcoefs = m_currentCovarianceCoefs[t.first];

        std::vector<std::vector<Stats::MVN> > &mvec = m_precalcMVN[t.first];
        mvec = std::vector<std::vector<Stats::MVN> >(m_components);

        for (int r = 0; r < m_components; r++)
            for (int delta = 0; delta < m_steps; delta++)
                mvec[r].push_back(Stats::MVN(
                                      means[r].Rows(delta + 1, delta + m_dataLength),
                                      ExpandSmallPrecision(t.first, covcoefs[r], delta)));
    }

    BOOST_LOG_TRIVIAL(debug) << "Precalc (computing normalisation factors)";
    m_precalcZ.clear();
    for (std::size_t i = 0; i < m_data.cbegin()->second.size(); i++)
    {
        std::vector<double> logterms(m_steps);

        for (int delta = 0; delta < m_steps; delta++)
            logterms[delta] = Stats::lognormal(delta, m_deltaMean, m_deltaPrecision);

        for (const auto &t : m_data)
        {
            std::vector<std::vector<Stats::MVN> > &mvec = m_precalcMVN[t.first];
            std::vector<double> &mixcoefs = m_currentMixingCoefsFull[t.first];

            for (int delta = 0; delta < m_steps; delta++)
            {
                double lse_r = -std::numeric_limits<double>::infinity();

                for (int r = 0; r < m_components; r++)
                    lse_r = Stats::logsumexp(lse_r, mvec[r][delta].LogP(t.second[i]) + std::log(mixcoefs[r]));

                logterms[delta] += lse_r;
            }
        }

        m_precalcZ.push_back(Stats::logsumexp(logterms));
    }

    m_precalcValid = true;
}

double ProfileMixtureProbability::LogProbability()
{
    BOOST_LOG_TRIVIAL(debug) << "Calculating log probability";

    if (!m_precalcValid)
    {
        BOOST_LOG_TRIVIAL(debug) << "Pracalc is not up to date";
        Precalc();
        BOOST_LOG_TRIVIAL(debug) << "Precalc returned";
    }

    double value = 0.0;

    for (const auto &t : m_currentMixingCoefsFull)
    {
        ColumnVector theta(m_components - 1);
        for (int i = 0; i < m_components - 1; i++)
            theta(i + 1) = t.second[i];

        ColumnVector alpha(m_components);
        alpha << m_mixingAlpha.data();

        value += Stats::Dirichlet(alpha).LogP(theta);

        std::vector<ColumnVector> &means = m_currentComponentMeans[t.first];
        std::vector<ColumnVector> &covcoefs = m_currentCovarianceCoefs[t.first];

        for (int j = 0; j < m_components; j++)
            value += Stats::MVNWishart(m_profilePriorMeans[t.first][j], m_profileN0, m_profileAlpha,
                    m_profileBetas[t.first][j]).LogP(means[j], ExpandPrecision(t.first, covcoefs[j]));
    }

    for (const auto &z : m_precalcZ)
        value += z;

    double tmpzs = 0.0;
    for (const auto &z : m_precalcZ)
        tmpzs += z;

    return value;
}

void ProfileMixtureProbability::Derivatives(
        std::unordered_map<std::string, std::vector<double> > &dLPdMixing,
        std::unordered_map<std::string, std::vector<ColumnVector> > &dLPdMeans,
        std::unordered_map<std::string, std::vector<ColumnVector> > &dLPdCovarianceCoefs)
{
    BOOST_LOG_TRIVIAL(debug) << "Calculating derivatives";

    dLPdMixing.clear();
    dLPdMeans.clear();

    if (!m_precalcValid)
    {
        BOOST_LOG_TRIVIAL(debug) << "Pracalc is not up to date";
        Precalc();
        BOOST_LOG_TRIVIAL(debug) << "Precalc returned";
    }

    std::size_t subjects = m_data.cbegin()->second.size();

    std::unordered_map<std::string, std::vector<std::vector<std::vector<double> > > > datalogps;
    for (const auto &s : m_data)
    {
        datalogps[s.first].resize(m_components);

        for (std::size_t j = 0; j < m_components; j++)
        {
            datalogps[s.first][j].resize(m_steps);

            for (int delta = 0; delta < m_steps; delta++)
            {
                datalogps[s.first][j][delta].resize(subjects);

                for (std::size_t i = 0; i < subjects; i++)
                    datalogps[s.first][j][delta][i] = m_precalcMVN[s.first][j][delta].LogP(s.second[i]);
            }
        }
    }

    std::unordered_map<std::string, std::vector<Matrix> > dLPdCovariance;

    for (const auto &s : m_data)
    {
        dLPdMixing[s.first] = std::vector<double>(m_components - 1, 0.0);
        dLPdMeans[s.first] = std::vector<ColumnVector>(m_components);
        dLPdCovariance[s.first] = std::vector<Matrix>(m_components);

        BOOST_LOG_TRIVIAL(debug) << " -> Derivatives for " << s.first;

        for (int j = 0; j < m_components; j++)
        {
            BOOST_LOG_TRIVIAL(debug) << " --> Component " << j;

            if (j < m_components - 1)
                dLPdMixing[s.first][j] = (1.0 - m_mixingAlpha[m_components - 1]) / m_currentMixingCoefsFull[s.first][m_components - 1]
                        + (m_mixingAlpha[j] - 1.0) / m_currentMixingCoefsFull[s.first][j];

            ColumnVector deltamu = m_currentComponentMeans[s.first][j] - m_profilePriorMeans[s.first][j];

            dLPdMeans[s.first][j] = -m_profileN0 * (deltamu.t()
                    * ExpandPrecision(s.first, m_currentCovarianceCoefs[s.first][j])).AsColumn();

            Matrix fullprec = ExpandPrecision(s.first, m_currentCovarianceCoefs[s.first][j]);
            dLPdCovariance[s.first][j] = (m_refLength / 2.0 - m_profileAlpha) * fullprec
                    + fullprec * (m_profileBetas[s.first][j] + m_profileN0 / 2.0 * deltamu * deltamu.t()) * fullprec;

            std::vector<Matrix> precs(m_steps);
            for (int delta = 0; delta < m_steps; delta++)
                precs[delta] = ExpandSmallPrecision(s.first, m_currentCovarianceCoefs[s.first][j], delta);

            for (std::size_t i = 0; i < subjects; i++)
            {
                for (int delta = 0; delta < m_steps; delta++)
                {
                    double logsum = - m_precalcZ[i] + Stats::lognormal(delta, m_deltaMean, m_deltaPrecision);

                    for (const auto &t : m_data)
                    {
                        if (s.first != t.first)
                        {
                            double lse_r = -std::numeric_limits<double>::infinity();

                            for (int r = 0; r < m_components; r++)
                                lse_r = Stats::logsumexp(lse_r, std::log(m_currentMixingCoefsFull[t.first][r])
                                        + datalogps[t.first][r][delta][i]);

                            logsum += lse_r;
                        }
                    }

                    if (j < m_components - 1)
                        dLPdMixing[s.first][j] += std::exp(logsum + datalogps[s.first][j][delta][i])
                                - std::exp(logsum + datalogps[s.first][m_components - 1][delta][i]);

                    logsum += std::log(m_currentMixingCoefsFull[s.first][j]) + datalogps[s.first][j][delta][i];
                    
                    ColumnVector deltax = s.second[i] - m_currentComponentMeans[s.first][j].Rows(1 + delta, m_dataLength + delta);
                    
                    dLPdMeans[s.first][j].Rows(1 + delta, m_dataLength + delta) += std::exp(logsum)
                            * (deltax.t() * precs[delta]).AsColumn();
                    
                    dLPdCovariance[s.first][j].SymSubMatrix(1 + delta, m_dataLength + delta)
                            += std::exp(logsum) * 0.5 * (precs[delta] * deltax * deltax.t() * precs[delta] - precs[delta]);
                }
            }
        }
    }

    dLPdCovarianceCoefs.clear();

    // NB: This also deals with symmetry of covariance matrix (i.e. dependent elements)
    for (const auto &s : dLPdCovariance)
    {
        dLPdCovarianceCoefs[s.first] = std::vector<ColumnVector>(m_components);
        for (int j = 0; j < m_components; j++)
        {
            Matrix gpg = m_fullG[s.first] * s.second[j] * m_fullG[s.first];
            dLPdCovarianceCoefs[s.first][j] = ColumnVector(gpg.Nrows());

            for (int i = 1; i <= gpg.Nrows(); i++)
                dLPdCovarianceCoefs[s.first][j](i) = gpg(i, i);
        }
    }
}

std::vector<std::string> ProfileMixtureProbability::SortModalities() const
{
    std::vector<std::string> names;
    for (const auto &t : m_data)
        names.push_back(t.first);

    std::sort(names.begin(), names.end());

    return names;
}

void ProfileMixtureProbability::UnpackX(const std::vector<double> &x,
                                        std::unordered_map<std::string, std::vector<double> > &mixing,
                                        std::unordered_map<std::string, std::vector<ColumnVector> > &means,
                                        std::unordered_map<std::string, std::vector<ColumnVector> > &covcoefs) const
{
    BOOST_LOG_TRIVIAL(debug) << "Unpacking point";

    std::size_t modalities = m_data.size();

    if (x.size() !=  modalities * (m_components - 1 + 2 * m_components * m_refLength))
        throw ProfileModel::ModelException("Packed vector has incorrect length");

    std::vector<std::string> names = SortModalities();

    mixing = std::unordered_map<std::string, std::vector<double> >();
    means = std::unordered_map<std::string, std::vector<ColumnVector> >();
    covcoefs = std::unordered_map<std::string, std::vector<ColumnVector> >();

    for (std::size_t k = 0; k < modalities; k++)
    {
        mixing[names[k]] = std::vector<double>(m_components - 1);
        for (int j = 0; j < m_components - 1; j++)
            mixing[names[k]][j] = x[k * (m_components - 1) + j];

        means[names[k]] = std::vector<ColumnVector>(m_components);
        covcoefs[names[k]] = std::vector<ColumnVector>(m_components);
        for (int j = 0; j < m_components; j++)
        {
            means[names[k]][j] = ColumnVector(m_refLength);
            covcoefs[names[k]][j] = ColumnVector(m_refLength);
            for (int i = 0; i < m_refLength; i++)
            {
                means[names[k]][j](i + 1) = x[modalities * (m_components - 1) + (k * m_components + j) * m_refLength + i];
                covcoefs[names[k]][j](i + 1) = x[modalities * (m_components - 1) + ((k + modalities) * m_components + j) * m_refLength + i];
            }
        }
    }
}

std::vector<double> ProfileMixtureProbability::PackX(const std::unordered_map<string, std::vector<double> > &mixing,
                                                     const std::unordered_map<string, std::vector<ColumnVector> > &means,
                                                     const std::unordered_map<string, std::vector<ColumnVector> > &covcoefs) const
{
    BOOST_LOG_TRIVIAL(debug) << "Packing point";

    std::size_t modalities = m_data.size();

    std::vector<double> x(modalities * (m_components - 1 + 2 * m_components * m_refLength));

    std::vector<std::string> names = SortModalities();
    for (std::size_t k = 0; k < modalities; k++)
    {
        for (int j = 0; j < m_components - 1; j++)
            x[k * (m_components - 1) + j] = mixing.find(names[k])->second[j];

        for (int j = 0; j < m_components; j++)
        {
            for (int i = 0; i < m_refLength; i++)
            {
                x[modalities * (m_components - 1) + (k * m_components + j) * m_refLength + i] = means.find(names[k])->second[j](i + 1);
                x[modalities * (m_components - 1) + ((modalities + k) * m_components + j) * m_refLength + i] = covcoefs.find(names[k])->second[j](i + 1);
            }
        }
    }

    return x;
}

void ProfileMixtureProbability::UnpackPoint(const std::vector<double> &x)
{
    std::unordered_map<std::string, std::vector<double> > mixing;
    std::unordered_map<std::string, std::vector<ColumnVector> > means;
    std::unordered_map<std::string, std::vector<ColumnVector> > covcoefs;

    UnpackX(x, mixing, means, covcoefs);

    NewMixingCoefs(mixing);
    NewComponentMeans(means);
    NewCovarianceCoefs(covcoefs);
}

std::vector<double> ProfileMixtureProbability::PackDerivatives()
{
    std::unordered_map<std::string, std::vector<double> > mixingderiv;
    std::unordered_map<std::string, std::vector<ColumnVector> > meansderiv;
    std::unordered_map<std::string, std::vector<ColumnVector> > covcoefsderiv;

    Derivatives(mixingderiv, meansderiv, covcoefsderiv);

    return PackX(mixingderiv, meansderiv, covcoefsderiv);
}

double ProfileMixtureProbability::Objective(const std::vector<double> &x, std::vector<double> &grad, void *f_data)
{
    BOOST_LOG_TRIVIAL(debug) << "Entering cost function";

    double result;

    try
    {
        auto *obj = static_cast<ProfileMixtureProbability *>(f_data);

        obj->UnpackPoint(x);
        result = obj->LogProbability();

        BOOST_LOG_TRIVIAL(info) << "Log probability: " << result;

        if (!grad.empty())
        {
            BOOST_LOG_TRIVIAL(debug) << " + derivatives";
            grad = obj->PackDerivatives();

            // Derivative check
            if (false)
            {
                RowVector g(grad.size());
                g << grad.data();

                std::vector<double> fdgrad(x.size());

                double eps = 1e-3;

                for (std::size_t i = 0; i < x.size(); i++)
                {
                    std::vector<double> x2 = x;
                    x2[i] += eps;
                    obj->UnpackPoint(x2);
                    double prob2 = obj->LogProbability();
                    fdgrad[i] = (prob2 - result) / eps;
                    BOOST_LOG_TRIVIAL(debug) << "Gradient check (dim " << i << ") - Analytical: " << grad[i]
                                                << ", finite difference: " << fdgrad[i];
                }

                BOOST_LOG_TRIVIAL(debug) << "Differences between finite differences and analytical derivatives:";
                for (std::size_t i = 0; i < fdgrad.size(); i++)
                    BOOST_LOG_TRIVIAL(debug) << "(fdgrad[" << i << "] - grad[" << i << "]) / grad[" << i << "] = "
                                             << (fdgrad[i] - grad[i]) / grad[i];
            }
        }
    }
    catch (const Exception &e)
    {
        BOOST_LOG_TRIVIAL(error) << "Exception in cost function (rethrowing): " << e.what();
        throw;
    }
    catch (const std::exception &e)
    {
        BOOST_LOG_TRIVIAL(error) << "Exception in cost function (rethrowing): " << e.what();
        throw;
    }

//    if (std::isnan(result))
//    {
//        BOOST_LOG_TRIVIAL(warning) << "Log probability is NaN - halting optimisation";
//        throw nlopt::forced_stop();
//    }

    return result;
}

std::vector<double> ProfileMixtureProbability::LowerBounds() const
{
    std::size_t modalities = m_data.size();

    std::vector<double> lb(modalities * (m_components - 1 + 2 * m_components * m_refLength), -HUGE_VAL);

    for (std::size_t i = 0; i < modalities * (m_components - 1); i++)
        lb[i] = 1e-10;

    for (std::size_t i = modalities * (m_components - 1 + m_components * m_refLength); i < lb.size(); i++)
        lb[i] = 1e-8;

    return lb;
}

std::vector<double> ProfileMixtureProbability::UpperBounds() const
{
    std::size_t modalities = m_data.size();

    std::vector<double> ub(modalities * (m_components - 1 + 2 * m_components * m_refLength), HUGE_VAL);

    for (std::size_t i = 0; i < modalities * (m_components - 1); i++)
        ub[i] = 1.0 - 1e-10;

    return ub;
}

void ProfileMixtureProbability::MixingInequality(unsigned int m, double *result, unsigned n,
                                                   const double *x, double *grad, void *f_data)
{
    BOOST_LOG_TRIVIAL(debug) << "Entering inequality constraint function for mixing coefficients";

    const auto *obj = static_cast<ProfileMixtureProbability *>(f_data);

    // n is checked by UnpackX() below
    if (m != obj->m_data.size())
        throw ProfileMixtures::ModelException("Constraint gradient output array has incorrect size");

    std::unordered_map<std::string, std::vector<double> > mixing;
    std::unordered_map<std::string, std::vector<ColumnVector> > means;
    std::unordered_map<std::string, std::vector<ColumnVector> > covcoefs;

    std::vector<double> xvec(x, x + n);
    obj->UnpackX(xvec, mixing, means, covcoefs);

    std::vector<std::string> sortedmodalities = obj->SortModalities();

    std::unordered_map<std::string, std::vector<double> > mixingzeros;
    std::unordered_map<std::string, std::vector<ColumnVector> > meanszeros;
    std::unordered_map<std::string, std::vector<ColumnVector> > covcoefszeros;

    ColumnVector zerocv(means.cbegin()->second.at(0).Nrows());
    zerocv = 0.0;

    for (const auto &modality : sortedmodalities)
    {
        mixingzeros[modality] = std::vector<double>(mixing[modality].size(), 0.0);

        for (int c = 0; c < means[modality].size(); c++)
        {
            meanszeros[modality].push_back(zerocv);
            covcoefszeros[modality].push_back(zerocv);
        }
    }

    for (const auto &modality : sortedmodalities)
    {
        const std::vector<double> &coefs = mixing[modality];
        *result++ = std::accumulate(coefs.cbegin(), coefs.cend(), 0.0) - 1.0;

        if (grad != nullptr)
        {
            std::unordered_map<std::string, std::vector<double> > mixinggrad = mixingzeros;
            mixinggrad[modality] = std::vector<double>(mixing[modality].size(), 1.0);

            std::vector<double> constraintgrad = obj->PackX(mixinggrad, meanszeros, covcoefszeros);

            for (const auto v : constraintgrad)
                *grad++ = v;
        }
    }
}

std::vector<double> ProfileMixtureProbability::GetInitialGuess() const
{
    std::unordered_map<std::string, std::vector<double> > mixing;

    for (auto t : m_profilePriorMeans)
        for (int j = 0; j < m_components - 1; j++)
            mixing[t.first].push_back(1.0 / m_components);

    std::vector<double> x0 = PackX(mixing, m_profilePriorMeans, m_profilePriorCovarianceCoefs);

//    std::random_device rd;
//    std::mt19937 mt(rd());
//    std::uniform_real_distribution<> dist(0.98, 1.02);
//    for (auto &it : x0)
//        it *= dist(mt);

    return x0;
}
