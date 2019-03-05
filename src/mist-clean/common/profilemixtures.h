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

#ifndef PROFILEMIXTURES_H
#define PROFILEMIXTURES_H

#include "profilemodel.h"
#include "profilepriors.h"
#include "stats.h"
#include "plotting.h"
#include "newmat.h"
#include <vector>

class ProfileMixtures;
BOOST_CLASS_EXPORT_KEY(ProfileMixtures)

class ProfileMixtureProbability;

namespace boost
{
namespace serialization
{
    template<class Archive>
    inline void save_construct_data(Archive &ar, const ProfileMixtures *m, const unsigned int version);
}
}

class ProfileMixtures : public ProfileModel
{
public:
    ProfileMixtures(const std::vector<string> &modalitynames, int reflength, int datalength, int components);

    // NB: Using default copy constructor / operator=

    virtual boost::shared_ptr<ProfileModel> Clone() const;

    virtual void Train(const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &trainingdata);

    virtual std::vector<double> GetDeltaLikelihoods(const std::unordered_map<std::string, NEWMAT::ColumnVector> &data,
                                                    bool usedeltaprior) const;

    // NB: SetComponentPriorMean() applies the smoothing matrix to the mean profile; the result of GetComponentPriorMean()
    // will not be identical to what is passed here
    void SetComponentPriorMean(const std::string &modality, int component, boost::shared_ptr<const ProfilePrior> mean, double spacing);
    void SetComponentPriorCovarianceCoefs(const std::string &modality, int component, boost::shared_ptr<const ProfilePrior> coefs, double spacing);
    void SetSmoothness(const std::string &modality, double smoothness);
    void SetN0(int n0);
    void SetAlpha(double alpha);
    void SetDeltaStdev(double stdev);
    void SetMixingAlpha(const std::vector<double> &mixingalpha);

    void SetFTolerance(double ftol);
    void SetFToleranceIgnoreCount(double count);
    void SetMaxEvaluations(int maxeval);
    void SetMaxRetries(int maxretries);

    void SetUseMCMC(bool usemcmc);

    int GetNumberOfComponents() const;

    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > GetComponentPriorMeans() const;
    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > GetComponentPriorCovarianceCoefs() const;

    std::unordered_map<std::string, std::vector<double> > GetFullMixingCoefs() const;
    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > GetComponentMeans() const;
    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > GetComponentCovarianceCoefs() const;

    std::unordered_map<std::string, std::vector<std::vector<double> > > GetFullMixingCoefsSamples() const;
    std::unordered_map<std::string, std::vector<std::vector<NEWMAT::ColumnVector> > > GetComponentMeansSamples() const;
    std::unordered_map<std::string, std::vector<std::vector<NEWMAT::ColumnVector> > > GetComponentCovarianceCoefsSamples() const;

    std::unordered_map<std::string, double> GetSmoothness() const;

    static std::vector<double> FullMixingCoefs(const std::vector<double> &coefs);

    virtual void WritePlots(const Plotting::plotfunc &pfunc, int vertex,
                            std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > *profiles = nullptr) const;

private:
    int m_components;
    double m_deltaMean;
    double m_deltaPrecision;
    std::vector<double> m_mixingAlpha;
    double m_alpha;
    int m_n0;

    double m_ftol;
    int m_ftolIgnore;
    int m_maxEval;
    int m_maxRetries;

    bool m_useMCMC;
    int m_mcmcIterations;
    int m_mcmcBurnIn;
    int m_mcmcThin;
    double m_mcmcScaleMixing;
    double m_mcmcScaleMeans;
    double m_mcmcScaleCovarianceCoefs;

    struct ModalityParameters
    {
        // Smoothness is sigma in points
        double smoothness;

        // These are set in the constructor
        std::vector<NEWMAT::ColumnVector> componentPriorMeans;
        std::vector<NEWMAT::ColumnVector> componentPriorCovarianceCoefs;

        // Variables below are set during training
        std::vector<double> mixingCoefs;
        std::vector<NEWMAT::ColumnVector> componentMeans;
        std::vector<NEWMAT::ColumnVector> componentCovarianceCoefs;

        std::vector<std::vector<double> > mixingCoefsSamples;
        std::vector<std::vector<NEWMAT::ColumnVector> > componentMeansSamples;
        std::vector<std::vector<NEWMAT::ColumnVector> > componentCovarianceCoefsSamples;

        NEWMAT::Matrix G;
        // Need to store these to speed up GetDeltaLikelihoods()
        std::vector<NEWMAT::Matrix> verySmallGpinv;

    private:
        friend class boost::serialization::access;
        template<class Archive>
        void serialize(Archive &ar, const unsigned int version);
    };

    std::unordered_map<std::string, ModalityParameters> m_modalityParameters;

    std::vector<double> Optimise(ProfileMixtureProbability &pmp);
    std::vector<std::vector<double> > GenerateSamples(ProfileMixtureProbability &pmp, std::vector<double> current);

    friend class boost::serialization::access;
    template<class Archive>
    friend void boost::serialization::save_construct_data(Archive &ar, const ProfileMixtures *m, const unsigned int version);
    template<class Archive>
    void save(Archive &ar, const unsigned int version) const;
    template<class Archive>
    void load(Archive &ar, const unsigned int version);

    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

template<class Archive>
void ProfileMixtures::ModalityParameters::serialize(Archive &ar, const unsigned int version)
{
    BOOST_LOG_TRIVIAL(debug) << "(De)serialising ProfileMixtures::ModalityParameters";

    ar & smoothness;
    ar & componentPriorMeans;
    ar & componentPriorCovarianceCoefs;
    ar & mixingCoefs;
    ar & componentMeans;
    ar & componentCovarianceCoefs;
    ar & mixingCoefsSamples;
    ar & componentMeansSamples;
    ar & componentCovarianceCoefsSamples;
}

template<class Archive>
void ProfileMixtures::save(Archive &ar, const unsigned int version) const
{
    BOOST_LOG_TRIVIAL(debug) << "Serialising ProfileMixtures";

    ar << boost::serialization::base_object<ProfileModel>(*this);

    ar << m_mixingAlpha;
    ar << m_alpha;
    ar << m_n0;
    ar << m_deltaPrecision;

    ar << m_ftol;
    ar << m_ftolIgnore;
    ar << m_maxEval;
    ar << m_maxRetries;

    ar << m_useMCMC;
    ar << m_mcmcIterations;
    ar << m_mcmcBurnIn;
    ar << m_mcmcThin;
    ar << m_mcmcScaleMixing;
    ar << m_mcmcScaleMeans;
    ar << m_mcmcScaleCovarianceCoefs;

    ar << m_modalityParameters;
}

template<class Archive>
void ProfileMixtures::load(Archive &ar, const unsigned int version)
{
    BOOST_LOG_TRIVIAL(debug) << "Deserialising ProfileMixtures";

    ar >> boost::serialization::base_object<ProfileModel>(*this);

    ar >> m_mixingAlpha;
    ar >> m_alpha;
    ar >> m_n0;
    ar >> m_deltaPrecision;

    ar >> m_ftol;
    ar >> m_ftolIgnore;
    ar >> m_maxEval;
    ar >> m_maxRetries;

    ar >> m_useMCMC;
    ar >> m_mcmcIterations;
    ar >> m_mcmcBurnIn;
    ar >> m_mcmcThin;
    ar >> m_mcmcScaleMixing;
    ar >> m_mcmcScaleMeans;
    ar >> m_mcmcScaleCovarianceCoefs;

    ar >> m_modalityParameters;

    for (const auto &mp : m_modalityParameters)
        SetSmoothness(mp.first, mp.second.smoothness);
}

namespace boost
{
namespace serialization
{
    template<class Archive>
    inline void save_construct_data(Archive &ar, const ProfileMixtures *m, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for ProfileMixtures";

        int modalities = m->m_modalityParameters.size();
        ar << modalities;
        for (const auto &t : m->m_modalityParameters)
            ar << t.first;

        ar << m->m_refLength;
        ar << m->m_dataLength;
        ar << m->m_components;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, ProfileMixtures *m, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for ProfileMixtures";

        std::vector<std::string> modalitynames;
        int modalities;

        ar >> modalities;
        while (modalities--)
        {
            std::string modality;
            ar >> modality;
            modalitynames.push_back(modality);
        }

        int reflength;
        int datalength;
        int components;

        ar >> reflength;
        ar >> datalength;
        ar >> components;

        BOOST_LOG_TRIVIAL(debug) << "Calling placement new for ProfileMixtures";
        ::new(m)ProfileMixtures(modalitynames, reflength, datalength, components);
    }
}
}

// -----------------------------------------------------------------------------------------------------------------------------

class ProfileMixtureProbability
{
public:
    ProfileMixtureProbability(const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &data,
                              const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &profilepriormeans,
                              const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &profilepriorcovcoefs,
                              int profilen0, int profilealpha, double deltamean, double deltaprecision,
                              const std::vector<double> &mixingalpha,
                              const std::unordered_map<std::string, double> &stdevsmooth);

    static NEWMAT::Matrix GetSmoothingMatrix(double sigmasq, std::size_t size);

    int GetProfileN0() const;
    double GetProfileAlpha() const;
    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > GetProfilePriorMeans() const;
    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > GetProfilePriorCovarianceCoefs() const;

    // Final component should not be specified (i.e. this is theta)
    void NewMixingCoefs(const std::unordered_map<std::string, std::vector<double> > &coefs);
    void NewComponentMeans(const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &means);
    void NewCovarianceCoefs(const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &coefs);

    void UnpackX(const std::vector<double> &x,
                 std::unordered_map<std::string, std::vector<double> > &mixing,
                 std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &means,
                 std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &covcoefs) const;
    std::vector<double> PackX(const std::unordered_map<std::string, std::vector<double> > &mixing,
                              const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &means,
                              const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &covcoefs) const;
    void UnpackPoint(const std::vector<double> &x);
    std::vector<double> PackDerivatives();

    double LogProbability();
    void Derivatives(
            std::unordered_map<std::string, std::vector<double> > &dLPdMixing,
            std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &dLPdMeans,
            std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &dLPdCovarianceCoefs);

    NEWMAT::Matrix ExpandPrecision(const std::string &modalilty, const NEWMAT::ColumnVector &coefs);
    NEWMAT::Matrix ExpandCovariance(const std::string &modalilty, const NEWMAT::ColumnVector &coefs);
    NEWMAT::Matrix ExpandSmallPrecision(const std::string &modalilty, const NEWMAT::ColumnVector &coefs, int delta);
    NEWMAT::Matrix ExpandSmallCovariance(const std::string &modalilty, const NEWMAT::ColumnVector &coefs, int delta);

    static double Objective(const std::vector<double> &x, std::vector<double> &grad, void *f_data);
    std::vector<double> LowerBounds() const;
    std::vector<double> UpperBounds() const;
    static void MixingInequality(unsigned int m, double *result, unsigned n, const double *x, double *grad, void *f_data);
    std::vector<double> GetInitialGuess() const;

private:
    int m_dataLength;
    int m_refLength;
    int m_components;
    int m_subjects;
    int m_steps;
    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > m_data;

    double m_deltaMean;
    double m_deltaPrecision;
    std::vector<double> m_mixingAlpha;

    int m_profileN0;
    double m_profileAlpha;

    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > m_profilePriorMeans;
    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > m_profilePriorCovarianceCoefs;
    std::unordered_map<std::string, std::vector<NEWMAT::Matrix> > m_profileBetas;

    std::unordered_map<std::string, std::vector<double> > m_currentMixingCoefsFull;
    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > m_currentComponentMeans;
    std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > m_currentCovarianceCoefs;

    std::unordered_map<std::string, NEWMAT::Matrix> m_fullG;
    std::unordered_map<std::string, NEWMAT::Matrix> m_fullGinv;
    std::unordered_map<std::string, std::vector<NEWMAT::Matrix> > m_smallG;
    std::unordered_map<std::string, std::vector<NEWMAT::Matrix> > m_smallGpinv;

    bool m_precalcValid = false;
    std::vector<double> m_precalcZ;
    std::unordered_map<std::string, std::vector<std::vector<Stats::MVN> > > m_precalcMVN;

    void Precalc();
    std::vector<std::string> SortModalities() const;
};

#endif // PROFILEMIXTURES_H
