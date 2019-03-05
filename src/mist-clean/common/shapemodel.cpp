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

#include "shapemodel.h"
#include "stats.h"
#include "newimageio.h"
#include <boost/log/trivial.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <algorithm>
#include <memory>

#define WANT_STREAM
#include "newmatio.h"
#include "boost/make_shared.hpp"

using namespace NEWMAT;
using namespace NEWIMAGE;

ShapeModel::ShapeModel(const Shape &shape, const std::vector<std::string> &modalitynames, int profilepoints, double profilespacing,
                       bool usenormalisationmasks) :
    m_shape(new Shape(shape)),
    m_modalityNames(modalitynames),
    m_profilePoints(profilepoints),
    m_profileSpacing(profilespacing),
    m_useNormalisationMasks(usenormalisationmasks)
{
    for (const auto &mn : modalitynames)
    {
        SetNormalise(mn, NormalisationMode::None);
        SetFilter(mn, boost::dynamic_pointer_cast<ProfileFilter>(boost::make_shared<IdentityFilter>()));
    }
}

ShapeModel::~ShapeModel()
{
}

std::vector<std::string> ShapeModel::GetModalityNames() const
{
    return m_modalityNames;
}

int ShapeModel::GetProfilePoints() const
{
    return m_profilePoints;
}

double ShapeModel::GetProfileSpacing() const
{
    return m_profileSpacing;
}

boost::shared_ptr<const Shape> ShapeModel::GetShape() const
{
    return m_shape;
}

std::unordered_map<std::string, double> ShapeModel::GetMeanIntensities() const
{
    return m_meanIntensities;
}

boost::shared_ptr<ProfileModel> ShapeModel::GetVertexModel(int vertex) const
{
    if (vertex < 0 || vertex >= m_vertexModels.size())
        throw ModelException("Cannot return vertex model; vertex number is out of range");

    return m_vertexModels[vertex];
}

void ShapeModel::SetMaxIterations(int maxiter)
{
    m_maxiter = maxiter;
}

bool ShapeModel::GetUseNormalisationMasks() const
{
    return m_useNormalisationMasks;
}

void ShapeModel::SetNormalise(string modality, NormalisationMode mode)
{
    m_normalise[modality] = mode;
}

ShapeModel::NormalisationMode ShapeModel::GetNormalise(std::string modality) const
{
    CheckModality(modality);

    return m_normalise.find(modality)->second;
}

ColumnVector ShapeModel::NormaliseProfile(string modality, const ColumnVector &profile, double mean) const
{
    if (m_meanIntensities.find(modality) == m_meanIntensities.end())
        throw ModelException("NormaliseProfile() called before setting m_meanIntensities for the specified modality");

    ColumnVector normalised;

    switch (m_normalise.find(modality)->second)
    {
        default:
        case NormalisationMode::None:
            BOOST_LOG_TRIVIAL(trace) << "Normalisation mode for " << modality << ": None";
            normalised = profile;
            break;

        case NormalisationMode::Additive:
            BOOST_LOG_TRIVIAL(trace) << "Normalisation mode for " << modality << ": Additive";
            normalised = profile - mean + m_meanIntensities.find(modality)->second;
            break;

        case NormalisationMode::Multiplicative:
            BOOST_LOG_TRIVIAL(trace) << "Normalisation mode for " << modality << ": Multiplicative";
            normalised = profile / mean * m_meanIntensities.find(modality)->second;
            break;
    }

    return normalised;
}

void ShapeModel::SetFilter(string modality, boost::shared_ptr<ProfileFilter> filter)
{
    m_filters[modality] = filter;
}

boost::shared_ptr<ProfileFilter> ShapeModel::GetFilter(string modality) const
{
    return m_filters.find(modality)->second;
}

void ShapeModel::CheckModality(const string &modalityname) const
{
    if (std::find(m_modalityNames.cbegin(), m_modalityNames.cend(), modalityname) == m_modalityNames.cend())
        throw ModelException("Specified modality does not exist in shape model");
}

void ShapeModel::CheckFilenames(const std::unordered_map<std::string, std::vector<std::string> > &trainingdata,
                                const std::vector<boost::shared_ptr<const Transformation> > &transformations,
                                const std::vector<std::string> &normalisationmasks) const
{
    if (m_useNormalisationMasks && normalisationmasks.empty())
        throw ModelException("The model requires exclusion masks for normalisation, but none were specified");

    if (!m_useNormalisationMasks && !normalisationmasks.empty())
        throw ModelException("Exclusion masks were specified for normalisation, but the model does not use these");

    ModelException maperror("Training data should be provided for all modalities that were included when constructing the model");

    if (trainingdata.size() != m_modalityNames.size())
        throw maperror;

    for (const auto &t : trainingdata)
    {
        if (std::find(m_modalityNames.cbegin(), m_modalityNames.cend(), t.first) == m_modalityNames.cend())
            throw maperror;

        if (t.second.size() != trainingdata.cbegin()->second.size())
            throw ModelException("The number of images provided is not the same for all modalities");
    }

    if (transformations.size() != trainingdata.cbegin()->second.size())
        throw ModelException("Incorrect number of transformations");

    if (!normalisationmasks.empty() && normalisationmasks.size() != trainingdata.cbegin()->second.size())
        throw ModelException("Incorrect number of exclusion masks for normalisation (need either none or one for each subject)");
}

void ShapeModel::CreateVertexModels(boost::shared_ptr<const ProfileModel> example)
{
    m_vertexModels.clear();

    for (int i = 0; i < m_shape->GetNumberOfVertices(); i++)
        m_vertexModels.push_back(example->Clone());
}

void ShapeModel::Train(const std::unordered_map<std::string, std::vector<std::string> > &trainingdata,
                       const std::vector<boost::shared_ptr<const Transformation> > &transformations,
                       const std::vector<std::string> &normalisationmasks)
{
    BOOST_LOG_TRIVIAL(debug) << "Training shape model (base class)";

    if (m_vertexModels.size() != m_shape->GetNumberOfVertices())
        throw ModelException("CreateVertexModels() must be called before Train()");

    CheckFilenames(trainingdata, transformations, normalisationmasks);

    std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > allprofiles
            = SampleAllVertices(trainingdata, transformations, true, normalisationmasks);

    for (int i = 0; i < m_shape->GetNumberOfVertices(); i++)
    {
        BOOST_LOG_TRIVIAL(info) << "Training vertex (vertex " << i << ")";
        TrainVertexImplementation(i, allprofiles);
    }

    TrainShapeImplementation(allprofiles);
}

std::vector<boost::shared_ptr<ProfileModel> > ShapeModel::TrainVertices(const std::unordered_map<string, std::vector<string> > &trainingdata,
                                                                        const std::vector<boost::shared_ptr<const Transformation> > &transformations,
                                                                        const std::vector<string> &normalisationmasks,
                                                                        std::vector<int> verts,
                                                                        const Plotting::plotfunc &pfunc)
{
    if (m_vertexModels.size() != m_shape->GetNumberOfVertices())
        throw ModelException("CreateVertexModels() must be called before TrainVertices()");

    CheckFilenames(trainingdata, transformations, normalisationmasks);

    std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > profiles =
            SampleAllVertices(trainingdata, transformations, true, normalisationmasks);

    std::vector<boost::shared_ptr<ProfileModel> > pms;
    for (auto vert : verts)
    {
        BOOST_LOG_TRIVIAL(info) << "Training vertex (vertex " << vert << ")";
        boost::shared_ptr<ProfileModel> pm = TrainVertexImplementation(vert, profiles);

        if (pfunc)
        {
            std::unordered_map<std::string, std::vector<ColumnVector> > vertprofiles;
            for (const auto &t : profiles)
                for (const auto &p : t.second[vert])
                    vertprofiles[t.first].push_back(p);

            pm->WritePlots(pfunc, vert, &vertprofiles);
        }

        pms.push_back(pm);
    }

    return pms;
}

void ShapeModel::LoadVertexModelsAndTrain(const std::vector<string> &archivefiles,
                                          const std::unordered_map<string, std::vector<string> > &trainingdata,
                                          const std::vector<boost::shared_ptr<const Transformation> > &transformations,
                                          const std::vector<std::string> &normalisationmasks)
{
    CheckFilenames(trainingdata, transformations, normalisationmasks);

    m_vertexModels.clear();

    for (const auto &af : archivefiles)
    {
        BOOST_LOG_TRIVIAL(info) << "Loading profile model from " << af;
        std::ifstream istr(af);
        boost::archive::text_iarchive in(istr);
        boost::shared_ptr<ProfileModel> pm;
        in >> pm;
        m_vertexModels.push_back(pm);
    }

    std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > allprofiles = SampleAllVertices(trainingdata, transformations, true, normalisationmasks);

    TrainShapeImplementation(allprofiles);
}

ColumnVector ShapeModel::Fit(const std::unordered_map<string, string> &data,
                             const boost::shared_ptr<const Transformation> &transformation,
                             const std::string &normalisationmask,
                             ColumnVector &cilower, ColumnVector &ciupper, double cilevel)
{
    BOOST_LOG_TRIVIAL(info) << "Fitting shape to new dataset";

    if (m_vertexModels.size() != m_shape->GetNumberOfVertices())
        throw ModelException("Cannot fit shape model with incorrect number of vertex models");

    if (m_useNormalisationMasks && normalisationmask.empty())
        throw ModelException("The model requires an exclusion mask for normalisation");

    if (!m_useNormalisationMasks && !normalisationmask.empty())
        throw ModelException("An exclusion mask was specified for normalisation, but the model cannot use one");

    std::unordered_map<std::string, std::vector<ColumnVector> > profiles;

    for (const auto &t : data)
    {
        if (std::find(m_modalityNames.cbegin(), m_modalityNames.cend(), t.first) == m_modalityNames.cend())
            throw ModelException(std::string("Modality ") + t.first + " does not exist in model");

        BOOST_LOG_TRIVIAL(info) << "Loading modality " << t.first;

        volume<double> vol;
        read_volume(vol, t.second);

        std::unique_ptr<volume<double> > mask = nullptr;
        if (!normalisationmask.empty())
        {
            mask.reset(new volume<double>);
            read_volume(*mask, normalisationmask);
        }

        double xmin, xmax, ymin, ymax, zmin, zmax;
        GetShapeExtent(transformation, xmin, xmax, ymin, ymax, zmin, zmax);
        double mean = GetMeanIntensity(vol, mask.get(), xmin, xmax, ymin, ymax, zmin, zmax);

        for (int vert = 0; vert < m_shape->GetNumberOfVertices(); vert++)
        {
            BOOST_LOG_TRIVIAL(debug) << "Vertex " << vert;

            ColumnVector profile = SampleProfile(vert, m_profilePoints, m_profileSpacing, vol, transformation);

            profile = NormaliseProfile(t.first, profile, mean);
            profile = m_filters[t.first]->Filter(profile);

            profiles[t.first].push_back(profile);
        }
    }

    return FitImplementation(profiles, cilower, ciupper, cilevel);
}

std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > ShapeModel::SampleAllVertices(
        const std::unordered_map<string, std::vector<string> > &vols,
        const std::vector<boost::shared_ptr<const Transformation> > &transformations,
        bool computemeans,
        const std::vector<std::string> &normalisationmasks)
{
    std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > allprofiles;

    for (const auto &v : vols)
    {
        std::vector<std::vector<ColumnVector> > profiles(m_shape->GetNumberOfVertices());
        std::vector<double> means;

        for (std::size_t i = 0; i < v.second.size(); i++)
        {
            volume<double> vol;
            read_volume(vol, v.second[i]);

            std::unique_ptr<volume<double> > mask = nullptr;
            if (!normalisationmasks.empty())
            {
                mask.reset(new volume<double>);
                read_volume(*mask, normalisationmasks[i]);
            }

            double xmin, xmax, ymin, ymax, zmin, zmax;
            GetShapeExtent(transformations[i], xmin, xmax, ymin, ymax, zmin, zmax);
            means.push_back(GetMeanIntensity(vol, mask.get(), xmin, xmax, ymin, ymax, zmin, zmax));

            BOOST_LOG_TRIVIAL(info) << "Modality " << v.first << ", subject " << i << ": Mean = " << means[i]
                << " (bounds: " << xmin << ", " << xmax << ", " << ymin << ", " << ymax << ", " << zmin << ", " << zmax << ")";

            for (int vert = 0; vert < m_shape->GetNumberOfVertices(); vert++)
                profiles[vert].push_back(SampleProfile(vert, m_profilePoints, m_profileSpacing, vol, transformations[i]));

            transformations[i]->ReleaseMemory();
        }

        if (computemeans)
        {
            m_meanIntensities[v.first] = std::accumulate(means.cbegin(), means.cend(), 0.0) / means.size();
            BOOST_LOG_TRIVIAL(info) << "Modality " << v.first << ": Mean = " << m_meanIntensities[v.first];
        }
        else if (m_meanIntensities.find(v.first) == m_meanIntensities.end())
            throw ModelException("Cannot normalise as mean intensities have not been computed");

        for (auto &p : profiles)
            for (std::size_t i = 0; i < p.size(); i++)
                p[i] = NormaliseProfile(v.first, p[i], means[i]);

        for (auto &p : profiles)
            for (std::size_t i = 0; i < p.size(); i++)
                p[i] = m_filters[v.first]->Filter(p[i]);

        allprofiles[v.first] = profiles;
    }

    return allprofiles;
}

double ShapeModel::GetMeanIntensity(volume<double> vol, volume<double> *exclusionmask,
                                    double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) const
{
    ColumnVector voxelsize(3);
    voxelsize << vol.xdim() << vol.ydim() << vol.zdim();

    int x1 = static_cast<int>(xmin / voxelsize(1) + 0.5);
    int x2 = static_cast<int>(xmax / voxelsize(1) + 0.5);
    int y1 = static_cast<int>(ymin / voxelsize(2) + 0.5);
    int y2 = static_cast<int>(ymax / voxelsize(2) + 0.5);
    int z1 = static_cast<int>(zmin / voxelsize(3) + 0.5);
    int z2 = static_cast<int>(zmax / voxelsize(3) + 0.5);

    double sum = 0.0;
    long int count = 0;

    for (int z = z1 ; z <= z2; z++)
    {
        for (int y = y1; y <= y2; y++)
        {
            for (int x = x1; x <= x2; x++)
            {
                if (exclusionmask == nullptr || (*exclusionmask)(x, y, z) <= 0.0)
                {
                    sum += vol(x, y, z);
                    count++;
                }
            }
        }
    }

    return sum / count;
}

void ShapeModel::GetShapeExtent(const boost::shared_ptr<const Transformation> &transformation,
                                double &xmin, double &xmax, double &ymin, double &ymax, double &zmin, double &zmax) const
{
    xmin = ymin = zmin = std::numeric_limits<double>::infinity();
    xmax = ymax = zmax = -std::numeric_limits<double>::infinity();

    for (int i = 0; i < m_shape->GetNumberOfVertices(); i++)
    {
        ColumnVector tp = transformation->InverseTransformPoint(m_shape->GetPoint(i));

        if (tp(1) < xmin)
            xmin = tp(1);

        if (tp(1) > xmax)
            xmax = tp(1);

        if (tp(2) < ymin)
            ymin = tp(2);

        if (tp(2) > ymax)
            ymax = tp(2);

        if (tp(3) < zmin)
            zmin = tp(3);

        if (tp(3) > zmax)
            zmax = tp(3);
    }
}

ColumnVector ShapeModel::SampleProfile(int vertex, int profilepoints, double profilespacing,
                                       const NEWIMAGE::volume<double> &vol, const boost::shared_ptr<const Transformation> &transformation) const
{
    // Do we need to set an interpolation method? Probably assume it has been set before ..

    BOOST_LOG_TRIVIAL(debug) << "Sampling profile at vertex " << vertex;

    ColumnVector voxelsize(3);
    voxelsize << vol.xdim() << vol.ydim() << vol.zdim();

    ColumnVector point = m_shape->GetPoint(vertex);
    ColumnVector normal = m_shape->GetNormal(vertex);

    ColumnVector result(profilepoints);
    double cent = (profilepoints - 1) * profilespacing / 2.0;

    for (int i = 0; i < profilepoints; i++)
    {
        ColumnVector p = point + (i * profilespacing - cent) * normal;
        ColumnVector tp = transformation->InverseTransformPoint(p);

        result(i + 1) = vol.interpolate(tp(1) / voxelsize(1), tp(2) / voxelsize(2), tp(3) / voxelsize(3));
    }

    return result;
}

void ShapeModel::WritePlots(const Plotting::plotfunc &pfunc) const
{
    if (m_vertexModels.size() != m_shape->GetNumberOfVertices())
        throw ModelException("Can't create plots before creating vertex models");

    for (int vert = 0; vert < m_shape->GetNumberOfVertices(); vert++)
        m_vertexModels[vert]->WritePlots(pfunc, vert);
}

std::vector<std::vector<double> > ShapeModel::GetAllDeltaLikelihoods(const std::unordered_map<string, std::vector<ColumnVector> > &data)
{
    BOOST_LOG_TRIVIAL(info) << "Getting likelihoods for deltas";

    std::vector<std::vector<double> > profilelogprobs;
    for (int v = 0; v < m_shape->GetNumberOfVertices(); v++)
    {
        std::unordered_map<std::string, ColumnVector> vertexdata;
        for (const auto &t : data)
            vertexdata[t.first] = t.second[v];

        BOOST_LOG_TRIVIAL(debug) << "Calling GetDeltaLikelihoods for vertex " << v;
        profilelogprobs.push_back(m_vertexModels[v]->GetDeltaLikelihoods(vertexdata, false));
    }

    return profilelogprobs;
}

boost::shared_ptr<ProfileModel> ShapeModel::TrainVertexImplementation(
        int vert, const std::unordered_map<string, std::vector<std::vector<ColumnVector> > > &trainingdata)
{
    BOOST_LOG_TRIVIAL(debug) << "Training vertex " << vert << " (default implementation)";

    boost::shared_ptr<ProfileModel> pm = m_vertexModels[vert];

    std::unordered_map<std::string, std::vector<ColumnVector> > verttd;

    // NB: This is an additional normalisation step to make scaling between modalities consistent to improve numerical stability.
    // It is transparent and unrelated to volume normalisation, which is what is controlled using the CLI options
    for (const auto &td : trainingdata)
    {
        for (const auto &v : td.second[vert])
            verttd[td.first].push_back(v);

        ColumnVector mean(m_profilePoints);
        mean = 0.0;
        for (const auto &vd : td.second)
            mean += Stats::mean(vd) / m_shape->GetNumberOfVertices();

        double scalarmean = mean.Sum() / m_profilePoints;

        double ssq = 0.0;
        for (const auto &vd : td.second)
            for (const auto &sd : vd)
                ssq += (sd - scalarmean).SumSquare();

        double scalarstdev = std::sqrt(ssq / m_profilePoints / m_shape->GetNumberOfVertices()
                                         / td.second.cbegin()->cbegin()->Nrows());

        BOOST_LOG_TRIVIAL(info) << "Normalisation of " << td.first << " -> Mean: " << scalarmean
                                   << ", Stdev: " << scalarstdev <<  " (used during optimisation only)";

        pm->SetNormalisation(td.first, scalarmean, scalarstdev);
    }

    pm->Train(verttd);

    return pm;
}

std::vector<int> ShapeModel::UpdateDeltas(const std::vector<int> &current,
                                         const std::vector<std::vector<double> > &profilelogprobs,
                                         int &diffs) const
{
    std::vector<int> updated(current);
    diffs = 0;

    for (int v = 0; v < current.size(); v++)
    {
        std::vector<double> logprobs = GetConditional(updated, v);

        auto plp = profilelogprobs[v].cbegin();
        for (auto lp = logprobs.begin(); lp != logprobs.end();)
            *lp++ += *plp++;

        updated[v] = std::max_element(logprobs.cbegin(), logprobs.cend()) - logprobs.cbegin();
        diffs += std::abs(updated[v] - current[v]);
    }

    return updated;
}

ColumnVector ShapeModel::FitImplementation(const std::unordered_map<std::string, std::vector<ColumnVector> > &data,
                                           ColumnVector &cilower, ColumnVector &ciupper, double cilevel)
{
    std::vector<std::vector<double> > profilelogprobs = GetAllDeltaLikelihoods(data);

    int verts = profilelogprobs.size();
    int steps = profilelogprobs[0].size();

    std::vector<int> current(verts);
    for (int i = 0; i < verts; i++)
        current[i] = std::max_element(profilelogprobs[i].cbegin(), profilelogprobs[i].cend()) - profilelogprobs[i].cbegin();

    int iter = 0;
    int diffs;
    do
    {
        current = UpdateDeltas(current, profilelogprobs, diffs);

        BOOST_LOG_TRIVIAL(info) << "ICM iteration " << iter << ", sum of absolute differences: " << diffs;

        if (++iter == m_maxiter)
        {
            BOOST_LOG_TRIVIAL(warning) << "Maximum number of ICM iterations reached";
            break;
        }
    }
    while (diffs);

    ColumnVector result(verts);
    for (int i = 0; i < verts; i++)
        result(i + 1) = m_profileSpacing * ((steps - 1) / 2.0 - current[i]);

    cilower = ColumnVector(verts);
    cilower = - (steps - 1) / 2.0 * m_profileSpacing;
    ciupper = ColumnVector(verts);
    ciupper = (steps - 1) / 2.0 * m_profileSpacing;

    return result;
}
