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

#include "mrfshapemodel.h"
#include <boost/log/trivial.hpp>
#include <algorithm>
#include <numeric>

using namespace NEWMAT;

BOOST_CLASS_EXPORT_IMPLEMENT(MRFShapeModel)

MRFShapeModel::MRFShapeModel(const Shape &shape, const std::vector<std::string> &modalitynames, int profilepoints, double profilespacing,
                             bool usenormalisationmasks) :
    ShapeModel(shape, modalitynames, profilepoints, profilespacing, usenormalisationmasks),
    m_weight(1.0),
    m_meanFraction(0.5),
    m_mean(shape.GetNumberOfVertices(), 0)
{
    BOOST_LOG_TRIVIAL(debug) << "Initialising MRF shape model";
}

void MRFShapeModel::SetWeight(double weight)
{
    m_weight = weight;
}

void MRFShapeModel::SetMeanFraction(double meanfrac)
{
    m_meanFraction = meanfrac;
}

double MRFShapeModel::TrianglePotential(int delta1, int delta2, int delta3) const
{
    double trimean = (delta1 + delta2 + delta3) / 3.0;

    double d1 = m_profileSpacing * (delta1 - trimean);
    double d2 = m_profileSpacing * (delta2 - trimean);
    double d3 = m_profileSpacing * (delta3 - trimean);

    return d1 * d1 + d2 * d2 + d3 * d3;
}

double MRFShapeModel::MeanPotential(int delta, int mean) const
{
    double dm = m_profileSpacing * (delta - mean);

    return dm * dm;
}

std::vector<int> MRFShapeModel::UpdateMean(const std::vector<int> &currentmean,
                                           const std::vector<std::vector<int> > &currentdeltas,
                                           int &diffs) const
{
    // NB: The updated mean doesn't actually depend on the current one, but we need the latter to find the number of changes
    std::vector<int> updatedmean(currentmean.size(), 0.0);
    int steps = m_vertexModels[0]->GetNumberOfSteps();
    diffs = 0;

    for (int vert = 0; vert < currentmean.size(); vert++)
    {
        double maxlogprob = - std::numeric_limits<double>::infinity();

        for (int d = 0; d < steps; d++)
        {
            double logprob = 0.0;

            for (int i = 0; i < currentdeltas.size(); i++)
                logprob -= m_meanFraction * m_weight * MeanPotential(d, currentdeltas[i][vert]);

            if (logprob > maxlogprob)
            {
                maxlogprob = logprob;
                updatedmean[vert] = d;
            }
        }
        
        diffs += std::abs(updatedmean[vert] - currentmean[vert]);
    }

    return updatedmean;
}

void MRFShapeModel::TrainShapeImplementation(const std::unordered_map<string, std::vector<std::vector<ColumnVector> > > &trainingdata)
{
    int verts = m_shape->GetNumberOfVertices();
    m_mean = std::vector<int>(verts, m_vertexModels[0]->GetNumberOfSteps() / 2);

    if (m_meanFraction)
    {
        BOOST_LOG_TRIVIAL(info) << "Training MRF shape model";

        int zsize = trainingdata.cbegin()->second.cbegin()->size();

        BOOST_LOG_TRIVIAL(info) << "Getting likelihoods for deltas (" << zsize << " subjects)";

        std::vector<std::vector<std::vector<double> > > traininglogprobs(zsize);
        for (int vert = 0; vert < verts; vert++)
        {
            BOOST_LOG_TRIVIAL(info) << "Vertex " << vert;

            for (int i = 0; i < zsize; i++)
            {
                std::unordered_map<std::string, ColumnVector> data;
                for (const auto &t : trainingdata)
                    data[t.first] = t.second[vert][i];

                // TODO: Might be better to set usedeltaprior to false here
                traininglogprobs[i].push_back(m_vertexModels[vert]->GetDeltaLikelihoods(data, true));
            }
        }

        BOOST_LOG_TRIVIAL(info) << "Fitting shape model";

        std::vector<std::vector<int> > current(zsize, std::vector<int>(verts, 0));
        for (int i = 0; i < zsize; i++)
            for (int vert = 0; vert < verts; vert++)
                current[i][vert] = std::max_element(traininglogprobs[i][vert].cbegin(), traininglogprobs[i][vert].cend())
                                    - traininglogprobs[i][vert].cbegin();

        int totaldiffs;
        int iter = 0;

        do
        {
            totaldiffs = 0;

            for (int i = 0; i < zsize; i++)
            {
                int ddiffs;
                current[i] = UpdateDeltas(current[i], traininglogprobs[i], ddiffs);
                totaldiffs += ddiffs;
            }

            int mdiffs;
            m_mean = UpdateMean(m_mean, current, mdiffs);
            totaldiffs += mdiffs;

            BOOST_LOG_TRIVIAL(info) << "Training shape model; iteration " << iter
                                    << ", sum of absolute differences: " << totaldiffs << " (mean only: " << mdiffs << ")";

            // TODO: Add option!
            if (++iter == 100)
            {
                BOOST_LOG_TRIVIAL(warning) << "Maximum number of iterations reached when training shape model";
                break;
            }
        }
        while (totaldiffs);
    }
    else
        BOOST_LOG_TRIVIAL(info) << "No training required for MRF shape model as mean fraction is zero";
}

std::vector<double> MRFShapeModel::GetConditional(const std::vector<int> &deltas, int vert) const
{
    // TODO: Change name of mean to something more meaningful?

    int steps = m_vertexModels[0]->GetNumberOfSteps();
    const std::vector<std::vector<vtkIdType> > tris = m_shape->GetPolysForVertex(vert);
    double neighbourweight = (1.0 - m_meanFraction) * m_weight;

    std::vector<double> logprobs(steps, 0.0);

    for (int d = 0; d < steps; d++)
    {
        logprobs[d] -= m_meanFraction * m_weight * MeanPotential(d, m_mean[vert]);

        for (const auto &tri : tris)
        {
            if (tri.size() != 3)
                throw ModelException("The ICM implementation assumes that all polygons in the mesh are triangles, but this is not the case");

            if (tri[0] == vert)
                logprobs[d] -= neighbourweight * TrianglePotential(d, deltas[tri[1]], deltas[tri[2]]);
            else if (tri[1] == vert)
                logprobs[d] -= neighbourweight * TrianglePotential(d, deltas[tri[0]], deltas[tri[2]]);
            else
                logprobs[d] -= neighbourweight * TrianglePotential(d, deltas[tri[0]], deltas[tri[1]]);
        }
    }

    return logprobs;
}
