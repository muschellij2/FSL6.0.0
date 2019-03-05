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

#include "gibbsshapemodel.h"

GibbsShapeModel::GibbsShapeModel(const Shape &shape, const std::vector<std::string> &modalitynames, int profilepoints, double profilespacing,
                                 bool usenormalisationmasks) :
    ShapeModel(shape, modalitynames, profilepoints, profilespacing, usenormalisationmasks)
{
    std::random_device rd;
    Reseed(rd());
}

void GibbsShapeModel::Reseed(std::mt19937::result_type seed)
{
    BOOST_LOG_TRIVIAL(info) << "Reseeding MT for shape model. Seed = " << seed;
    m_randomEngine.seed(seed);
}

void GibbsShapeModel::UseGibbs(bool use)
{
    m_useGibbs = use;
}

void GibbsShapeModel::SetGibbsIterations(int iterations)
{
    if (iterations <= m_gibbsBurnIn)
        throw ModelException("Total number of iterations must be larger than number of burn-in iterations");

    m_gibbsIters = iterations;
}

void GibbsShapeModel::SetGibbsBurnIn(int burnin)
{
    if (burnin >= m_gibbsIters)
        throw ModelException("Number of burn-in iterations must be smaller than total number of iterations");

    m_gibbsBurnIn = burnin;
}

ColumnVector GibbsShapeModel::FitImplementation(const std::unordered_map<string, std::vector<ColumnVector> > &data,
                                              ColumnVector &cilower, ColumnVector &ciupper, double cilevel)
{
    // TODO: Add m_useGibbs to class + serialization code
    if (m_useGibbs)
        return DoGibbs(data, cilower, ciupper, cilevel);

    return ShapeModel::FitImplementation(data, cilower, ciupper, cilevel);
}

std::vector<double> GibbsShapeModel::GetConditional(const std::vector<int> &deltas, int vert) const
{
    int steps = m_vertexModels[0]->GetNumberOfSteps();

    ColumnVector p(deltas.size());

    for (int i = 0; i < deltas.size(); i++)
        p(i + 1) = m_profileSpacing * ((steps - 1) / 2.0 - deltas[i]);

    return GetConditionalContinuous(p, vert);
}

ColumnVector GibbsShapeModel::DoGibbs(const std::unordered_map<string, std::vector<ColumnVector> > &data,
                                      ColumnVector &cilower, ColumnVector &ciupper, double cilevel)
{
    std::vector<std::vector<double> > profilelogprobs = GetAllDeltaLikelihoods(data);
    int steps = m_vertexModels[0]->GetNumberOfSteps();

    std::vector<ColumnVector> chain(m_gibbsIters - m_gibbsBurnIn);
    ColumnVector p(m_shape->GetNumberOfVertices());
    p = 0.0;

    for (int v = 0; v < m_shape->GetNumberOfVertices(); v++)
    {
        const std::vector<double> &dl = profilelogprobs[v];

        p(v + 1) = m_profileSpacing * ((steps - 1) / 2.0 - (std::max_element(dl.begin(), dl.end()) - dl.begin()));
        BOOST_LOG_TRIVIAL(debug) << "Using initial displacement " << p(v + 1) << " for vertex " << v;
    }

    for (int iter = 0; iter < m_gibbsIters; iter++)
    {
        BOOST_LOG_TRIVIAL(info) << "Gibbs iteration " << iter
                                << (iter < m_gibbsBurnIn ? " (burn-in)" : "") <<  " ...";

        for (int v = 0; v < m_shape->GetNumberOfVertices(); v++)
        {
            BOOST_LOG_TRIVIAL(debug) << "Old p(" << v + 1 << ") = " << p(v + 1) << " (Vertex " << v << ")";

            std::vector<double> logprobs = GetConditionalContinuous(p, v);

            for (auto lp = logprobs.begin(), plp = profilelogprobs[v].begin(); lp != logprobs.end();)
                *lp++ += *plp++;

            p(v + 1) = GenerateSample(logprobs);
            BOOST_LOG_TRIVIAL(debug) << "New p(" << v + 1 << ") = " << p(v + 1) << " (Vertex " << v << ")";
        }

        if (iter >= m_gibbsBurnIn)
            chain[iter - m_gibbsBurnIn] = p;
    }

    std::vector<std::vector<int> > hists(m_shape->GetNumberOfVertices(), std::vector<int>(steps, 0));
    for (const auto &t : chain)
    {
        for (int v = 0; v < m_shape->GetNumberOfVertices(); v++)
        {
            double val = - t(v + 1) / m_profileSpacing + (steps - 1) / 2.0;

            int bin = static_cast<int>(val + 0.5);

            if (bin < 0)
                bin = 0;

            if (bin >= steps)
                bin = steps - 1;

            hists[v][bin]++;
        }
    }

    ColumnVector displacements(m_shape->GetNumberOfVertices());
    ColumnVector cil(m_shape->GetNumberOfVertices());
    ColumnVector ciu(m_shape->GetNumberOfVertices());

    displacements = 0.0;
    for (int v = 0; v < m_shape->GetNumberOfVertices(); v++)
    {
        BOOST_LOG_TRIVIAL(debug) << "Histogram for vertex " << v << ":";

        int sum = 0;
        int maxval = 0;
        int maxind = 0;
        int ciuind = -1;
        int cilind = 0;

        for (std::size_t i = 0; i < steps; i++)
        {
            BOOST_LOG_TRIVIAL(debug) << "Bin " << i << ": " << hists[v][i];

            if (hists[v][i] > maxval)
            {
                maxval = hists[v][i];
                maxind = i;
            }

            sum += hists[v][i];

            if (sum < (m_gibbsIters - m_gibbsBurnIn) * (0.5 - cilevel / 2))
                ciuind = i;

            if (cilind == 0 && sum > (m_gibbsIters - m_gibbsBurnIn) * (0.5 + cilevel / 2))
                cilind = i;
        }

        displacements(v + 1) = m_profileSpacing * ((steps - 1) / 2.0 - maxind);
        cil(v + 1) = m_profileSpacing * ((steps - 1) / 2.0 - cilind - 0.5);
        ciu(v + 1) = m_profileSpacing * ((steps - 1) / 2.0 - ciuind - 0.5);
    }

    cilower = cil;
    ciupper = ciu;

    return displacements;
}

double GibbsShapeModel::GenerateSample(const std::vector<double> &logprobs)
{
    int steps = m_vertexModels[0]->GetNumberOfSteps();

    auto max = std::max_element(logprobs.cbegin(), logprobs.cend());
    std::vector<double> probs(logprobs.size());
    std::transform(logprobs.cbegin(), logprobs.cend(), probs.begin(), [&](double p){ return std::exp(p - *max); });
    std::vector<double> cumprobs(probs.size() + 1);
    cumprobs[0] = 0.0;
    std::partial_sum(probs.cbegin(), probs.cend(), cumprobs.begin() + 1);

    std::uniform_real_distribution<> dist(0.0, *(cumprobs.cend() - 1));
    double rn = dist(m_randomEngine);
    auto upper = std::find_if(cumprobs.cbegin(), cumprobs.cend(), [&](double p){ return p >= rn; });

    double sample;
    if (upper != cumprobs.cbegin())
    {
        auto lower = upper - 1;
        sample = (lower - cumprobs.cbegin()) + (rn - *lower) / (*upper - *lower);
    }
    else
        sample = 0.0;

    // Needed to center cumprobs probability regions around probs bin centres
    sample -= 0.5;

    return m_profileSpacing * ((steps - 1) / 2.0 - sample);
}
