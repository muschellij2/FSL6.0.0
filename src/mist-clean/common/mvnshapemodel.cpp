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

#include "newmat.h"

#include "mvnshapemodel.h"

#include "stats.h"
#include <boost/log/trivial.hpp>
#include <cstring>

#define WANT_STREAM
#include "newmatio.h"

BOOST_CLASS_EXPORT_IMPLEMENT(MVNShapeModel)

using namespace NEWMAT;

MVNShapeModel::MVNShapeModel(const Shape &shape, const std::vector<std::string> &modalitynames, int profilepoints, double profilespacing,
                             bool usenormalisationmasks)
    : GibbsShapeModel(shape, modalitynames, profilepoints, profilespacing, usenormalisationmasks)
{
    BOOST_LOG_TRIVIAL(debug) << "Initialising MVN shape model";

    m_alpha = (shape.GetNumberOfVertices() - 1.0) / 2.0 + 1.0;
    m_n0 = 1;
    SetShapeCovariancePrior(1.0, 1.0);
}

void MVNShapeModel::SetShapeN0(int n0)
{
    m_n0 = n0;
}

void MVNShapeModel::SetShapeAlpha(double alpha)
{
    double alphamin = (m_shape->GetNumberOfVertices() - 1) / 2.0 + 1.0;

    if (alpha < alphamin)
        throw ModelException(std::string("Alpha should be at least ") + std::to_string(alphamin));

    m_alpha = alpha;
}

void MVNShapeModel::SetShapeCovariancePrior(double stdev, double smoothness)
{
    m_betaStdev = stdev;
    m_betaSmoothness = smoothness;
}

Matrix MVNShapeModel::GetBeta() const
{
    BOOST_LOG_TRIVIAL(debug) << "Building smoothness prior with width (smoothness) " << m_betaSmoothness
                             << " and height (stdev) " << m_betaStdev;

    Matrix D = m_shape->GetDistanceMatrix();

    Matrix G(D.Nrows(), D.Ncols());
    for (int i = 1; i <= D.Nrows(); i++)
        for (int j = 1; j <= D.Nrows(); j++)
            G(i, j) = std::exp(- std::pow(D(i, j) / m_betaSmoothness, 2.0) / 2.0);

    Matrix beta = m_alpha * m_betaStdev * m_betaStdev * G;

    return beta;
}

void MVNShapeModel::UseShapeVariability(bool use)
{
    m_useShapeVariability = use;
}

void MVNShapeModel::TrainShapeImplementation(const std::unordered_map<string, std::vector<std::vector<ColumnVector> > > &trainingdata)
{
    int verts = m_shape->GetNumberOfVertices();

    Matrix beta = GetBeta();

    BOOST_LOG_TRIVIAL(info) << "Training MVN shape model with parameters:\n"
            << "\tUseShapeVariability = " << m_useShapeVariability << "\n"
            << "\tN0 = " << m_n0 << "\n"
            << "\tAlpha = " << m_alpha << "\n"
            << "\tCovariance prior stdev = " << m_betaStdev << "\n"
            << "\tCovariance prior smoothness = " << m_betaSmoothness << "\n"
            << "\tBeta diagonal = " << beta(1, 1)
            << " ... " << beta(verts / 2, verts / 2)
            << " ... " << beta(verts, verts) << "\n";

    int zsize = 0;
    ColumnVector zmean(m_shape->GetNumberOfVertices());
    zmean = 0.0;
    Matrix zcov(m_shape->GetNumberOfVertices(), m_shape->GetNumberOfVertices());
    zcov = 0.0;

    if (m_useShapeVariability)
    {
        zsize = trainingdata.cbegin()->second.cbegin()->size();
        int steps = m_vertexModels[0]->GetNumberOfSteps();

        std::vector<ColumnVector> trainingdisplacements(zsize, ColumnVector(m_shape->GetNumberOfVertices()));
        for (int vert = 0; vert < m_shape->GetNumberOfVertices(); vert++)
        {
            for (int i = 0; i < zsize; i++)
            {
                std::unordered_map<std::string, ColumnVector> data;
                for (const auto &t : trainingdata)
                    data[t.first] = t.second[vert][i];

                std::vector<double> logprobs = m_vertexModels[vert]->GetDeltaLikelihoods(data, true);

                trainingdisplacements[i](vert + 1) = m_profileSpacing
                        * ((steps - 1) / 2.0 - (std::max_element(logprobs.cbegin(), logprobs.cend()) - logprobs.cbegin()));
            }
        }

        // As mentioned in class definition, mu0 is zero.
        zmean = Stats::mean(trainingdisplacements);
        zcov = Stats::cov(trainingdisplacements);
    }

    m_stMu = zsize * zmean / (m_n0 + zsize);

    double alphan = m_alpha + 0.5 * (zsize - m_shape->GetNumberOfVertices() + 1);
    Matrix betan = beta + 0.5 * ((zsize - 1.0) * zcov + m_n0 * zsize / (m_n0 + zsize) * (zmean * zmean.t()));

    m_stLambda = (m_n0 + zsize) / (m_n0 + zsize + 1.0) * alphan * betan.i();
    m_stAlpha = 2 * alphan;
}

std::vector<double> MVNShapeModel::GetConditionalContinuous(const ColumnVector &x, int vert) const
{
    // See Bernardo+Smith p140
    int k = m_shape->GetNumberOfVertices();
    int steps = m_vertexModels[0]->GetNumberOfSteps();
    double condalpha = m_stAlpha + k - 1;

    ColumnVector d(k - 1);
    for (std::size_t i = 1; i < vert + 1; i++)
        d(i) = x(i) - m_stMu(i);

    for (std::size_t i = vert + 1; i < k; i++)
        d(i) = x(i + 1) - m_stMu(i + 1);

    RowVector lambda12(k - 1);

    for (int i = 0; i < vert; i++)
        lambda12(i + 1) = m_stLambda(vert + 1, i + 1);

    for (int i = vert; i < k - 1; i++)
        lambda12(i + 1) = m_stLambda(vert + 1, i + 2);

    double lambda11 = m_stLambda(vert + 1, vert + 1);
    // This is d.t() * lambda11Schur:
    RowVector dtL(k - 1);
    dtL = 0.0;

    // TODO: Exploit symmetry!

    for (int i = 1; i <= vert; i++)
        for (int j = 1; j <= vert; j++)
            dtL(j) += d(i) * (m_stLambda(i, j) - lambda12(i) * lambda12(j) / lambda11);

    for (int i = vert; i < k - 1; i++)
        for (int j = 0; j < vert; j++)
            dtL(j + 1) += d(i + 1) * (m_stLambda(i + 2, j + 1) - lambda12(i + 1) * lambda12(j + 1) / lambda11);

    for (int i = 0; i < vert; i++)
        for (int j = vert; j < k - 1; j++)
            dtL(j + 1) += d(i + 1) * (m_stLambda(i + 1, j + 2) - lambda12(i + 1) * lambda12(j + 1) / lambda11);

    for (int i = vert; i < k - 1; i++)
        for (int j = vert; j < k - 1; j++)
            dtL(j + 1) += d(i + 1) * (m_stLambda(i + 2, j + 2) - lambda12(i + 1) * lambda12(j + 1) / lambda11);

    double condmu = m_stMu(vert + 1) - (lambda12 * d).AsScalar() / lambda11;
    double condlambda = lambda11 * condalpha / (m_stAlpha + (dtL * d).AsScalar());

    BOOST_LOG_TRIVIAL(debug) << "condmu = " << condmu;
    BOOST_LOG_TRIVIAL(debug) << "condlambda = " << condlambda;

    std::vector<double> logprobs(steps);
    for (int delta = 0; delta < steps; delta++)
    {
        logprobs[delta] = Stats::logstudent(m_profileSpacing * ((steps - 1) / 2.0 - delta), condmu, condlambda, condalpha);

        BOOST_LOG_TRIVIAL(debug) << "delta = " << delta << " --> shape: " << logprobs[delta];
    }

    return logprobs;
}
