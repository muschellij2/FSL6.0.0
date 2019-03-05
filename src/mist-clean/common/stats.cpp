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

#include "stats.h"
#include "newmatap.h"
#include <cmath>
#include <algorithm>
#include <boost/log/trivial.hpp>

using namespace NEWMAT;

namespace Stats
{
    double logsumexp(double a, double b)
    {
        if (a > b)
            return a + log(1 + exp(b - a));
        else
            return b + log(1 + exp(a - b));
    }

    double logsumexp(const vector<double> &x)
    {
        double max = *std::max_element(x.cbegin(), x.cend());
        double s = 0.0;

        for (auto ix : x)
            s += exp(ix - max);

        return max + log(s);
    }

    ColumnVector mean(const std::vector<ColumnVector> &x)
    {
        int k = x[0].Nrows();
        ColumnVector sum(k);
        sum = 0.0;

        for (auto ix : x)
        {
            if (ix.Nrows() != k)
                throw StatsException("All vectors should have the same length");

            sum += ix;
        }

        return sum / x.size();
    }

    Matrix cov(const std::vector<ColumnVector> &x)
    {
        int k = x[0].Nrows();
        ColumnVector m = Stats::mean(x);

        Matrix ssq(k, k);
        ssq = 0.0;

        for (const auto &ix : x)
        {
            if (ix.Nrows() != k)
                throw StatsException("All vectors should have the same length");

            ColumnVector d = ix - m;
            ssq += d * d.t();
        }

        return ssq / (x.size() - 1);
    }

    Matrix pinv(const Matrix &x, double tolfrac, int &rank)
    {
        Matrix U, V;
        DiagonalMatrix D;

        SVD(x, D, U, V);
        double tol = tolfrac * D(1, 1);
        rank = x.Ncols();

        for (int i = 1; i <= x.Ncols(); i++)
        {
            if (D(i, i) >= tol)
                D(i, i) = 1.0 / D(i, i);
            else
            {
                D(i, i) = 0.0;
                rank--;
            }
        }

        return V * D * U.t();
    }

    double lognormal(double x, double mu, double lambda)
    {
        double c = 0.5 * std::log(lambda) - 0.5 * std::log(2 * M_PI);
        double d = x - mu;

        return c - 0.5 * lambda * d * d;
    }

    double logstudent(double x, double mu, double lambda, double alpha)
    {
        double c = std::lgamma(0.5 * (alpha + 1)) - std::lgamma(0.5 * alpha) + 0.5 * std::log(lambda / alpha / M_PI);
        double d = x - mu;

        return c - 0.5 * (alpha + 1) * std::log(1 + lambda / alpha * d * d);
    }

    double VectorDistribution::P(const ColumnVector &x) const
    {
        if (x.Nrows() != m_k)
            throw StatsException("Vector x has incorrect number of dimensions");

        return ComputeP(x);
    }

    double VectorDistribution::LogP(const ColumnVector &x) const
    {
        if (x.Nrows() != m_k)
            throw StatsException("Vector x has incorrect number of dimensions");

        return ComputeLogP(x);
    }

    double VectorDistribution::ComputeP(ColumnVector const &x) const
    {
        return std::exp(ComputeLogP(x));
    }

    double MatrixDistribution::P(const Matrix &x) const
    {
        if (x.Nrows() != m_k || x.Ncols() != m_k)
            throw StatsException("Matrix x should have same dimensions as beta");

        return ComputeP(x);
    }

    double MatrixDistribution::LogP(const Matrix &x) const
    {
        if (x.Nrows() != m_k || x.Ncols() != m_k)
            throw StatsException("Matrix x should have same dimensions as beta");

        return ComputeLogP(x);
    }

    double MatrixDistribution::ComputeP(const Matrix &x) const
    {
        return std::exp(ComputeLogP(x));
    }

    MVN::MVN(const NEWMAT::ColumnVector &mu, const NEWMAT::Matrix &lambda) :
        m_mu(mu),
        m_lambda(lambda)
    {
        if (mu.Nrows() != lambda.Nrows())
            throw StatsException("Number of dimensions of mu and lambda do not agree");

        if (lambda.Nrows() != lambda.Nrows())
            throw StatsException("Lambda is not square");

        m_k = mu.Nrows();

        m_logConstant = 0.5 * lambda.LogDeterminant().LogValue() - m_k / 2.0 * std::log(2 * M_PI);
    }

    double MVN::ComputeLogP(const ColumnVector &x) const
    {
        ColumnVector diff = x - m_mu;

        return m_logConstant - 0.5 * (diff.t() * m_lambda * diff).AsScalar();
    }

    MVStudent::MVStudent(const ColumnVector &mu, const Matrix &lambda, double alpha) :
        m_mu(mu),
        m_lambda(lambda),
        m_alpha(alpha)
    {
        if (mu.Nrows() != lambda.Nrows())
            throw StatsException("Number of dimensions of mu and lambda do not agree");

        if (lambda.Nrows() != lambda.Nrows())
            throw StatsException("Lambda is not square");

        m_k = mu.Nrows();

        m_logConstant = std::lgamma(0.5 * (alpha + m_k)) - std::lgamma(0.5 * alpha) - 0.5 * m_k * std::log(alpha * M_PI)
                + 0.5 * lambda.LogDeterminant().LogValue();
    }

    double MVStudent::ComputeLogP(const ColumnVector &x) const
    {
        ColumnVector diff = x - m_mu;
        return m_logConstant - 0.5 * (m_alpha + m_k) * std::log(1 + 1 / m_alpha * (diff.t() * m_lambda * diff).AsScalar());
    }

    Wishart::Wishart(double alpha, const Matrix &beta) :
        m_alpha(alpha),
        m_beta(beta)
    {
        if (beta.Nrows() != beta.Nrows())
            throw StatsException("Beta is not square");

        m_k = beta.Nrows();

        m_logConstant = -m_k * (m_k - 1) / 4.0 * std::log(M_PI) + alpha * beta.LogDeterminant().LogValue();
        for (int i = 1; i <= m_k; i++)
            m_logConstant -= std::lgamma(0.5 * (2 * alpha + 1 - i));
    }

    double Wishart::ComputeLogP(const Matrix &x) const
    {
        return m_logConstant + (m_alpha - (m_k + 1) / 2.0) * x.LogDeterminant().LogValue()
                - SP(m_beta, x).Sum();
    }

    MVNWishart::MVNWishart(const ColumnVector& mu, double lambda, double alpha, const Matrix& beta) :
        m_mu(mu),
        m_lambda(lambda),
        m_wishart(alpha, beta)
    {

    }

    double MVNWishart::P(const ColumnVector &x, const Matrix &y) const
    {
        return std::exp(LogP(x, y));
    }

    double MVNWishart::LogP(const ColumnVector &x, const Matrix &y) const
    {
        return MVN(m_mu, m_lambda * y).LogP(x) + m_wishart.LogP(y);
    }

    Dirichlet::Dirichlet(const ColumnVector &alpha) :
        m_alpha(alpha)
    {
        if (alpha.Nrows() < 1)
            throw StatsException("Alpha should have at least one element");

        m_k = alpha.Nrows() - 1;
        m_logConstant = std::lgamma(alpha.Sum());
        for (int i = 1; i <= alpha.Nrows(); i++)
            m_logConstant -= std::lgamma(alpha(i));
    }

    double Dirichlet::ComputeLogP(const ColumnVector &x) const
    {
        double val = m_logConstant + (m_alpha(m_k + 1) - 1) * std::log(1 - x.Sum());
        for (int i = 1; i <= m_k; i++)
            val += (m_alpha(i) - 1.0) * std::log(x(i));

        return val;
    }
}
