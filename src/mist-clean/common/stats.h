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

#ifndef STATS_H
#define STATS_H

#include "newmat.h"
#include <vector>
#include <stdexcept>

namespace Stats
{
    double logsumexp(double a, double b);
    double logsumexp(const std::vector<double> &x);
    NEWMAT::ColumnVector mean(const std::vector<NEWMAT::ColumnVector> &x);
    NEWMAT::Matrix cov(const std::vector<NEWMAT::ColumnVector> &x);
    NEWMAT::Matrix pinv(const NEWMAT::Matrix &x, double tolfrac, int &rank);

    double lognormal(double x, double mu, double lambda);
    double logstudent(double x, double mu, double lambda, double alpha);

    class StatsException : public std::logic_error
    {
    public:
        StatsException(const std::string &what) : logic_error(what) { }
    };

    class VectorDistribution
    {
    public:
        double P(const NEWMAT::ColumnVector &x) const;
        double LogP(const NEWMAT::ColumnVector &x) const;

        virtual ~VectorDistribution() { }

    protected:
        // Length of x
        int m_k = 0;

        virtual double ComputeP(const NEWMAT::ColumnVector &x) const;
        virtual double ComputeLogP(const NEWMAT::ColumnVector &x) const = 0;
    };

    class MatrixDistribution
    {
    public:
        double P(const NEWMAT::Matrix &x) const;
        double LogP(const NEWMAT::Matrix &x) const;

        virtual ~MatrixDistribution() { }

    protected:
        // Length of x
        int m_k = 0;

        virtual double ComputeP(const NEWMAT::Matrix &x) const;
        virtual double ComputeLogP(const NEWMAT::Matrix &x) const = 0;
    };

    class MVN : public VectorDistribution
    {
    public:
        MVN(const NEWMAT::ColumnVector &mu, const NEWMAT::Matrix &lambda);
        virtual double ComputeLogP(const NEWMAT::ColumnVector &x) const;

    private:
        NEWMAT::ColumnVector m_mu;
        NEWMAT::Matrix m_lambda;
        double m_logConstant;
    };

    class MVStudent : public VectorDistribution
    {
    public:
        MVStudent(const NEWMAT::ColumnVector &mu, const NEWMAT::Matrix &lambda, double alpha);
        virtual double ComputeLogP(const NEWMAT::ColumnVector &x) const;

    private:
        NEWMAT::ColumnVector m_mu;
        NEWMAT::Matrix m_lambda;
        double m_alpha;
        double m_logConstant;
    };

    class Wishart : public MatrixDistribution
    {
    public:
        Wishart(double alpha, const NEWMAT::Matrix &beta);
        virtual double ComputeLogP(const NEWMAT::Matrix &x) const;

    private:
        double m_alpha;
        NEWMAT::Matrix m_beta;
        double m_logConstant;
    };

    class MVNWishart
    {
    public:
        MVNWishart(const NEWMAT::ColumnVector &mu, double lambda, double alpha, const NEWMAT::Matrix &beta);
        double P(const NEWMAT::ColumnVector &x, const NEWMAT::Matrix &y) const;
        double LogP(const NEWMAT::ColumnVector &x, const NEWMAT::Matrix &y) const;

    private:
        NEWMAT::ColumnVector m_mu;
        double m_lambda;
        Wishart m_wishart;
    };

    class Dirichlet : public VectorDistribution
    {
    public:
        Dirichlet(const NEWMAT::ColumnVector &alpha);
        virtual double ComputeLogP(const NEWMAT::ColumnVector &x) const;

    private:
        NEWMAT::ColumnVector m_alpha;
        double m_logConstant;
    };
}
#endif // STATS_H
