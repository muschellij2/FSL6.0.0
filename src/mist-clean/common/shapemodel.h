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

#ifndef SHAPEMODEL_H
#define SHAPEMODEL_H

#include "serialisation.h"
#include "shape.h"
#include "newimage.h"
#include "profilemodel.h"
#include "profilefilters.h"
#include "transformation.h"
#include <boost/shared_ptr.hpp>

#include <unordered_map>
#include <vector>
#include <memory>
#include <random>

class ShapeModel
{
public:
    class ModelException : public std::logic_error
    {
    public:
        ModelException(const std::string &what) : logic_error(what) { };
    };

    enum class NormalisationMode
    {
        None,
        Additive,
        Multiplicative
    };

    ShapeModel(const Shape &shape, const std::vector<std::string> &modalitynames, int profilepoints, double profilespacing,
               bool usenormalisationmasks);
    virtual ~ShapeModel();

    std::vector<std::string> GetModalityNames() const;
    int GetProfilePoints() const;
    double GetProfileSpacing() const;
    boost::shared_ptr<const Shape> GetShape() const;
    boost::shared_ptr<ProfileModel> GetVertexModel(int vertex) const;
    std::unordered_map<std::string, double> GetMeanIntensities() const;
    void SetMaxIterations(int maxiter);

    bool GetUseNormalisationMasks() const;
    void SetNormalise(std::string modality, NormalisationMode mode);
    NormalisationMode GetNormalise(std::string modality) const;
    NEWMAT::ColumnVector NormaliseProfile(std::string modality, const NEWMAT::ColumnVector &profile, double mean) const;

    void SetFilter(std::string modality, boost::shared_ptr<ProfileFilter> filter);
    boost::shared_ptr<ProfileFilter> GetFilter(std::string modality) const;

    void CreateVertexModels(boost::shared_ptr<const ProfileModel> example);

    void Train(const std::unordered_map<std::string, std::vector<std::string> > &trainingdata,
               const std::vector<boost::shared_ptr<const Transformation> > &transformations,
               const std::vector<std::string> &normalisationmasks);

    NEWMAT::ColumnVector Fit(const std::unordered_map<std::string, std::string> &data,
                             const boost::shared_ptr<const Transformation> &transformation,
                             const std::string &normalisationmask,
                             NEWMAT::ColumnVector &cilower, NEWMAT::ColumnVector &ciupper,
                             double cilevel);

    std::vector<boost::shared_ptr<ProfileModel> > TrainVertices(const std::unordered_map<string, std::vector<string> > &trainingdata,
                                                                const std::vector<boost::shared_ptr<const Transformation> > &transformations,
                                                                const std::vector<std::string> &normalisationmasks,
                                                                std::vector<int> verts,
                                                                const Plotting::plotfunc &pfunc = nullptr);

    void LoadVertexModelsAndTrain(const std::vector<std::string> &archivefiles,
                                  const std::unordered_map<string, std::vector<string> > &trainingdata,
                                  const std::vector<boost::shared_ptr<const Transformation> > &transformations,
                                  const std::vector<std::string> &normalisationmasks);

    std::unordered_map<std::string, std::vector<std::vector<ColumnVector> > > SampleAllVertices(
            const std::unordered_map<std::string, std::vector<std::string> > &vols,
            const std::vector<boost::shared_ptr<const Transformation> > &transformations,
            bool computemeans,
            const std::vector<std::string> &normalisationmasks);

    virtual void WritePlots(const Plotting::plotfunc &pfunc) const;

    // TODO: Make these protected again and sort out normalisation in PlotWindow
    NEWMAT::ColumnVector SampleProfile(int vertex, int profilepoints, double profilespacing,
                                       const NEWIMAGE::volume<double> &vol, const boost::shared_ptr<const Transformation> &transformation) const;

    double GetMeanIntensity(NEWIMAGE::volume<double> vol, NEWIMAGE::volume<double> *exclusionmask,
                            double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) const;
    void GetShapeExtent(const boost::shared_ptr<const Transformation> &transformation,
                        double &xmin, double &xmax, double &ymin, double &ymax, double &zmin, double &zmax) const;

protected:
    boost::shared_ptr<Shape> m_shape;

    std::vector<std::string> m_modalityNames;

    int m_profilePoints = 0;
    double m_profileSpacing = 0;

    int m_maxiter = 1;

    std::vector<boost::shared_ptr<ProfileModel> > m_vertexModels;

    void CheckModality(const std::string &modalityname) const;
    void CheckFilenames(const std::unordered_map<std::string, std::vector<std::string> > &trainingdata,
                                    const std::vector<boost::shared_ptr<const Transformation> > &transformations,
                                    const std::vector<std::string> &normalisationmasks) const;

    std::vector<std::vector<double> > GetAllDeltaLikelihoods(const std::unordered_map<string, std::vector<NEWMAT::ColumnVector> > &data);

    std::vector<int> UpdateDeltas(const std::vector<int> &current,
                                  const std::vector<std::vector<double> > &profilelogprobs,
                                  int &diffs) const;

    virtual std::vector<double> GetConditional(const std::vector<int> &deltas, int vert) const = 0;

    virtual NEWMAT::ColumnVector FitImplementation(const std::unordered_map<std::string, std::vector<NEWMAT::ColumnVector> > &data,
                                                   NEWMAT::ColumnVector &cilower, NEWMAT::ColumnVector &ciupper, double cilevel);

    virtual boost::shared_ptr<ProfileModel> TrainVertexImplementation(int vert,
            const std::unordered_map<std::string, std::vector<std::vector<NEWMAT::ColumnVector> > > &trainingdata);

    virtual void TrainShapeImplementation(
            const std::unordered_map<std::string, std::vector<std::vector<NEWMAT::ColumnVector> > > &trainingdata) = 0;

private:
    std::unordered_map<std::string, NormalisationMode> m_normalise;
    std::unordered_map<std::string, double> m_meanIntensities;
    bool m_useNormalisationMasks;

    std::unordered_map<std::string, boost::shared_ptr<ProfileFilter> > m_filters;

    ShapeModel(const ShapeModel &);
    ShapeModel& operator=(const ShapeModel &);

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version);
};

template<class Archive>
void ShapeModel::serialize(Archive &ar, const unsigned int version)
{
    BOOST_LOG_TRIVIAL(debug) << "(De)serialising ShapeModel";
    // NB: Constructor parameters are serialised in non-virtual derived classes!
    ar & m_vertexModels;
    ar & m_normalise;
    ar & m_meanIntensities;
    ar & m_filters;
    ar & m_maxiter;
}


#endif // SHAPEMODEL_H
