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

#ifndef MRFSHAPEMODEL_H
#define MRFSHAPEMODEL_H

#include "shapemodel.h"

class MRFShapeModel;
BOOST_CLASS_EXPORT_KEY(MRFShapeModel)

class MRFShapeModel : public ShapeModel
{
public:
    MRFShapeModel(const Shape &shape, const std::vector<std::string> &modalitynames, int profilepoints, double profilespacing,
                  bool usenormalisationmasks);

    void SetWeight(double weight);
    void SetMeanFraction(double meanfrac);

protected:
    virtual void TrainShapeImplementation(const std::unordered_map<std::string, std::vector<std::vector<NEWMAT::ColumnVector> > > &trainingdata);

    virtual std::vector<double> GetConditional(const std::vector<int> &deltas, int vert) const;

private:
    double m_weight;
    double m_meanFraction;
    std::vector<int> m_mean;

    double TrianglePotential(int delta1, int delta2, int delta3) const;
    double MeanPotential(int delta, int mean) const;

    std::vector<int> UpdateMean(const std::vector<int> &currentmean,
                                const std::vector<std::vector<int> > &currentdeltas,
                                int &diffs) const;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version);
};

template<class Archive>
void MRFShapeModel::serialize(Archive &ar, const unsigned int version)
{
    BOOST_LOG_TRIVIAL(debug) << "(De)serialising MRFShapeModel";

    ar & boost::serialization::base_object<ShapeModel>(*this);

    ar & m_weight;
    ar & m_meanFraction;
    ar & m_mean;
}


namespace boost
{
namespace serialization
{
    template<class Archive>
    inline void save_construct_data(Archive &ar, const MRFShapeModel *m, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving contruction info for MRFShapeModel";

        const std::vector<std::string> modalitynames = m->GetModalityNames();
        int points = m->GetProfilePoints();
        double spacing = m->GetProfileSpacing();
        bool usenormalisationmasks = m->GetUseNormalisationMasks();
        boost::shared_ptr<const Shape> shape = m->GetShape();

        ar << modalitynames;
        ar << points;
        ar << spacing;
        ar << usenormalisationmasks;
        ar << shape;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, MRFShapeModel *m, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for MRFShapeModel";

        std::vector<std::string> modalitynames;
        int profilepoints;
        double profilespacing;
        bool usenormalisationmasks;
        boost::shared_ptr<Shape> shape;

        ar >> modalitynames;
        ar >> profilepoints;
        ar >> profilespacing;
        ar >> usenormalisationmasks;
        ar >> shape;

        ::new(m)MRFShapeModel(*shape, modalitynames, profilepoints, profilespacing, usenormalisationmasks);
    }
}
}

#endif // MRFSHAPEMODEL_H

