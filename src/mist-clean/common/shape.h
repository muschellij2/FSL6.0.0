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

#ifndef SHAPE_H
#define SHAPE_H

#include "transformation.h"
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <memory>
#include <string>
#include <vector>
#include <boost/log/trivial.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include "newimage.h"

#include <newmat.h>


class Shape
{
public:
    class ShapeException : public std::runtime_error
    {
    public:
        ShapeException(const std::string& what) : runtime_error(what) { }
    };

    Shape(vtkSmartPointer<vtkPolyData> polydata, const std::string& displayname);
    Shape(const std::string& filename, const std::string& displayname);
    Shape(const std::string& filename, boost::shared_ptr<const Transformation> &transformation,
          const std::string& displayname);

    Shape(const Shape &other);
    Shape(const Shape &other, boost::shared_ptr<const Transformation> &transformation);

    Shape& operator=(const Shape &other);

    std::string GetDisplayName() const;
    vtkSmartPointer<vtkPolyData> GetPolyData() const;
    vtkIdType GetNumberOfVertices() const;
    NEWMAT::ColumnVector GetPoint(vtkIdType vertex) const;
    NEWMAT::ColumnVector GetNormal(vtkIdType vertex) const;
    NEWMAT::ColumnVector GetBounds() const;
    NEWMAT::ColumnVector GetCenter() const;

    std::vector<std::vector<vtkIdType> > GetPolysForVertex(vtkIdType vertex);

    // NB: Distances are in #edges, not mm!
    NEWMAT::Matrix GetDistanceMatrix() const;

    void SetDisplacement(vtkIdType vertex, double amount);
    double GetDisplacement(vtkIdType vertex) const;

    Shape RemoveSelfIntersections(double stepsize) const;

    void ToVolume(NEWIMAGE::volume<char> &vol, bool usevoxelvertices);

    void WritePolyData(const std::string &filename) const;

    virtual ~Shape() { };

protected:
    void Copy(const Shape &other);

private:
    // The VTK code and documentation are slightly unclear and inconsistent about when vtkPolyData::BuildLinks()
    // is called and when the links are invalidated, so I am adding this just to be sure
    bool m_linksBuilt;

    vtkSmartPointer<vtkPolyData> m_polyData;
    std::vector<double> m_displacements;
    std::string m_displayName;

    boost::shared_ptr<const Transformation> m_transformation;
    std::vector<NEWMAT::ColumnVector> m_untransformedPoints;
    std::vector<NEWMAT::ColumnVector> m_untransformedNormals;

    vtkSmartPointer<vtkPolyData> ReadPolyData(const std::string &filename);
    void Init(vtkSmartPointer<vtkPolyData> polydata);
    static vtkSmartPointer<vtkPolyData> ComputeInitialNormals(vtkSmartPointer<vtkPolyData>);

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version);
};

template<class Archive>
void Shape::serialize(Archive &ar, const unsigned int version)
{
    BOOST_LOG_TRIVIAL(debug) << "(De)serialising Shape";

    ar & m_displacements;
}

namespace boost
{
namespace serialization
{
    template<class Archive>
    inline void save_construct_data(Archive &ar, const Shape *m, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Saving construction info for Shape";

        vtkIdType npoints = m->GetPolyData()->GetNumberOfPoints();
        ar << npoints;
        for (vtkIdType i = 0; i < npoints; i++)
        {
            double p[3];
            m->GetPolyData()->GetPoint(i, p);
            ar << p;
        }

        vtkIdType npolys = m->GetPolyData()->GetNumberOfPolys();
        ar << npolys;
        auto idl = vtkSmartPointer<vtkIdList>::New();
        m->GetPolyData()->GetPolys()->InitTraversal();
        while (m->GetPolyData()->GetPolys()->GetNextCell(idl))
        {
            vtkIdType nids = idl->GetNumberOfIds();
            ar << nids;
            for (vtkIdType i = 0; i < nids; i++)
            {
                vtkIdType id = idl->GetId(i);
                ar << id;
            }
        }

        std::string displayname = m->GetDisplayName();
        ar << displayname;
    }

    template<class Archive>
    inline void load_construct_data(Archive &ar, Shape *m, const unsigned int version)
    {
        BOOST_LOG_TRIVIAL(debug) << "Loading construction info for Shape";

        BOOST_LOG_TRIVIAL(debug) << "Deserialising polydata";

        int npoints;
        ar >> npoints;
        BOOST_LOG_TRIVIAL(debug) << " ... points: " << npoints;
        auto points = vtkSmartPointer<vtkPoints>::New();
        points->SetNumberOfPoints(npoints);
        for (vtkIdType i = 0; i < npoints; i++)
        {
            double p[3];
            ar >> p;
            points->SetPoint(i, p);
        }

        int npolys;
        ar >> npolys;
        BOOST_LOG_TRIVIAL(debug) << " ... polys: " << npolys;
        auto polys = vtkSmartPointer<vtkCellArray>::New();
        for (int i = 0; i < npolys; i++)
        {
            vtkIdType nids;
            ar >> nids;
            auto idl = vtkSmartPointer<vtkIdList>::New();
            idl->SetNumberOfIds(nids);
            for (int j = 0; j < nids; j++)
            {
                vtkIdType id;
                ar >> id;
                idl->SetId(j, id);
            }

            polys->InsertNextCell(idl);
        }

        auto polydata = vtkSmartPointer<vtkPolyData>::New();
        polydata->SetPoints(points);
        polydata->SetPolys(polys);

        BOOST_LOG_TRIVIAL(debug) << "Done deserialising polydata";

        std::string displayname;
        ar >> displayname;

        ::new(m)Shape(polydata, displayname);
    }
}
}

#endif // SHAPE_H
