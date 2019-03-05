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

#include "shape.h"
#include <vtkCellArray.h>
#include <vtkImageData.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSelectEnclosedPoints.h>
#include <iostream>
#include <sstream>
#include <set>
#include <functional>
#include "boost/make_shared.hpp"

using namespace NEWMAT;

Shape::Shape(vtkSmartPointer<vtkPolyData> polydata, const std::string& displayname)
    : m_linksBuilt(false),
      m_displayName(displayname),
      m_transformation(boost::make_shared<IdentityTransformation>())
{
    Init(polydata);
}

Shape::Shape(const std::string& filename, const std::string& displayname)
    : m_linksBuilt(false),
      m_displayName(displayname),
      m_transformation(boost::make_shared<IdentityTransformation>())
{
    Init(ReadPolyData(filename));
}

Shape::Shape(const std::string& filename, boost::shared_ptr<const Transformation> &transformation,
             const std::string& displayname)
    : m_linksBuilt(false),
      m_displayName(displayname)
{
    m_transformation = transformation;
    Init(ReadPolyData(filename));
}

Shape::Shape(const Shape &other)
    : m_linksBuilt(false)
{
    Copy(other);
}

Shape::Shape(const Shape &other, boost::shared_ptr<const Transformation> &transformation)
    : m_linksBuilt(false)
{
    auto points = vtkSmartPointer<vtkPoints>::New();
    const auto &otherup = other.m_untransformedPoints;
    points->SetNumberOfPoints(otherup.size());
    for (vtkIdType i = 0; i < otherup.size(); i++)
        points->SetPoint(i, otherup[i](1), otherup[i](2), otherup[i](3));

    auto polys = vtkSmartPointer<vtkCellArray>::New();
    polys->DeepCopy(other.m_polyData->GetPolys());

    auto newpd = vtkSmartPointer<vtkPolyData>::New();
    newpd->SetPoints(points);
    newpd->SetPolys(polys);

    m_transformation = transformation;
    Init(newpd);

    for (int i = 0; i < GetNumberOfVertices(); i++)
        SetDisplacement(i, other.m_displacements[i]);

    m_displayName = other.m_displayName;
}

Shape& Shape::operator=(const Shape &other)
{
    // Not sure if VTK's deep copy will accept copying from self
    if (&other != this)
        Copy(other);

    m_linksBuilt = false;

    return *this;
}

void Shape::Copy(const Shape &other)
{
    m_displacements = other.m_displacements;
    m_displayName = other.m_displayName;
    m_polyData = vtkSmartPointer<vtkPolyData>::New();
    m_polyData->DeepCopy(other.m_polyData);
    m_transformation = other.m_transformation;
    m_untransformedPoints = other.m_untransformedPoints;
    m_untransformedNormals = other.m_untransformedNormals;
}

void Shape::Init(vtkSmartPointer<vtkPolyData> polydata)
{
    m_untransformedPoints.clear();
    m_untransformedNormals.clear();

    // The original normals are needed when setting deformations ..
    vtkSmartPointer<vtkPolyData> orignormals = ComputeInitialNormals(polydata);
    for (int i = 0; i < orignormals->GetNumberOfPoints(); i++)
    {
        ColumnVector n(3);
        orignormals->GetPointData()->GetNormals()->GetTuple(i, n.Store());
        m_untransformedNormals.push_back(n);
    }

    for (int i = 0; i < polydata->GetNumberOfPoints(); i++)
    {
        ColumnVector p(3);
        polydata->GetPoints()->GetPoint(i, p.Store());
        m_untransformedPoints.push_back(p);
        ColumnVector tp = m_transformation->InverseTransformPoint(p);
        polydata->GetPoints()->SetPoint(i, tp.Store());
    }

    m_polyData = ComputeInitialNormals(polydata);

    for (int i = 0; i < m_polyData->GetNumberOfPoints(); i++)
        m_displacements.push_back(0.0);
}

vtkSmartPointer<vtkPolyData> Shape::ReadPolyData(const std::string &filename)
{
    BOOST_LOG_TRIVIAL(debug) << "Loading " << filename;

    auto reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    if (int errorcode = reader->GetErrorCode())
    {
        std::ostringstream what;
        what << "vtkPolyDataReader() failed with error code " << errorcode << ".";
        throw ShapeException(what.str());
    }

    // We just want the points and polys, no scalars/normals/etc that might be in the file
    vtkSmartPointer<vtkPolyData> origpd = reader->GetOutput();
    auto points = vtkSmartPointer<vtkPoints>::New();
    auto polys = vtkSmartPointer<vtkCellArray>::New();
    auto newpd = vtkSmartPointer<vtkPolyData>::New();
    points->DeepCopy(origpd->GetPoints());
    polys->DeepCopy(origpd->GetPolys());
    newpd->SetPoints(points);
    newpd->SetPolys(polys);

    return newpd;
}

void Shape::WritePolyData(const string &filename) const
{
    BOOST_LOG_TRIVIAL(debug) << "Writing polydata to " << filename;

    // ASCII output is the default
    auto writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(m_polyData);
    if (!writer->Write())
        throw ShapeException(std::string("Writing VTK file ") + filename + " failed");
}

vtkSmartPointer<vtkPolyData> Shape::ComputeInitialNormals(vtkSmartPointer<vtkPolyData> polydata)
{
    // These are called initial normals because they are not updated after calling
    // SetDeformation(); this is the intended behaviour.

    auto normalfilter = vtkSmartPointer<vtkPolyDataNormals>::New();
    normalfilter->SetInputData(polydata); // VTK6
    // normalfilter->SetInput(polydata); // VTK5
    normalfilter->ComputeCellNormalsOff();
    normalfilter->ComputePointNormalsOn();
    normalfilter->SplittingOff();
    normalfilter->ConsistencyOff(); // Polys may still be reordered as we are using AutoOrient (points will not change)
    normalfilter->AutoOrientNormalsOn();
    normalfilter->Update();

    return normalfilter->GetOutput();
}

vtkSmartPointer<vtkPolyData> Shape::GetPolyData() const
{
    return m_polyData;
}

ColumnVector Shape::GetPoint(vtkIdType vertex) const
{
    ColumnVector point(3);
    point << m_polyData->GetPoint(vertex);
    return point;
}

ColumnVector Shape::GetNormal(vtkIdType vertex) const
{
    ColumnVector normal(3);
    normal << m_polyData->GetPointData()->GetNormals()->GetTuple3(vertex);
    return normal;
}

ColumnVector Shape::GetBounds() const
{
    ColumnVector bounds(6);
    bounds << m_polyData->GetBounds();
    return bounds;
}

ColumnVector Shape::GetCenter() const
{
    ColumnVector center(3);
    center << m_polyData->GetCenter();
    return center;
}

vtkIdType Shape::GetNumberOfVertices() const
{
    return m_polyData->GetNumberOfPoints();
}

std::string Shape::GetDisplayName() const
{
    return m_displayName;
}

Matrix Shape::GetDistanceMatrix() const
{
    std::vector<std::set<vtkIdType> > neighbours(GetNumberOfVertices());

    m_polyData->GetPolys()->InitTraversal();
    auto idl = vtkSmartPointer<vtkIdList>::New();
    while (m_polyData->GetPolys()->GetNextCell(idl))
    {
        std::set<vtkIdType> poly;
        for (vtkIdType p = 0; p < idl->GetNumberOfIds(); p++)
            poly.insert(idl->GetId(p));

        for (auto s = poly.cbegin(); s != poly.cend(); s++)
            for (auto t = poly.cbegin(); t != poly.cend(); t++)
                if (s != t)
                    neighbours[*s].insert(*t);
    }

    Matrix D(GetNumberOfVertices(), GetNumberOfVertices());
    D = 0.0;

    for (vtkIdType i = 0; i < GetNumberOfVertices(); i++)
    {
        std::function<void(std::set<vtkIdType>, int)> recurse
                = [&](std::set<vtkIdType> points, int depth)
                {
                    std::set<vtkIdType> newpoints;

                    for (vtkIdType s : points)
                        for (vtkIdType p : neighbours[s])
                            if (p != i && D(i + 1, p + 1) == 0.0)
                                newpoints.insert(p);

                    for (auto p : newpoints)
                        D(i + 1, p + 1) = depth;

                    if (newpoints.size())
                        recurse(newpoints, ++depth);
                };

        recurse(std::set<vtkIdType>({ i }), 1);
    }

    return D;
}

std::vector<std::vector<vtkIdType> > Shape::GetPolysForVertex(vtkIdType vertex)
{
    // GetPointCells calls this by itself if necessary - although the documentations says it doesn't?
    if (!m_linksBuilt)
    {
        m_polyData->BuildLinks();
        m_linksBuilt = true;
    }

    std::vector<std::vector<vtkIdType> > polys;

    vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
    m_polyData->GetPointCells(vertex, cells);

    for (vtkIdType i = 0; i < cells->GetNumberOfIds(); i++)
    {
        vtkIdType npoints;
        vtkIdType *points;

        m_polyData->GetCellPoints(cells->GetId(i), npoints, points);

        polys.push_back(std::vector<vtkIdType>(points, points + npoints));
    }

    return polys;
}

void Shape::ToVolume(NEWIMAGE::volume<char> &vol, bool usevoxelvertices)
{
    vol = 0;
    vol.setDisplayMaximumMinimum(1.0, 0.0);

    if (usevoxelvertices)
    {
        auto imd = vtkImageData::New();
        imd->SetExtent(0, vol.xsize(), 0, vol.ysize(), 0, vol.zsize());
        imd->SetSpacing(vol.xdim(), vol.ydim(), vol.zdim());
        imd->SetOrigin(-0.5 * vol.xdim(), -0.5 * vol.ydim(), -0.5 * vol.zdim());

        auto filt = vtkSelectEnclosedPoints::New();
        filt->SetInputData(imd);
        filt->SetSurfaceData(m_polyData);
        filt->SetTolerance(0.00001);
        filt->Update();

        vtkSmartPointer<vtkDataArray> sel = filt->GetOutput()->GetPointData()->GetArray("SelectedPoints");

        int nextx = 1;
        int nexty = vol.xsize() + 1;
        int nextz = (vol.ysize() + 1) * (vol.xsize() + 1);

        vtkIdType id = 0;
        for (int k = 0; k < vol.zsize(); k++)
        {
            for (int j = 0; j < vol.ysize(); j++)
            {
                for (int i = 0; i < vol.xsize(); i++)
                {
                    vol(i, j, k) = sel->GetTuple1(id) + sel->GetTuple1(id + nextx)
                                 + sel->GetTuple1(id + nexty) + sel->GetTuple1(id + nexty + nextx)
                                 + sel->GetTuple1(id + nextz) + sel->GetTuple1(id + nextz + nextx)
                                 + sel->GetTuple1(id + nextz + nexty) + sel->GetTuple1(id + nextz + nexty + nextx) > 0.0;
                    id++;
                }
                id += nextx;
            }
            id += nexty;
        }
        id += nextz;
    }
    else
    {
        auto imd = vtkImageData::New();
        imd->SetExtent(0, vol.xsize() - 1, 0, vol.ysize() - 1, 0, vol.zsize() - 1);
        imd->SetSpacing(vol.xdim(), vol.ydim(), vol.zdim());

        auto filt = vtkSelectEnclosedPoints::New();
        filt->SetInputData(imd);
        filt->SetSurfaceData(m_polyData);
        filt->SetTolerance(0.00001);
        filt->Update();

        vtkSmartPointer<vtkDataArray> sel = filt->GetOutput()->GetPointData()->GetArray("SelectedPoints");

        vtkIdType id = 0;
        for (int k = 0; k < vol.zsize(); k++)
            for (int j = 0; j < vol.ysize(); j++)
                for (int i = 0; i < vol.xsize(); i++)
                    vol(i, j, k) = sel->GetTuple1(id++) > 0.0;
    }
}

void Shape::SetDisplacement(vtkIdType vertex, double amount)
{
    m_linksBuilt = false;

    m_displacements[vertex] = amount;

    ColumnVector newpoint = m_transformation->InverseTransformPoint(
            m_untransformedPoints[vertex] + amount * m_untransformedNormals[vertex]);

    m_polyData->GetPoints()->SetPoint(vertex, newpoint.Store());
    m_polyData->GetPoints()->Modified();
}

double Shape::GetDisplacement(vtkIdType vertex) const
{
    return m_displacements[vertex];
}

Shape Shape::RemoveSelfIntersections(double stepsize) const
{
    Shape result(*this);

    for (int vert = 0; vert < result.GetNumberOfVertices(); vert++)
        result.SetDisplacement(vert, 0.0);

    double maxdisp = std::abs(*std::max_element(
                        m_displacements.cbegin(), m_displacements.cend(),
                        [] (double a, double b) { return std::abs(a) < std::abs(b); }));
    int iters = static_cast<int>(std::ceil(maxdisp / stepsize));
    auto filt = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
    filt->SetTolerance(0.00001);

    // NB: This only updates the mesh after each full iteration - should be fine with small step size
    for (int iter = 0; iter < iters; iter++)
    {
        std::vector<double> newdisplacements(result.GetNumberOfVertices(), 0.0);

        filt->Initialize(result.m_polyData);

        for (int vert = 0; vert < result.GetNumberOfVertices(); vert++)
        {
            double fitted = m_displacements[vert];
            double current = result.m_displacements[vert];
            double newdisp;
            bool testoutside;

            if (fitted > 0.0)
            {
                newdisp = std::min(fitted, current + stepsize);
                testoutside = true;
            }
            else if (fitted < 0.0)
            {
                newdisp = std::max(fitted, current - stepsize);
                testoutside = false;
            }
            else
                continue;

            ColumnVector newpoint = result.m_transformation->InverseTransformPoint(
                        result.m_untransformedPoints[vert] + newdisp * m_untransformedNormals[vert]);

            if (testoutside == (filt->IsInsideSurface(newpoint(1), newpoint(2), newpoint(3)) == 0))
                newdisplacements[vert] = newdisp;
            else
                newdisplacements[vert] = current;
        }

        filt->Complete();

        for (int vert = 0; vert < result.GetNumberOfVertices(); vert++)
            result.SetDisplacement(vert, newdisplacements[vert]);
    }

    return result;
}
