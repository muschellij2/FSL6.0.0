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

#include "viewbase.h"
#include "viewshape.h"
#include "viewdata.h"
#include <cmath>
#include <newimageall.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPlane.h>
#include <vtkImageData.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkImageProperty.h>
#include <vtkLookupTable.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkCamera.h>
#include <vtkCutter.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkSphereSource.h>
#include <vtkGlyph3D.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkPointPicker.h>
#include <vtkImageInterpolator.h>

// NB: Need to call ResetCameras() after setting up or this will not show anything!

using namespace NEWMAT;

ViewBase::ViewBase()
{

}

void ViewBase::ResetDataAndShapes()
{
    for (auto it = m_views.begin(); it != m_views.end(); ++it)
    {
        it->renderWindow->GetRenderers()->RemoveAllItems();
        it->renderWindow->AddRenderer(vtkSmartPointer<vtkRenderer>::New());
        it->visibleVerts = nullptr;
        it->markerRenderer = nullptr;
    }

    m_viewData.reset();
    m_shapes.clear();
    m_markerShape.reset();
    m_markerPointId = 0;
}

void ViewBase::SetData(std::shared_ptr<ViewData> data)
{
    m_viewData = data;
    ContentsChanged();

    for (auto it = m_views.begin(); it != m_views.end(); ++it)
        SetMarkerVisibilityForView(*it);

    Render();
}

void ViewBase::SetOverlay(std::shared_ptr<ViewData> overlay)
{
    m_overlay = overlay;
    ContentsChanged();

    for (auto it = m_views.begin(); it != m_views.end(); ++it)
        SetMarkerVisibilityForView(*it);

    Render();
}

void ViewBase::AddShape(std::shared_ptr<ViewShape> shape)
{
    m_shapes.push_back(shape);
    ContentsChanged();
}

void ViewBase::SetMarker(std::shared_ptr<ViewShape> shape, vtkIdType pointid)
{
    m_markerShape = shape;
    m_markerPointId = pointid;

    ContentsChanged();
    MarkerPointIDChanged();
}

void ViewBase::SetMarkerTolerance(const double& value)
{
    m_markerTolerance = value;
    RedrawMarkers();
}

double ViewBase::GetMarkerTolerance() const
{
    return m_markerTolerance;
}

void ViewBase::SetScale(const double &scale)
{
    m_scale = scale;
}

void ViewBase::SetHideOverlays(HideMode mode)
{
    m_hideOverlays = mode;

    for (auto it = m_views.begin(); it != m_views.end(); ++it)
    {
        vtkSmartPointer<vtkRendererCollection> renderers = it->renderWindow->GetRenderers();
        // Layer 0 is background (never disabled), layer 1 contains shapes ...
        static_cast<vtkRenderer *>(renderers->GetItemAsObject(1))->SetDraw(m_hideOverlays != HideAllOverlays);

        // Other layers
        for (int i = 2; i < renderers->GetNumberOfItems(); i++)
            static_cast<vtkRenderer *>(renderers->GetItemAsObject(i))->SetDraw(m_hideOverlays == HideNone);
    }

    Render();
}

ViewBase::HideMode ViewBase::GetHideOverlays()
{
    return m_hideOverlays;
}

void ViewBase::ShowVertices(bool show)
{
    m_showVertices = show;
}

void ViewBase::Render()
{
    for (auto it = m_views.begin(); it != m_views.end(); ++it)
        if (it->renderWindow->IsDrawable()) // VTK6
           it->renderWindow->Render();
}

void ViewBase::ResetCameras()
{
    if (m_viewData)
    {
        for (auto it = m_views.begin(); it != m_views.end(); ++it)
            SetupCameraForView(*it);

        Render();
    }
}

void ViewBase::ContentsChanged()
{
    if (m_viewData)
    {
        for (auto it = m_views.begin(); it != m_views.end(); ++it)
        {
            SetupSceneForView(*it);
            SetupMarkersForView(*it);

            auto renderercollection = it->renderWindow->GetRenderers();
            renderercollection->InitTraversal();
            while(auto renderer = renderercollection->GetNextItem())
                renderer->SetActiveCamera(it->camera);

        }

        Render();
    }
}

void ViewBase::PlanesChanged()
{
    if (m_viewData)
        for (auto it = m_views.begin(); it != m_views.end(); ++it)
            SetMarkerVisibilityForView(*it);

    Render();
}

void ViewBase::MarkerPointIDChanged()
{
    for (auto it = m_views.begin(); it != m_views.end(); ++it)
        SetMarkerVisibilityForView(*it);
}

void ViewBase::RedrawMarkers()
{
    // Full redraw is workaround for VTK 5.10 - setting modified on (points in) polydata
    // doesn't seem to cause the cutter to update ..
    // ContentsChanged(); // VTK5

    PlanesChanged();
}

void ViewBase::CreateView(const ColumnVector& normal, const ColumnVector& viewup,
                          vtkSmartPointer<vtkCamera> camera, vtkSmartPointer<vtkRenderWindow> renderwindow)
{
    // Dummy renderer is added to avoid screen corruption until we have a ViewData
    // object and can properly initialise the vtkRenderWindows

    auto plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetNormal(normal.Store());
    if (camera == nullptr)
        camera = vtkSmartPointer<vtkCamera>::New();
    camera->SetViewUp(viewup.Store());
    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetActiveCamera(camera);
    renderwindow->AddRenderer(renderer);

    View newview = {plane, renderwindow, camera, nullptr, nullptr};
    m_views.push_back(newview);
}

void ViewBase::SetupSceneForView(View& vw)
{
    vw.renderWindow->GetRenderers()->RemoveAllItems();

    auto imagemapper = vtkSmartPointer<vtkImageResliceMapper>::New();
    imagemapper->SetInputData(m_viewData->GetImageData()); // VTK6
    // imagemapper->SetInput(m_viewData->GetImageData()); // VTK5
    imagemapper->SetSlicePlane(vw.plane);
    auto interpolator = vtkSmartPointer<vtkImageInterpolator>::New();
    interpolator->SetInterpolationModeToCubic();
    imagemapper->SetInterpolator(interpolator);
    auto imageslice = vtkSmartPointer<vtkImageSlice>::New();
    imageslice->SetMapper(imagemapper);
    imageslice->GetProperty()->SetLookupTable(m_viewData->GetLut());
    imageslice->GetProperty()->UseLookupTableScalarRangeOn();
    // imageslice->GetProperty()->SetInterpolationTypeToNearest();
    auto imagerenderer = vtkSmartPointer<vtkRenderer>::New();
    imagerenderer->AddActor(imageslice);
    imagerenderer->SetLayer(0);
    vw.renderWindow->AddRenderer(imagerenderer);

    auto shaperenderer = vtkSmartPointer<vtkRenderer>::New();
    shaperenderer->SetLayer(1);
    vw.renderWindow->AddRenderer(shaperenderer);

    // Same layer as shape
    if (m_overlay)
    {
        auto overlaymapper = vtkSmartPointer<vtkImageResliceMapper>::New();
        overlaymapper->SetInputData(m_overlay->GetImageData());
        overlaymapper->SetSlicePlane(vw.plane);
        auto overlayslice = vtkSmartPointer<vtkImageSlice>::New();
        overlayslice->SetMapper(overlaymapper);
        overlayslice->GetProperty()->SetLookupTable(m_overlay->GetLut());
        overlayslice->GetProperty()->UseLookupTableScalarRangeOn();
        overlayslice->GetProperty()->SetInterpolationTypeToNearest();
        overlayslice->GetProperty()->SetOpacity(0.5);
        shaperenderer->AddActor(overlayslice);
    }

    for (auto it = m_shapes.begin(); it != m_shapes.end(); ++it)
    {
        auto shapecutter = vtkSmartPointer<vtkCutter>::New();
        shapecutter->SetCutFunction(vw.plane);
        shapecutter->SetInputData((*it)->GetPolyData()); // VTK6
//        shapecutter->SetInput((*it)->GetPolyData()); // VTK5
        auto shapemapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        shapemapper->SetInputConnection(shapecutter->GetOutputPort());
        auto shapeactor = vtkSmartPointer<vtkActor>::New();
        shapeactor->SetMapper(shapemapper);
        shapeactor->GetProperty()->LightingOff();
        shapeactor->GetProperty()->SetColor((*it)->GetColor().data());
        shapeactor->GetProperty()->SetLineWidth((*it)->GetLineWidth());
        shapeactor->GetProperty()->SetOpacity((*it)->GetOpacity());
        shaperenderer->AddActor(shapeactor);
    }

    vw.renderWindow->SetNumberOfLayers(2);
}

void ViewBase::SetupMarkersForView(View& vw)
{
    if (m_markerShape)
    {
        vw.visibleVerts = vtkSmartPointer<vtkPolyData>::New();

        auto glyphsource = vtkSmartPointer<vtkSphereSource>::New();
        glyphsource->SetRadius(0.3);
        auto vertfilter = vtkSmartPointer<vtkGlyph3D>::New();
        vertfilter->SetInputData(vw.visibleVerts); // VTK6
//        vertfilter->SetInput(vw.visibleVerts); // VTK5
        vertfilter->SetSourceConnection(glyphsource->GetOutputPort());
        vertfilter->SetColorModeToColorByScalar();
        vertfilter->SetRange(0.0, 1.0);
        vertfilter->ScalingOff();
        vertfilter->OrientOff();
        auto vertmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        vertmapper->SetInputConnection(vertfilter->GetOutputPort());
        auto vertactor = vtkSmartPointer<vtkActor>::New();
        vertactor->GetProperty()->SetColor(0.0, 0.0, 1.0);
        vertactor->SetMapper(vertmapper);
        vw.markerRenderer = vtkSmartPointer<vtkRenderer>::New();
        vw.markerRenderer->SetLayer(vw.renderWindow->GetNumberOfLayers());
        vw.markerRenderer->AddActor(vertactor);
        vw.renderWindow->SetNumberOfLayers(vw.renderWindow->GetNumberOfLayers() + 1);
        vw.renderWindow->AddRenderer(vw.markerRenderer);
    }
    else
    {
        vw.visibleVerts = nullptr;
        vw.markerRenderer = nullptr;
    }
}

void ViewBase::SetMarkerVisibilityForView(View& vw)
{
    if (m_markerShape && vw.visibleVerts)
    {
        auto visiblepoints = vtkSmartPointer<vtkPoints>::New();
        auto glyphcolors = vtkSmartPointer<vtkDoubleArray>::New();
        visiblepoints->Allocate(m_markerShape->GetNumberOfVertices());
        glyphcolors->Allocate(m_markerShape->GetNumberOfVertices());

        for (vtkIdType i = 0; i < m_markerShape->GetNumberOfVertices(); i++)
        {
            double point[3];
            m_markerShape->GetPolyData()->GetPoints()->GetPoint(i, point);
            double dist = vw.plane->EvaluateFunction(point);
            if (std::abs(dist) <= m_markerTolerance &&
                    (m_showVertices || (i == m_markerPointId)))
            {
                visiblepoints->InsertNextPoint(point);

                if (i == m_markerPointId)
                    glyphcolors->InsertNextTuple1(0.0); // Selected vertex - red
                else if (dist < 0)
                    glyphcolors->InsertNextTuple1(1.0); // Below plane - blue
                else
                    glyphcolors->InsertNextTuple1(0.5); // Above plane - green
            }
        }

        vw.visibleVerts->SetPoints(visiblepoints);
        vw.visibleVerts->GetPointData()->SetScalars(glyphcolors);
        vw.visibleVerts->Modified();
    }
}

void ViewBase::SetupCameraForView(View& vw)
{
    ColumnVector center(3);
    center << m_viewData->GetImageData()->GetCenter();
    ColumnVector normal(3);
    normal << vw.plane->GetNormal();
    ColumnVector position = center + normal;
    vw.camera->SetFocalPoint(center.Store());
    vw.camera->SetPosition(position.Store());
    vw.camera->ParallelProjectionOn();
    vw.renderWindow->GetRenderers()->GetFirstRenderer()->ResetCamera();
    vw.camera->Zoom(m_scale);
}

vtkIdType ViewBase::Pick(View& vw, const double& x, const double& y)
{
    if (m_markerShape)
    {
        auto picker = vtkSmartPointer<vtkPointPicker>::New();
        // TODO: Store actor in View? This seems a bit fragile ..
        picker->AddPickList(vw.markerRenderer->GetActors()->GetLastActor());
        picker->PickFromListOn();
        picker->SetTolerance(0.01);

        if (picker->Pick(x, y, 0.0, vw.markerRenderer))
        {
            vtkIdType point = m_markerShape->GetPolyData()->FindPoint(picker->GetPickPosition());
            return point;
        }
    }

    return -1;
}
