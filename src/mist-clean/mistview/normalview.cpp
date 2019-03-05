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

#include "normalview.h"
#include "viewshape.h"
#include "viewdata.h"
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkRenderWindow.h>
#include <vtkPlane.h>
#include <vtkLineSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkRendererCollection.h>
#include <vtkCamera.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkAxesActor.h>

using namespace NEWMAT;

NormalView::NormalView(vtkSmartPointer<vtkCamera> camera, vtkSmartPointer<vtkRenderWindow> renderwindow)
{
    ColumnVector initnormal(3);
    initnormal << 1.0 << 0.0 << 0.0;
    ColumnVector initviewup(3);
    initviewup << 0.0 << 0.0 << 1.0;
    CreateView(initnormal, initviewup, camera, renderwindow);
}

vtkSmartPointer<vtkRenderWindow> NormalView::GetView()
{
    return m_views[0].renderWindow;
}

vtkIdType NormalView::Pick(const double &x, const double &y)
{
    return ViewBase::Pick(m_views[0], x, y);
}

//void NormalView::SetupOrientationWidget(vtkSmartPointer<vtkRenderWindowInteractor> interactor)
//{
//    auto axes = vtkSmartPointer<vtkAxesActor>::New();
//    axes->SetPosition(100, 100, 100);
//    m_orientationWidget = vtkSmartPointer<vtkOrientationMarkerWidgetFixed>::New();
//    m_orientationWidget->SetOrientationMarker(axes);
//    m_orientationWidget->SetInteractor(interactor);
//}

void NormalView::ContentsChanged()
{
    ViewBase::ContentsChanged();

    if (m_markerShape)
    {
        ColumnVector point = m_markerShape->GetPoint(m_markerPointId);
        ColumnVector normal = m_markerShape->GetNormal(m_markerPointId);
        ColumnVector center = m_markerShape->GetCenter();
        ColumnVector delta = point - center;
        delta /= std::sqrt(delta(1) * delta(1) + delta(2) * delta(2) + delta(3) * delta(3));
        ColumnVector planenormal(3);
        planenormal << normal(2) * delta(3) - normal(3) * delta(2)
                    << normal(3) * delta(1) - normal(1) * delta(3)
                    << normal(1) * delta(2) - normal(2) * delta(1);

        m_views[0].plane->SetOrigin(point.Store());
        m_views[0].plane->SetNormal(planenormal.Store());
        m_views[0].camera->SetViewUp(normal.Store());

        PlanesChanged();
    }
}

void NormalView::SetupSceneForView(View& vw)
{
    ViewBase::SetupSceneForView(vw);

    if (m_markerShape)
    {
        ColumnVector point =  m_markerShape->GetPoint(m_markerPointId);
        ColumnVector normal = m_markerShape->GetNormal(m_markerPointId);
        ColumnVector point1 = point - 5 * normal;
        ColumnVector point2 = point + 5 * normal;

        auto normalsource = vtkSmartPointer<vtkLineSource>::New();
        normalsource->SetPoint1(point1.Store());
        normalsource->SetPoint2(point2.Store());
        auto normalmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        normalmapper->SetInputConnection(normalsource->GetOutputPort());
        auto normalactor = vtkSmartPointer<vtkActor>::New();
        normalactor->SetMapper(normalmapper);
        normalactor->GetProperty()->SetColor(1.0, 0.0, 0.0);
        auto normalrenderer = vtkSmartPointer<vtkRenderer>::New();
        normalrenderer->SetLayer(vw.renderWindow->GetNumberOfLayers());
        normalrenderer->AddActor(normalactor);
        vw.renderWindow->SetNumberOfLayers(vw.renderWindow->GetNumberOfLayers() + 1);
        vw.renderWindow->AddRenderer(normalrenderer);
    }
}

void NormalView::SetupCameraForView(View& vw)
{
    ColumnVector center = m_markerShape->GetCenter();
    ColumnVector normal(3);
    normal << vw.plane->GetNormal();
    ColumnVector position = center + normal;
    vw.camera->SetFocalPoint(center.Store());
    vw.camera->SetPosition(position.Store());
    vw.camera->ParallelProjectionOn();
    ColumnVector bounds(6);
    bounds << m_markerShape->GetPolyData()->GetBounds();
    vw.camera->SetParallelScale(std::max(std::max(bounds(2) - bounds(1), bounds(4) - bounds(3)), bounds(6) - bounds(5)));
    vw.camera->Zoom(m_scale);
}

//void NormalView::SetupMarkersForView(View &vw)
//{
//    ViewBase::SetupMarkersForView(vw);

//    if (m_orientationWidget)
//    {
////        m_orientationWidget->SetEnabled(0);
//        m_orientationWidget->SetEnabled(1);
//        m_orientationWidget->InteractiveOff();
//    }
//}
