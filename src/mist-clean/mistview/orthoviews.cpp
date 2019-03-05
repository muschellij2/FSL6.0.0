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

#include "orthoviews.h"
#include "viewshape.h"
#include "viewdata.h"
#include <newimageall.h>
#include <newmat.h>
#include <vtkPlane.h>
#include <vtkImageData.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkCursor3D.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkParallelCoordinatesInteractorStyle.h>
#include <vtkCamera.h>
#include <array>

using namespace NEWMAT;

OrthoViews::OrthoViews(vtkSmartPointer<vtkRenderWindow> sagittal, vtkSmartPointer<vtkRenderWindow> coronal,
                       vtkSmartPointer<vtkRenderWindow> transverse)
    : m_origin(3)
{
    ColumnVector sagittalnormal(3);
    sagittalnormal << 1.0 << 0.0 << 0.0;
    ColumnVector sagittalviewup(3);
    sagittalviewup << 0.0 << 0.0 << 1.0;
    CreateView(sagittalnormal, sagittalviewup, nullptr, sagittal);

    ColumnVector coronalnormal(3);
    coronalnormal << 0.0 << -1.0 << 0.0;
    ColumnVector coronalviewup(3);
    coronalviewup << 0.0 << 0.0 << 1.0;
    CreateView(coronalnormal, coronalviewup, nullptr, coronal);

    ColumnVector transversenormal(3);
    transversenormal << 0.0 << 0.0 << 1.0;
    ColumnVector transverseviewup(3);
    transverseviewup << 0.0 << 1.0 << 0.0;
    CreateView(transversenormal, transverseviewup, nullptr, transverse);

    ContentsChanged();

    m_cursor = vtkSmartPointer<vtkCursor3D>::New();
    m_cursor->OutlineOff();
    m_cursor->XShadowsOff();
    m_cursor->YShadowsOff();
    m_cursor->ZShadowsOff();
    m_origin = 0.0;
    OriginChanged();
}


void OrthoViews::OriginChanged()
{
    m_cursor->SetFocalPoint(m_origin.Store());

    for (auto it = m_views.begin(); it != m_views.end(); ++it)
        it->plane->SetOrigin(m_origin.Store());

    PlanesChanged();
}

ColumnVector OrthoViews::GetOrigin() const
{
    return m_origin;
}

void OrthoViews::SetOrigin(const ColumnVector& origin)
{
    m_origin = origin;
    OriginChanged();
}

ColumnVector OrthoViews::GetBounds() const
{
    ColumnVector bounds(6);
    if (m_viewData)
        bounds << m_viewData->GetImageData()->GetBounds();
    else
        bounds = 0.0;

    return bounds;
}

vtkSmartPointer<vtkRenderWindow> OrthoViews::GetSagittalView() const
{
    return m_views[0].renderWindow;
}

vtkSmartPointer<vtkRenderWindow> OrthoViews::GetCoronalView() const
{
    return m_views[1].renderWindow;
}

vtkSmartPointer<vtkRenderWindow> OrthoViews::GetTransverseView() const
{
    return m_views[2].renderWindow;
}

vtkIdType OrthoViews::PickSagittal(const double& x, const double& y)
{
    return ViewBase::Pick(m_views[0], x, y);
}

vtkIdType OrthoViews::PickCoronal(const double& x, const double& y)
{
    return ViewBase::Pick(m_views[1], x, y);
}

vtkIdType OrthoViews::PickTransverse(const double& x, const double& y)
{
    return ViewBase::Pick(m_views[2], x, y);
}

void OrthoViews::ContentsChanged()
{
    ViewBase::ContentsChanged();

    if (m_viewData)
        m_cursor->SetModelBounds(m_viewData->GetImageData()->GetBounds());
}

void OrthoViews::MarkerPointIDChanged()
{
    if (m_markerShape)
    {
        ColumnVector position(3);
        position << m_markerShape->GetPolyData()->GetPoint(m_markerPointId);
        m_origin = position;
        OriginChanged();
    }
}

void OrthoViews::SetupSceneForView(View& vw)
{
    ViewBase::SetupSceneForView(vw);

    auto cursormapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cursormapper->SetInputConnection(m_cursor->GetOutputPort());
    auto cursoractor = vtkSmartPointer<vtkActor>::New();
    cursoractor->SetMapper(cursormapper);
    cursoractor->GetProperty()->SetColor(0.0, 1.0, 0.0);
    auto cursorrenderer = vtkSmartPointer<vtkRenderer>::New();
    cursorrenderer->SetLayer(vw.renderWindow->GetNumberOfLayers());
    cursorrenderer->AddActor(cursoractor);
    vw.renderWindow->SetNumberOfLayers(vw.renderWindow->GetNumberOfLayers() + 1);
    vw.renderWindow->AddRenderer(cursorrenderer);
}
