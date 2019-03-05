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

#ifndef VIEWBASE_H
#define VIEWBASE_H

#include <memory>
#include <vector>
#include <vtkSmartPointer.h>
#include <vtkType.h>

#include <newmat.h>

class vtkLookupTable;
class vtkPlane;
class vtkPolyData;
class vtkRenderer;
class vtkRenderWindow;
class vtkCamera;
class ViewShape;
class ViewData;

class ViewBase
{
public:
    enum HideMode
    {
        HideNone,
        HideMarkers,
        HideAllOverlays
    };

    ViewBase();
    void ResetDataAndShapes();
    void AddShape(std::shared_ptr<ViewShape> shape);
    void SetData(std::shared_ptr<ViewData> data);
    void SetOverlay(std::shared_ptr<ViewData> overlay);
    void SetMarker(std::shared_ptr<ViewShape> shape, vtkIdType pointid);
    void SetMarkerTolerance(const double& value);
    double GetMarkerTolerance() const;
    void SetHideOverlays(HideMode value);
    HideMode GetHideOverlays();
    void SetScale(const double& scale);
    void ShowVertices(bool show);
    // TODO: Handle this properly (needed to handle changes in shapes)
    void RedrawMarkers();
    void Render();
    void ResetCameras();

protected:
    struct View
    {
        vtkSmartPointer<vtkPlane> plane;
        vtkSmartPointer<vtkRenderWindow> renderWindow;
        vtkSmartPointer<vtkCamera> camera;
        vtkSmartPointer<vtkPolyData> visibleVerts;
        vtkSmartPointer<vtkRenderer> markerRenderer;
    };

    std::vector<View> m_views;
    std::shared_ptr<ViewData> m_viewData;
    std::shared_ptr<ViewData> m_overlay;
    std::vector<std::shared_ptr<ViewShape> > m_shapes;

    double m_scale = 1.0;
    HideMode m_hideOverlays = HideNone;

    bool m_showVertices = true;
    std::shared_ptr<ViewShape> m_markerShape = nullptr;
    vtkIdType m_markerPointId = 0;
    double m_markerTolerance = 0.5;

    void CreateView(const NEWMAT::ColumnVector& normal, const NEWMAT::ColumnVector& viewup,
                    vtkSmartPointer<vtkCamera> camera, vtkSmartPointer<vtkRenderWindow> renderwindow);
    virtual void ContentsChanged();
    virtual void PlanesChanged();
    virtual void MarkerPointIDChanged();
    virtual void SetupSceneForView(View& vw);
    virtual void SetupCameraForView(View& vw);
    virtual void SetupMarkersForView(View& vw);
    void SetMarkerVisibilityForView(View& vw);
    vtkIdType Pick(View& vw, const double& x, const double& y);

private:
    ViewBase(const ViewBase&);
    ViewBase& operator=(const ViewBase&);
};

#endif // VIEWBASE_H
