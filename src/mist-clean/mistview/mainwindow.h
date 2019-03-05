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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QLabel>
#include <QVTKWidget.h>
#include "transformation.h"
#include "orthoviews.h"
#include "normalview.h"
#include <vtkType.h>
#include <memory>
#include <string>
#include <vector>
#include "plotwindow.h"
#include "newmat.h"

class NoKeysQVTKWidget;
class vtkOrientationMarkerWidget;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void PickSagittal(double x, double y);
    void PickCoronal(double x, double y);
    void PickTransverse(double x, double y);
    void PickNormal(const std::shared_ptr<NormalView>& normalView, double x, double y);

    void SetOrigin(const NEWMAT::ColumnVector& origin);
    void SelectVertex(vtkIdType vertex);
    void SelectVolume(int index);

    static std::vector<double> ReadDisplacements(std::string filename, int vertices);

private slots:
    void on_sagittalSlider_valueChanged(int position);
    void on_coronalSlider_valueChanged(int position);
    void on_transverseSlider_valueChanged(int position);
    void on_vertexIdSpinBox_valueChanged(int value);
    void on_modalityComboBox_currentIndexChanged(int index);
    void on_adjustMmSpinBox_valueChanged(double value);
    void on_xMmSpinBox_valueChanged(double value);
    void on_yMmSpinBox_valueChanged(double value);
    void on_zMmSpinBox_valueChanged(double value);
    void on_windowMinSpinBox_valueChanged(double value);
    void on_windowMaxSpinBox_valueChanged(double value);
    void on_openAction_triggered();
    void on_closeAction_triggered();
    void on_saveAction_triggered();
    void on_loadAction_triggered();
    void on_showProfileFitAction_triggered();
    void on_showMultiShapeAction_triggered();
    void on_saveSlicesPNGAction_triggered();
    void on_setupFilenamesAction_triggered();
    void on_quitAction_triggered();

    void on_slabThicknessSpinBox_valueChanged(double value);

    void on_hideOverlaysCheckBox_stateChanged(int value);

private:
    struct Modality
    {
        std::string name;
        std::string filename;
    };

    Ui::MainWindow *ui;

    const int m_sliderSubDivision = 10;
    const int m_minimumViewWidth = 250;
    const int m_minimumViewHeight = 250;

    std::shared_ptr<OrthoViews> m_orthoViews;
    std::vector<std::shared_ptr<NormalView> > m_normalViews;
    std::vector<std::unique_ptr<NoKeysQVTKWidget> > m_normalViewWidgets;

    boost::shared_ptr<ShapeModel> m_model;

    std::vector<std::shared_ptr<ViewData> > m_volumes;
    std::shared_ptr<ViewShape> m_shape;
    std::shared_ptr<ViewShape> m_extraShape;
    std::shared_ptr<ViewData> m_overlay;
    boost::shared_ptr<const Transformation> m_transformation;

    std::vector<Modality> m_volumeFilenames;
    std::string m_modelFilename;
    std::string m_transformFilename;
    std::string m_displacementsBasename;
    std::string m_normExclusionFilename;
    std::string m_outputFilename;
    std::string m_extraShapeFilename;
    std::string m_overlayFilename;

    QLabel m_statusLabel;
    bool m_unsavedChanges = false;
    std::string m_selectedDirectory;

    PlotWindow m_plotWindow;

    virtual void keyPressEvent(QKeyEvent* event);

    void NormalViewInteractionEventHandler(vtkObject *caller, unsigned long eventId, void *callData);

    void LoadVolume(std::string filename, std::string name);
    void CloseSubject();
    void WriteDisplacements(std::string filename);
    bool SaveInteractive();
    void LoadInteractive();

    void GetFilenames(bool showdialog);

    void SetupInterface();
    void EnableInterface(bool enable);
    void CreateNormalViews();
    void DestroyNormalViews();
};

#endif // MAINWINDOW_H
