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

#include "multishapewindow.h"
#include "ui_multishapewindow.h"
#include "setupfilenameswindow.h"
#include "orthoviews.h"
#include "custominteractorstyle.h"
#include "nokeysqvtkwidget.h"
#include "mainwindow.h"
#include "viewshape.h"
#include "viewdata.h"
#include "vtkInteractorStyleSwitch.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindow.h"
#include "vtkSmartPointer.h"
#include <QMessageBox>

#define BLOCK_SIGNALS(object, call) { bool oldstate = object->blockSignals(true); object->call; object->blockSignals(oldstate); }

MultiShapeWindow::MultiShapeWindow(QWidget *parent, const std::vector<std::shared_ptr<ViewData> > &volumes,
                                   boost::shared_ptr<const Transformation> transformation, const std::string &directory) :
    QDialog(parent),
    ui(new Ui::MultiShapeWindow),
    m_volumes(volumes)
{
    ui->setupUi(this);

    m_sagittalWidget = new NoKeysQVTKWidget(this);
    m_sagittalWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    m_sagittalWidget->setMinimumSize(m_minimumViewWidth, m_minimumViewHeight);
    ui->sagittalLayout->insertWidget(0, m_sagittalWidget);

    m_coronalWidget = new NoKeysQVTKWidget(this);
    m_coronalWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    m_coronalWidget->setMinimumSize(m_minimumViewWidth, m_minimumViewHeight);
    ui->coronalLayout->insertWidget(0, m_coronalWidget);

    m_transverseWidget = new NoKeysQVTKWidget(this);
    m_transverseWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    m_transverseWidget->setMinimumSize(m_minimumViewWidth, m_minimumViewHeight);
    ui->transverseLayout->insertWidget(0, m_transverseWidget);

    m_orthoViews = make_shared<OrthoViews>(m_sagittalWidget->GetRenderWindow(),
                                           m_coronalWidget->GetRenderWindow(),
                                           m_transverseWidget->GetRenderWindow());
    m_orthoViews->SetScale(2.0);

    auto originfunc = [this] (const ColumnVector& origin)
                {
                    SetOrigin(origin);
                };

    auto sagittalinteractorstyle = vtkSmartPointer<CustomInteractorStyle>::New();
    sagittalinteractorstyle->SetOrthoViews(m_orthoViews);
    sagittalinteractorstyle->SetOrthoMode(CustomInteractorStyle::OrthoMode::Sagittal);
    sagittalinteractorstyle->SetOriginFunc(originfunc);
    sagittalinteractorstyle->SetDefaultRenderer(
                m_orthoViews->GetSagittalView()->GetRenderers()->GetFirstRenderer());
    m_sagittalWidget->GetInteractor()->SetInteractorStyle(sagittalinteractorstyle);

    auto coronalinteractorstyle = vtkSmartPointer<CustomInteractorStyle>::New();
    coronalinteractorstyle->SetOrthoViews(m_orthoViews);
    coronalinteractorstyle->SetOrthoMode(CustomInteractorStyle::OrthoMode::Coronal);
    coronalinteractorstyle->SetOriginFunc(originfunc);
    coronalinteractorstyle->SetDefaultRenderer(
                m_orthoViews->GetCoronalView()->GetRenderers()->GetFirstRenderer());
    m_coronalWidget->GetInteractor()->SetInteractorStyle(coronalinteractorstyle);

    auto transverseinteractorstyle = vtkSmartPointer<CustomInteractorStyle>::New();
    transverseinteractorstyle->SetOrthoViews(m_orthoViews);
    transverseinteractorstyle->SetOrthoMode(CustomInteractorStyle::OrthoMode::Transverse);
    transverseinteractorstyle->SetOriginFunc(originfunc);
    transverseinteractorstyle->SetDefaultRenderer(
                m_orthoViews->GetTransverseView()->GetRenderers()->GetFirstRenderer());
    m_transverseWidget->GetInteractor()->SetInteractorStyle(transverseinteractorstyle);

    if (m_volumes.size())
        m_orthoViews->SetData(m_volumes[0]);

    SetupFilenamesWindow setupform(this);

    for (int i = 0; i < setupform.GetNumberOfShapes(); i++)
    {
        try
        {
            std::string filename, displacementsbasename;
            float red, green, blue;
            setupform.GetShape(i, filename, displacementsbasename, red, green, blue);

            std::string fullname = directory + "/" + filename;
            std::shared_ptr<ViewShape> shape = std::make_shared<ViewShape>(fullname, transformation, fullname);
            shape->SetColor(red, green, blue);
            shape->SetLineWidth(5.0);

            std::string displacementsfilename = directory + "/" + displacementsbasename + "_displacements.txt";
            std::vector<double> displacements = MainWindow::ReadDisplacements(displacementsfilename, shape->GetNumberOfVertices());

            for (int v = 0; v < shape->GetNumberOfVertices(); v++)
                shape->SetDisplacement(v, displacements[v]);

            m_orthoViews->AddShape(shape);
        }
        catch (Shape::ShapeException& e)
        {
            std::ostringstream message;
            message << "Loading of shape " << i + 1 << " failed. The exception text was:"
                       << std::endl << e.what();
            QMessageBox::warning(this, "Error", message.str().c_str());
            return;
        }
        catch (std::ios_base::failure &e)
        {
            std::ostringstream message;
            message << "Loading of displacements failed for shape " << i + 1 << ". The exception text was:"
                       << std::endl << e.what();
            QMessageBox::warning(this, "Error", message.str().c_str());
            return;
        }
    }

    SetupInterface();

    if (m_volumes.size())
        SelectVolume(0);

    m_orthoViews->ResetCameras();
}

MultiShapeWindow::~MultiShapeWindow()
{
    m_orthoViews.reset();

    delete ui;
}

void MultiShapeWindow::SetupInterface()
{
    ColumnVector bounds = m_orthoViews->GetBounds();
    ColumnVector origin(3);
    origin << (bounds(2) - bounds(1)) / 2.0
           << (bounds(4) - bounds(3)) / 2.0
           << (bounds(6) - bounds(5)) / 2.0;
    m_orthoViews->SetOrigin(origin);

    BLOCK_SIGNALS(ui->sagittalSlider, setMinimum(bounds(1) * m_sliderSubDivision))
    BLOCK_SIGNALS(ui->sagittalSlider, setMaximum(bounds(2) * m_sliderSubDivision))
    BLOCK_SIGNALS(ui->sagittalSlider, setSliderPosition(origin(1) * m_sliderSubDivision))
    BLOCK_SIGNALS(ui->coronalSlider, setMinimum(bounds(3) * m_sliderSubDivision))
    BLOCK_SIGNALS(ui->coronalSlider, setMaximum(bounds(4) * m_sliderSubDivision))
    BLOCK_SIGNALS(ui->coronalSlider, setSliderPosition(origin(2) * m_sliderSubDivision))
    BLOCK_SIGNALS(ui->transverseSlider, setMinimum(bounds(5) * m_sliderSubDivision))
    BLOCK_SIGNALS(ui->transverseSlider, setMaximum(bounds(6) * m_sliderSubDivision))
    BLOCK_SIGNALS(ui->transverseSlider, setSliderPosition(origin(3) * m_sliderSubDivision))

    BLOCK_SIGNALS(ui->modalityComboBox, clear());

    for (auto it = m_volumes.begin(); it != m_volumes.end(); ++it)
        BLOCK_SIGNALS(ui->modalityComboBox, addItem(QString((*it)->GetDisplayName().c_str())));
}

void MultiShapeWindow::SetOrigin(const ColumnVector& origin)
{
    BLOCK_SIGNALS(ui->sagittalSlider, setValue(static_cast<int>(origin(1) * m_sliderSubDivision)))
    BLOCK_SIGNALS(ui->coronalSlider, setValue(static_cast<int>(origin(2) * m_sliderSubDivision)))
    BLOCK_SIGNALS(ui->transverseSlider,setValue(static_cast<int>(origin(3) * m_sliderSubDivision)))

    if (m_orthoViews)
        m_orthoViews->SetOrigin(origin);
}

void MultiShapeWindow::SelectVolume(int index)
{
    m_orthoViews->SetData(m_volumes[index]);
    BLOCK_SIGNALS(ui->modalityComboBox, setCurrentIndex(index))
}

void MultiShapeWindow::on_sagittalSlider_valueChanged(int position)
{
    ColumnVector origin = m_orthoViews->GetOrigin();
    origin(1) = static_cast<double>(position) /  m_sliderSubDivision;
    SetOrigin(origin);
}

void MultiShapeWindow::on_coronalSlider_valueChanged(int position)
{
    ColumnVector origin = m_orthoViews->GetOrigin();
    origin(2) = static_cast<double>(position) /  m_sliderSubDivision;
    SetOrigin(origin);
}

void MultiShapeWindow::on_transverseSlider_valueChanged(int position)
{
    ColumnVector origin = m_orthoViews->GetOrigin();
    origin(3) = static_cast<double>(position) /  m_sliderSubDivision;
    SetOrigin(origin);
}

void MultiShapeWindow::on_modalityComboBox_currentIndexChanged(int index)
{
    if (index >= 0)
        SelectVolume(index);
}

void MultiShapeWindow::on_hideOverlaysComboBox_stateChanged(int value)
{
    m_orthoViews->SetHideOverlays(value ? OrthoViews::HideAllOverlays : OrthoViews::HideNone);
}
