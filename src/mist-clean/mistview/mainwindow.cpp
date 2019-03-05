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

#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "setupfilenameswindow.h"
#include "multishapewindow.h"
#include "orthoviews.h"
#include "viewdata.h"
#include "viewshape.h"
#include "custominteractorstyle.h"
#include "nokeysqvtkwidget.h"
#include "transformation.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindow.h"
#include "vtkCamera.h"
#include "vtkPNGWriter.h"
#include "vtkWindowToImageFilter.h"
#include "newimageall.h"
#include <QFileDialog>
#include <QFileInfo>
#include <QKeyEvent>
#include <QMessageBox>
#include <QString>
#include <QVTKWidget.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "boost/make_shared.hpp"

#define BLOCK_SIGNALS(object, call) { bool oldstate = object->blockSignals(true); object->call; object->blockSignals(oldstate); }

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    m_plotWindow(this)
{
    ui->setupUi(this);

    setWindowFlags(Qt::Window | Qt::CustomizeWindowHint | Qt::WindowCloseButtonHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint);
    GetFilenames(false);

    m_statusLabel.setParent(this);
    statusBar()->addWidget(&m_statusLabel);

    auto sagittalwidget = new NoKeysQVTKWidget(this);
    sagittalwidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    sagittalwidget->setMinimumSize(m_minimumViewWidth, m_minimumViewHeight);
    ui->sagittalLayout->insertWidget(0, sagittalwidget);

    auto coronalwidget = new NoKeysQVTKWidget(this);
    coronalwidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    coronalwidget->setMinimumSize(m_minimumViewWidth, m_minimumViewHeight);
    ui->coronalLayout->insertWidget(0, coronalwidget);

    auto transversewidget = new NoKeysQVTKWidget(this);
    transversewidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    transversewidget->setMinimumSize(m_minimumViewWidth, m_minimumViewHeight);
    ui->transverseLayout->insertWidget(0, transversewidget);

    m_orthoViews = make_shared<OrthoViews>(sagittalwidget->GetRenderWindow(),
                                           coronalwidget->GetRenderWindow(),
                                           transversewidget->GetRenderWindow());
    m_orthoViews->SetScale(2.0);

    auto originfunc = [this] (const ColumnVector& origin)
                {
                    SetOrigin(origin);
                };

    auto sagittalinteractorstyle = vtkSmartPointer<CustomInteractorStyle>::New();
    sagittalinteractorstyle->SetOrthoViews(m_orthoViews);
    sagittalinteractorstyle->SetOrthoMode(CustomInteractorStyle::OrthoMode::Sagittal);
    sagittalinteractorstyle->SetOriginFunc(originfunc);
    sagittalinteractorstyle->SetPickerFunc(
                [this] (double x, double y)
                {
                    PickSagittal(x, y);
                });
    sagittalinteractorstyle->SetDefaultRenderer(
                m_orthoViews->GetSagittalView()->GetRenderers()->GetFirstRenderer());
    sagittalwidget->GetInteractor()->SetInteractorStyle(sagittalinteractorstyle);

    auto coronalinteractorstyle = vtkSmartPointer<CustomInteractorStyle>::New();
    coronalinteractorstyle->SetOrthoViews(m_orthoViews);
    coronalinteractorstyle->SetOrthoMode(CustomInteractorStyle::OrthoMode::Coronal);
    coronalinteractorstyle->SetOriginFunc(originfunc);
    coronalinteractorstyle->SetPickerFunc(
                [this] (double x, double y)
                {
                    PickCoronal(x, y);
                });
    coronalinteractorstyle->SetDefaultRenderer(
                m_orthoViews->GetCoronalView()->GetRenderers()->GetFirstRenderer());
    coronalwidget->GetInteractor()->SetInteractorStyle(coronalinteractorstyle);

    auto transverseinteractorstyle = vtkSmartPointer<CustomInteractorStyle>::New();
    transverseinteractorstyle->SetOrthoViews(m_orthoViews);
    transverseinteractorstyle->SetOrthoMode(CustomInteractorStyle::OrthoMode::Transverse);
    transverseinteractorstyle->SetOriginFunc(originfunc);
    transverseinteractorstyle->SetPickerFunc(
                [this] (double x, double y)
                {
                    PickTransverse(x, y);
                });
    transverseinteractorstyle->SetDefaultRenderer(
                m_orthoViews->GetTransverseView()->GetRenderers()->GetFirstRenderer());
    transversewidget->GetInteractor()->SetInteractorStyle(transverseinteractorstyle);

    EnableInterface(false);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::CreateNormalViews()
{
    auto camera = vtkSmartPointer<vtkCamera>::New();

    for (auto it = m_volumes.begin(); it != m_volumes.end(); ++it)
    {
        std::unique_ptr<NoKeysQVTKWidget> normalwidget(new NoKeysQVTKWidget(this));
        normalwidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        normalwidget->setMinimumSize(m_minimumViewWidth, m_minimumViewHeight);
        ui->normalViewsLayout->addWidget(normalwidget.get());

        auto normalview = make_shared<NormalView>(camera, normalwidget->GetRenderWindow());
        normalview->SetScale(2.0);

        auto normalinteractorstyle = vtkSmartPointer<CustomInteractorStyle>::New();
        normalinteractorstyle->SetOriginFunc(
                    [this] (const ColumnVector& origin)
                    {
                        SetOrigin(origin);
                    });
        normalinteractorstyle->SetPickerFunc(
                    [this, normalview] (double x, double y)
                    {
                        PickNormal(normalview, x, y);
                    });

        normalwidget->GetInteractor()->SetInteractorStyle(normalinteractorstyle);
        normalwidget->GetInteractor()->AddObserver(vtkCommand::EndInteractionEvent, this, &MainWindow::NormalViewInteractionEventHandler);
//        normalview->SetupOrientationWidget(normalwidget->GetInteractor());

        m_normalViews.push_back(normalview);
        m_normalViewWidgets.push_back(std::move(normalwidget));

        normalview->SetData(*it);
    }
}

void MainWindow::NormalViewInteractionEventHandler(vtkObject *caller, unsigned long eventId, void *callData)
{
    for (auto it = m_normalViewWidgets.begin(); it != m_normalViewWidgets.end(); ++it)
        if ((*it)->GetRenderWindow() != caller)
            (*it)->update();
}

void MainWindow::DestroyNormalViews()
{
    m_normalViewWidgets.clear();
    m_normalViews.clear();
}

void MainWindow::LoadVolume(std::string filename, std::string name)
{
    boost::shared_ptr<NEWIMAGE::volume<double> > vol =
            boost::make_shared<NEWIMAGE::volume<double> >();
    NEWIMAGE::read_volume(*vol, filename);
    m_volumes.push_back(std::make_shared<ViewData>(vol, name));
    m_plotWindow.SetVolume(name, vol);
}

void MainWindow::CloseSubject()
{
    m_plotWindow.Clear();

    m_orthoViews->ResetDataAndShapes();
    m_orthoViews->Render();
    DestroyNormalViews();

    m_model.reset();
    m_plotWindow.SetModel(nullptr);
    m_volumes.clear();
    m_shape.reset();
    m_extraShape.reset();
    m_overlay.reset();
    m_unsavedChanges = false;
    m_selectedDirectory = "";
    EnableInterface(false);
    m_statusLabel.setText("Data unloaded");
}

void MainWindow::SetupInterface()
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

    BLOCK_SIGNALS(ui->xMmSpinBox, setMinimum(bounds(1)))
    BLOCK_SIGNALS(ui->xMmSpinBox, setMaximum(bounds(2)))
    BLOCK_SIGNALS(ui->xMmSpinBox, setValue(origin(1)))
    BLOCK_SIGNALS(ui->yMmSpinBox, setMinimum(bounds(3)))
    BLOCK_SIGNALS(ui->yMmSpinBox, setMaximum(bounds(4)))
    BLOCK_SIGNALS(ui->yMmSpinBox, setValue(origin(2)))
    BLOCK_SIGNALS(ui->zMmSpinBox, setMinimum(bounds(5)))
    BLOCK_SIGNALS(ui->zMmSpinBox, setMaximum(bounds(6)))
    BLOCK_SIGNALS(ui->zMmSpinBox, setValue(origin(3)))

    auto current = ui->modalityComboBox->currentIndex();
    BLOCK_SIGNALS(ui->modalityComboBox, clear());

    for (auto it = m_volumes.begin(); it != m_volumes.end(); ++it)
        BLOCK_SIGNALS(ui->modalityComboBox, addItem(QString((*it)->GetDisplayName().c_str())));

    if (current < 0)
        current = 0;
    else if (current >= m_volumes.size())
        current = m_volumes.size();

    BLOCK_SIGNALS(ui->modalityComboBox, setCurrentIndex(current))

    BLOCK_SIGNALS(ui->windowMinSpinBox, setValue(m_volumes[current]->GetDisplayMin()))
    BLOCK_SIGNALS(ui->windowMaxSpinBox, setValue(m_volumes[current]->GetDisplayMax()))

    BLOCK_SIGNALS(ui->vertexIdSpinBox, setMinimum(0));
    if (m_shape)
    {
        BLOCK_SIGNALS(ui->vertexIdSpinBox, setMaximum(m_shape->GetNumberOfVertices() - 1))
        BLOCK_SIGNALS(ui->shapeNameLabel, setText(m_shape->GetDisplayName().c_str()))
        BLOCK_SIGNALS(ui->adjustMmSpinBox, setValue(m_shape->GetDisplacement(ui->vertexIdSpinBox->value())))
    }
    else
        BLOCK_SIGNALS(ui->vertexIdSpinBox, setMaximum(0))

    BLOCK_SIGNALS(ui->slabThicknessSpinBox, setValue(m_orthoViews->GetMarkerTolerance()));
    BLOCK_SIGNALS(ui->hideOverlaysCheckBox, setChecked(m_orthoViews->GetHideOverlays()));
}

void MainWindow::EnableInterface(bool enable)
{
    ui->sagittalSlider->setEnabled(enable);
    ui->coronalSlider->setEnabled(enable);
    ui->transverseSlider->setEnabled(enable);
    ui->xMmSpinBox->setEnabled(enable);
    ui->yMmSpinBox->setEnabled(enable);
    ui->zMmSpinBox->setEnabled(enable);
    ui->modalityComboBox->setEnabled(enable);
    ui->windowMinSpinBox->setEnabled(enable);
    ui->windowMaxSpinBox->setEnabled(enable);
    ui->vertexIdSpinBox->setEnabled(enable);
    ui->adjustMmSpinBox->setEnabled(enable);
    ui->slabThicknessSpinBox->setEnabled(enable);
    ui->hideOverlaysCheckBox->setEnabled(enable);
    ui->saveAction->setEnabled(enable);
    ui->loadAction->setEnabled(enable);
    ui->showProfileFitAction->setEnabled(enable);
    ui->showMultiShapeAction->setEnabled(enable);
}

void MainWindow::SetOrigin(const ColumnVector& origin)
{
    BLOCK_SIGNALS(ui->xMmSpinBox, setValue(origin(1)))
    BLOCK_SIGNALS(ui->sagittalSlider, setValue(static_cast<int>(origin(1) * m_sliderSubDivision)))
    BLOCK_SIGNALS(ui->yMmSpinBox, setValue(origin(2)))
    BLOCK_SIGNALS(ui->coronalSlider, setValue(static_cast<int>(origin(2) * m_sliderSubDivision)))
    BLOCK_SIGNALS(ui->zMmSpinBox, setValue(origin(3)))
    BLOCK_SIGNALS(ui->transverseSlider,setValue(static_cast<int>(origin(3) * m_sliderSubDivision)))

    if (m_orthoViews)
        m_orthoViews->SetOrigin(origin);
}

void MainWindow::SelectVertex(vtkIdType vertex)
{
    // Check range or set on widget
    SetOrigin(m_shape->GetPoint(vertex));
    ui->vertexIdSpinBox->setValue(vertex);
    ui->adjustMmSpinBox->setValue(m_shape->GetDisplacement(vertex));
    m_orthoViews->SetMarker(m_shape, vertex);
    for (auto it = m_normalViews.begin(); it != m_normalViews.end(); ++it)
    {
        (*it)->SetMarker(m_shape, vertex);
        (*it)->ResetCameras();
    }

    if (m_plotWindow.isVisible())
        m_plotWindow.ShowVertex(vertex, m_shape->GetDisplacement(vertex));
}

void MainWindow::SelectVolume(int index)
{
    m_orthoViews->SetData(m_volumes[index]);
    BLOCK_SIGNALS(ui->modalityComboBox, setCurrentIndex(index))
    BLOCK_SIGNALS(ui->windowMinSpinBox, setValue(m_volumes[index]->GetDisplayMin()))
    BLOCK_SIGNALS(ui->windowMaxSpinBox, setValue(m_volumes[index]->GetDisplayMax()))
}

void MainWindow::PickSagittal(double x, double y)
{
    vtkIdType pointid = m_orthoViews->PickSagittal(x, y);

    if (pointid > -1)
        SelectVertex(pointid);
}

void MainWindow::PickCoronal(double x, double y)
{
    vtkIdType pointid = m_orthoViews->PickCoronal(x, y);

    if (pointid > -1)
        SelectVertex(pointid);
}

void MainWindow::PickTransverse(double x, double y)
{
    vtkIdType pointid = m_orthoViews->PickTransverse(x, y);

    if (pointid > -1)
        SelectVertex(pointid);
}

void MainWindow::PickNormal(const std::shared_ptr<NormalView>& normalView, double x, double y)
{
    vtkIdType pointid = normalView->Pick(x, y);

    if (pointid > -1)
        SelectVertex(pointid);
}

void MainWindow::keyPressEvent(QKeyEvent* event)
{
    switch (event->key())
    {
        case Qt::Key_1:
            SelectVolume(std::min(0, static_cast<int>(m_volumes.size() - 1)));
            return;

        case Qt::Key_2:
            SelectVolume(std::min(1, static_cast<int>(m_volumes.size() - 1)));
            return;

        case Qt::Key_3:
            SelectVolume(std::min(2, static_cast<int>(m_volumes.size() - 1)));
            return;

        case Qt::Key_4:
            SelectVolume(std::min(3, static_cast<int>(m_volumes.size() - 1)));
            return;

        case Qt::Key_5:
            SelectVolume(std::min(4, static_cast<int>(m_volumes.size() - 1)));
            return;

        case Qt::Key_6:
            SelectVolume(std::min(5, static_cast<int>(m_volumes.size() - 1)));
            return;

        case Qt::Key_7:
            SelectVolume(std::min(6, static_cast<int>(m_volumes.size() - 1)));
            return;

        case Qt::Key_8:
            SelectVolume(std::min(7, static_cast<int>(m_volumes.size() - 1)));
            return;

        case Qt::Key_9:
            SelectVolume(std::min(8, static_cast<int>(m_volumes.size() - 1)));
            return;

        case Qt::Key_Minus:
            if (ui->adjustMmSpinBox->value() >= ui->adjustMmSpinBox->minimum() + ui->adjustMmSpinBox->singleStep())
               ui->adjustMmSpinBox->setValue(ui->adjustMmSpinBox->value() - ui->adjustMmSpinBox->singleStep());
            return;

        case Qt::Key_Equal:
            if (ui->adjustMmSpinBox->value() <= ui->adjustMmSpinBox->maximum() - ui->adjustMmSpinBox->singleStep())
              ui->adjustMmSpinBox->setValue(ui->adjustMmSpinBox->value() + ui->adjustMmSpinBox->singleStep());
            return;

        case Qt::Key_BracketLeft:
            if (ui->vertexIdSpinBox->value() > ui->vertexIdSpinBox->minimum())
                ui->vertexIdSpinBox->setValue(ui->vertexIdSpinBox->value() - 1);
            return;

        case Qt::Key_BracketRight:
            if (ui->vertexIdSpinBox->value() < ui->vertexIdSpinBox->maximum())
                ui->vertexIdSpinBox->setValue(ui->vertexIdSpinBox->value() + 1);
            return;

        case Qt::Key_Slash:
            ui->hideOverlaysCheckBox->toggle();
            return;

        default:
            ;
    }

    QMainWindow::keyPressEvent(event);
}

void MainWindow::on_sagittalSlider_valueChanged(int position)
{
    ColumnVector origin = m_orthoViews->GetOrigin();
    origin(1) = static_cast<double>(position) /  m_sliderSubDivision;
    SetOrigin(origin);
}

void MainWindow::on_coronalSlider_valueChanged(int position)
{
    ColumnVector origin = m_orthoViews->GetOrigin();
    origin(2) = static_cast<double>(position) /  m_sliderSubDivision;
    SetOrigin(origin);
}

void MainWindow::on_transverseSlider_valueChanged(int position)
{
    ColumnVector origin = m_orthoViews->GetOrigin();
    origin(3) = static_cast<double>(position) /  m_sliderSubDivision;
    SetOrigin(origin);
}

void MainWindow::on_xMmSpinBox_valueChanged(double value)
{
    ColumnVector origin = m_orthoViews->GetOrigin();
    origin(1) = value;
    SetOrigin(origin);
}

void MainWindow::on_yMmSpinBox_valueChanged(double value)
{
    ColumnVector origin = m_orthoViews->GetOrigin();
    origin(2) = value;
    SetOrigin(origin);
}

void MainWindow::on_zMmSpinBox_valueChanged(double value)
{
    ColumnVector origin = m_orthoViews->GetOrigin();
    origin(3) = value;
    SetOrigin(origin);
}

void MainWindow::on_vertexIdSpinBox_valueChanged(int value)
{
    SelectVertex(value);
}

void MainWindow::on_modalityComboBox_currentIndexChanged(int index)
{
    if (index >= 0)
        SelectVolume(index);
}

void MainWindow::on_adjustMmSpinBox_valueChanged(double amount)
{
    m_unsavedChanges = true;
    m_statusLabel.setText("There are unsaved changes");
    m_shape->SetDisplacement(ui->vertexIdSpinBox->value(), amount);

    m_orthoViews->RedrawMarkers();
    for (auto it = m_normalViews.begin(); it != m_normalViews.end(); ++it)
        (*it)->RedrawMarkers();
}

void MainWindow::on_windowMinSpinBox_valueChanged(double value)
{
    m_volumes[ui->modalityComboBox->currentIndex()]->SetDisplayMin(value);
    m_orthoViews->Render();
    for (auto it = m_normalViews.begin(); it != m_normalViews.end(); ++it)
        (*it)->RedrawMarkers();
}

void MainWindow::on_windowMaxSpinBox_valueChanged(double value)
{
    m_volumes[ui->modalityComboBox->currentIndex()]->SetDisplayMax(value);
    m_orthoViews->Render();
    for (auto it = m_normalViews.begin(); it != m_normalViews.end(); ++it)
        (*it)->RedrawMarkers();
}

void MainWindow::on_setupFilenamesAction_triggered()
{
    GetFilenames(true);
}

void MainWindow::on_openAction_triggered()
{
    QFileDialog dialog(this);
    dialog.setFileMode(QFileDialog::FileMode::Directory);

    if (dialog.exec() == QFileDialog::Accepted)
    {
        CloseSubject();

        std::string directory = dialog.selectedFiles().at(0).toStdString();

        std::string modelfile = directory + "/" + m_modelFilename;
        std::string transformfile = directory + "/" + m_transformFilename;
        std::string displacementsfile = directory + "/" + m_displacementsBasename + "_displacements.txt";
        std::string normexclusionfile = directory + "/" + m_normExclusionFilename;
        std::string extrashapefile = directory + "/" + m_extraShapeFilename;
        std::string overlayfilename = directory + "/" + m_overlayFilename;

        boost::shared_ptr<const Transformation> transformation;

        class LocalRethrowException : public std::runtime_error { using std::runtime_error::runtime_error; };

        try
        {
            if (m_volumeFilenames.empty())
                throw LocalRethrowException("Error: List of volumes is empty");

            if (NEWIMAGE::fsl_imageexists(transformfile))
            {
                try
                {
                    transformation = boost::make_shared<NonlinearTransformation>(transformfile, true);
                }
                catch (RBD_COMMON::Exception& e)
                {
                    throw LocalRethrowException(std::string("Couldn't load nonlinear transformation. The exception text was:\n") + e.what());
                }
            }
            else
            {
                try
                {
                    transformation = boost::make_shared<AffineTransformation>(transformfile);
                }
                catch (std::exception &e)
                {
                    throw LocalRethrowException(std::string("Couldn't load affine transformation (assuming affine as no warp exists with the specified filename). The exception text was:\n") + e.what());
                }
            }

            m_transformation = transformation;

            try
            {
                std::ifstream is(modelfile);
                boost::archive::text_iarchive ar(is);
                ar >> m_model;
            }
            catch (std::exception &e)
            {
                throw LocalRethrowException(std::string("Couldn't load model. The exception text was:\n") + e.what());
            }
        
            m_shape = std::make_shared<ViewShape>(*m_model->GetShape(), transformation);

            try
            {
                std::vector<double> displacements = ReadDisplacements(displacementsfile, m_shape->GetNumberOfVertices());
                // If the displacements file is too short, ReadDisplacements() will have thrown an exception.
                for (int i = 0; i < m_shape->GetNumberOfVertices(); i++)
                    m_shape->SetDisplacement(i, displacements[i]);
            }
            catch (std::exception &e)
            {
                throw LocalRethrowException(std::string("Couldn't load displacements. The exception text was:\n") + e.what());
            }

            if (!m_extraShapeFilename.empty())
            {
                try
                {
                    boost::shared_ptr<const Transformation> identity = boost::make_shared<IdentityTransformation>();
                    m_extraShape = std::make_shared<ViewShape>(extrashapefile, identity, extrashapefile);
                }
                catch (std::exception &e)
                {
                    throw LocalRethrowException(std::string("Couldn't load extra shape. The exception text was:\n") + e.what());
                }
            }

            if (!m_overlayFilename.empty())
            {
                try
                {
                    auto overlayvol = boost::make_shared<NEWIMAGE::volume<float> >();
                    NEWIMAGE::read_volume(*overlayvol, overlayfilename);
                    m_overlay = std::make_shared<ViewData>(overlayvol, "Overlay");
                    m_overlay->SetBinaryLut(0.33, 1.0);
                }
                catch (RBD_COMMON::Exception& e)
                {
                    throw LocalRethrowException(std::string("Couldn't load overlay. The exception text was:\n") + e.what());
                }
            }

            for (auto it = m_volumeFilenames.begin(); it != m_volumeFilenames.end(); ++it)
            {
                try
                {
                    LoadVolume(directory + "/" + it->filename, it->name);
                }
                catch (RBD_COMMON::Exception& e)
                {
                    throw LocalRethrowException(std::string("Loading of volume ") + it->filename + " failed. The exception text was:\n" + e.what());
                }
            }

            if (!m_normExclusionFilename.empty())
            {
                try
                {
                    boost::shared_ptr<NEWIMAGE::volume<double> > normexclusionmask =
                            boost::make_shared<NEWIMAGE::volume<double> >();
                    NEWIMAGE::read_volume(*normexclusionmask, normexclusionfile);

                    m_plotWindow.SetNormExclusionVolume(normexclusionmask);
                }
                catch (RBD_COMMON::Exception& e)
                {
                    throw LocalRethrowException(std::string("Loading of normalisation exclusion mask ") + normexclusionfile
                                                + " failed. The exception text was:\n" + e.what());
                }
            }
        }
        catch (LocalRethrowException &e)
        {
            QMessageBox::warning(this, "Error", e.what());
            m_statusLabel.setText("Load failed");

            return;
        }

        m_plotWindow.SetModel(m_model);
        m_plotWindow.SetTransformation(m_transformation);

        CreateNormalViews();

        m_selectedDirectory = directory;
        SelectVolume(0);
        SetupInterface();
        EnableInterface(true);

        m_shape->SetColor(1.0, 0.0, 0.0);
        m_shape->SetLineWidth(7.0);
        m_orthoViews->AddShape(m_shape);
        m_orthoViews->SetMarker(m_shape, 0);

        if (m_extraShape)
        {
            m_extraShape->SetColor(0.0, 0.0, 1.0);
            m_extraShape->SetLineWidth(7.0);
            m_orthoViews->AddShape(m_extraShape);
        }

        if (m_overlay)
            m_orthoViews->SetOverlay(m_overlay);

        for (auto it = m_normalViews.begin(); it != m_normalViews.end(); ++it)
        {
            (*it)->AddShape(m_shape);
            (*it)->SetMarker(m_shape, 0);

            if (m_extraShape)
                (*it)->AddShape(m_extraShape);

            if (m_overlay)
                (*it)->SetOverlay(m_overlay);
        }

        m_orthoViews->ResetCameras();

        for (auto it = m_normalViews.begin(); it != m_normalViews.end(); ++it)
            (*it)->ResetCameras();

        m_statusLabel.setText("Data loaded");
    }
}

void MainWindow::on_closeAction_triggered()
{
    CloseSubject();
}

void MainWindow::GetFilenames(bool showdialog)
{
    SetupFilenamesWindow setupform(this);

    if (!showdialog || (setupform.exec() == QDialog::Accepted))
    {
        m_modelFilename = setupform.GetModelFilename();
        m_transformFilename = setupform.GetTransformFilename();
        m_displacementsBasename = setupform.GetDisplacementsBasename();
        m_normExclusionFilename = setupform.GetNormExclusionFilename();
        m_outputFilename = setupform.GetOutputFilename();
        m_extraShapeFilename = setupform.GetExtraShapeFilename();
        m_overlayFilename = setupform.GetOverlayFilename();
        m_volumeFilenames.clear();
        for (int i = 0; i < setupform.GetNumberOfVolumes(); i++)
        {
            std::string name, filename;
            setupform.GetVolumeFilename(i, name, filename);
            m_volumeFilenames.push_back(Modality({name, filename}));
        }
    }
}

void MainWindow::WriteDisplacements(std::string filename)
{
    std::ofstream out;
    out.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    out.open(filename);
    out << std::fixed;
    for (vtkIdType i = 0; i < m_shape->GetNumberOfVertices(); i++)
        out << m_shape->GetDisplacement(i) << std::endl;
    out.close();
}

bool MainWindow::SaveInteractive()
{
    bool success = false;

    std::string proposal = m_selectedDirectory + "/adjustments.txt";
    // Native dialog won't show suggested filename due to a bug, so use the Qt one for now
    std::string filename = QFileDialog::getSaveFileName(
                this, "Save displacements", proposal.c_str(), "Text files (*.txt)", 0, QFileDialog::DontUseNativeDialog).toStdString();

    if (filename.length())
    try
    {
        WriteDisplacements(filename);
        success = true;
    }
    catch (std::ios_base::failure &e)
    {
        std::ostringstream message;
        message << "Save failed. The exception text was:"
                   << std::endl << e.what();
        QMessageBox::warning(this, "Error", message.str().c_str());
    }

    return success;
}

std::vector<double> MainWindow::ReadDisplacements(std::string filename, int vertices)
{
    std::ifstream in;
    in.exceptions(std::ofstream::failbit | std::ofstream::badbit);
    in.open(filename);
    std::vector<double> displacements;
    for (vtkIdType i = 0; i < vertices; i++)
    {
        double val;
        in >> val;
        displacements.push_back(val);
    }
    in.close();

    return displacements;
}

void MainWindow::LoadInteractive()
{
    if (m_unsavedChanges
            && (QMessageBox::question(this, "Quit?", "There are unsaved changes. Save these first?",
                                      QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes) == QMessageBox::Yes))
    {
        if (!SaveInteractive())
            return; // Save failed; don't continue loading data
    }

    m_unsavedChanges = false;

    std::string filename = QFileDialog::getOpenFileName(
                this, "Load displacements", m_selectedDirectory.c_str(), "Text files (*.txt)").toStdString();

    if (filename.length())
    {
        vector<double> displacements;
        try
        {
            displacements = ReadDisplacements(filename, m_shape->GetNumberOfVertices());

            for (vtkIdType i = 0; i < m_shape->GetNumberOfVertices(); i++)
                m_shape->SetDisplacement(i, displacements[i]);

            m_statusLabel.setText("Load complete");
        }
        catch (std::ios_base::failure &e)
        {
            std::ostringstream message;
            message << "Load failed. The exception text was:"
                       << std::endl << e.what();
            QMessageBox::warning(this, "Error", message.str().c_str());

            m_statusLabel.setText("Load failed");
        }
    }

    SelectVertex(ui->vertexIdSpinBox->value());
}

void MainWindow::on_saveAction_triggered()
{
    if (SaveInteractive())
    {
        m_unsavedChanges = false;
        m_statusLabel.setText("Save complete");
    }
}

void MainWindow::on_loadAction_triggered()
{
    LoadInteractive();
}

void MainWindow::on_quitAction_triggered()
{
    if (m_unsavedChanges
            && (QMessageBox::question(this, "Quit?", "There are unsaved changes. Save these first?",
                                      QMessageBox::Yes | QMessageBox::No, QMessageBox::Yes) == QMessageBox::Yes))
    {
        if (!SaveInteractive())
            return; // Save failed; don't quit
    }

    CloseSubject();

    m_orthoViews.reset();

    close();
}


void MainWindow::on_slabThicknessSpinBox_valueChanged(double value)
{
    m_orthoViews->SetMarkerTolerance(value);
    for (auto it = m_normalViews.begin(); it != m_normalViews.end(); ++it)
        (*it)->SetMarkerTolerance(value);
}

void MainWindow::on_hideOverlaysCheckBox_stateChanged(int value)
{
    ViewBase::HideMode mode;

    if (value)
        mode = ui->markersOnlyCheckBox->checkState() ? ViewBase::HideMarkers : ViewBase::HideAllOverlays;
    else
        mode = ViewBase::HideNone;

    m_orthoViews->SetHideOverlays(mode);
    for (auto it = m_normalViews.begin(); it != m_normalViews.end(); ++it)
        (*it)->SetHideOverlays(mode);
}

void MainWindow::on_showProfileFitAction_triggered()
{
    m_plotWindow.show();

    int vertex = ui->vertexIdSpinBox->value();
    m_plotWindow.ShowVertex(vertex, m_shape->GetDisplacement(vertex));
}

void MainWindow::on_showMultiShapeAction_triggered()
{
    MultiShapeWindow ms(this, m_volumes, m_transformation, m_selectedDirectory);

    ms.exec();
}

void MainWindow::on_saveSlicesPNGAction_triggered()
{
    QString basename = QFileDialog::getSaveFileName(this, "Base filename for slices", "", "PNG files (*.png)");

    int slices = 10;

    if (basename.length())
    {
        QFileInfo fileinfo(basename);

        ViewBase::HideMode oldhidemode = m_orthoViews->GetHideOverlays();
        ColumnVector oldorigin = m_orthoViews->GetOrigin();
        int oldvolume = ui->modalityComboBox->currentIndex();

        auto pngwriter = vtkSmartPointer<vtkPNGWriter>::New();
        auto imfilt = vtkSmartPointer<vtkWindowToImageFilter>::New();
        pngwriter->SetInputConnection(imfilt->GetOutputPort());

        ColumnVector bounds = m_shape->GetBounds();

        for (auto hm : std::vector<ViewBase::HideMode>({ViewBase::HideMarkers, ViewBase::HideAllOverlays}))
        {
            for (int i = 0; i < m_volumes.size(); i++)
            {
                SelectVolume(i);
                m_orthoViews->SetHideOverlays(hm);

                imfilt->SetInput(m_orthoViews->GetSagittalView());
                int imnum = 0;
                for (double x = bounds(1); x <= bounds(2); x += (bounds(2) - bounds(1)) / (slices - 1))
                {
                    ColumnVector origin(3);
                    origin << x << bounds(3) << bounds(5);
                    m_orthoViews->SetOrigin(origin);
                    imfilt->Modified();
                    std::string filename = (fileinfo.absolutePath() + "/" + fileinfo.baseName() + "_" + m_volumes[i]->GetDisplayName().c_str()
                                            + "_sag" + QString::number(imnum++) + (hm == ViewBase::HideAllOverlays ? "_nocontour" : "" ) + ".png").toStdString();
                    pngwriter->SetFileName(filename.c_str());
                    pngwriter->Write();
                    pngwriter->SetFileName(nullptr);
                }

                imfilt->SetInput(m_orthoViews->GetCoronalView());
                imnum = 0;
                for (double y = bounds(3); y <= bounds(4); y += (bounds(4) - bounds(3)) / (slices - 1))
                {
                    ColumnVector origin(3);
                    origin << bounds(1) << y << bounds(5);
                    m_orthoViews->SetOrigin(origin);
                    imfilt->Modified();
                    std::string filename = (fileinfo.absolutePath() + "/" + fileinfo.baseName() + "_" + m_volumes[i]->GetDisplayName().c_str()
                                            + "_cor" + QString::number(imnum++) + (hm == ViewBase::HideAllOverlays ? "_nocontour" : "" ) + ".png").toStdString();
                    pngwriter->SetFileName(filename.c_str());
                    pngwriter->Write();
                    pngwriter->SetFileName(nullptr);
                }

                imfilt->SetInput(m_orthoViews->GetTransverseView());
                imnum = 0;
                for (double z = bounds(5); z <= bounds(6); z += (bounds(6) - bounds(5)) / (slices - 1))
                {
                    ColumnVector origin(3);
                    origin << bounds(1) << bounds(3) << z;
                    m_orthoViews->SetOrigin(origin);
                    imfilt->Modified();
                    std::string filename = (fileinfo.absolutePath() + "/" + fileinfo.baseName() + "_" + m_volumes[i]->GetDisplayName().c_str()
                                            + "_tra" + QString::number(imnum++) + (hm == ViewBase::HideAllOverlays ? "_nocontour" : "" ) + ".png").toStdString();
                    pngwriter->SetFileName(filename.c_str());
                    pngwriter->Write();
                    pngwriter->SetFileName(nullptr);
                }
            }
        }

        m_orthoViews->SetHideOverlays(oldhidemode);
        m_orthoViews->SetOrigin(oldorigin);
        ui->modalityComboBox->setCurrentIndex(oldvolume);
    }
}


