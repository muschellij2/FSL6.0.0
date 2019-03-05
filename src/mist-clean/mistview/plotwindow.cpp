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

#include "plotwindow.h"
#include "ui_plotwindow.h"
#include "profilemixtures.h"
#include <boost/pointer_cast.hpp>
#include <boost/log/trivial.hpp>
#include <QVTKWidget.h>
#include <vtkSmartPointer.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>
#include <vtkPlotLine.h>
#include <vtkPlotStacked.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkTable.h>
#include <vtkFloatArray.h>
#include <vtkContextActor.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkAbstractContextItem.h>
#include <vtkColorSeries.h>
#include <vtkAxis.h>
#include <cmath>
#include <iostream>
#include <QShowEvent>
#include <QMessageBox>
#include <QFileDialog>

#define WANT_STREAM
#include "newmatio.h"

PlotWindow::PlotWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PlotWindow)
{
    ui->setupUi(this);
}

PlotWindow::~PlotWindow()
{
    delete ui;
}

void PlotWindow::showEvent(QShowEvent *event)
{
    QDialog::showEvent(event);
}

void PlotWindow::hideEvent(QHideEvent *event)
{
    if (!event->spontaneous())
    {
        m_contextViews.clear();
        m_model = nullptr;
    }

    QDialog::hideEvent(event);
}

void PlotWindow::SetVolume(std::string name, boost::shared_ptr<NEWIMAGE::volume<double> > vol)
{
    m_volumes[name] = vol;
}

void PlotWindow::SetNormExclusionVolume(boost::shared_ptr<NEWIMAGE::volume<double> > vol)
{
    m_normExclusionVolume = vol;
}

void PlotWindow::SetTransformation(boost::shared_ptr<const Transformation> &transformation)
{
    m_transformation = transformation;
}

void PlotWindow::ShowVertex(int vertex, double displacement)
{
    BOOST_LOG_TRIVIAL(debug) << "Showing model fit for vertex " << vertex;

    if (!m_model)
        return;

    try
    {
        auto pm = boost::dynamic_pointer_cast<const ProfileMixtures>(m_model->GetVertexModel(vertex));

        int components = pm->GetNumberOfComponents();
        int profilepoints = m_model->GetProfilePoints();
        double profilespacing = m_model->GetProfileSpacing();
        int steps = pm->GetNumberOfSteps();

        auto mixingcoefs = pm->GetFullMixingCoefs();
        auto means = pm->GetComponentMeans();
        auto covcoefs = pm->GetComponentCovarianceCoefs();
        auto smoothness = pm->GetSmoothness();

        auto names = m_model->GetModalityNames();

        for (auto name = names.cbegin(); name != names.cend(); name++)
        {
            Matrix G = ProfileMixtureProbability::GetSmoothingMatrix(smoothness[*name] * smoothness[*name], profilepoints + steps);

            for (int c = 0; c < components; c++)
            {
                Matrix C = G * covcoefs[*name][c].AsDiagonal() * G;

                auto xarr = vtkSmartPointer<vtkFloatArray>::New();
                xarr->SetNumberOfTuples(profilepoints + steps);
                xarr->SetName("x");
                auto y1arr = vtkSmartPointer<vtkFloatArray>::New();
                y1arr->SetNumberOfTuples(profilepoints + steps);
                y1arr->SetName("y1");
                auto y2arr = vtkSmartPointer<vtkFloatArray>::New();
                y2arr->SetNumberOfTuples(profilepoints + steps);
                y2arr->SetName("y2");
                auto y3arr = vtkSmartPointer<vtkFloatArray>::New();
                y3arr->SetNumberOfTuples(profilepoints + steps);
                y3arr->SetName("y3");

                for (int p = 0; p < profilepoints + steps; p++)
                {
                    xarr->SetTuple1(p, (p - (profilepoints + steps) / 2.0) * profilespacing);

                    double mean = means[*name][c](p + 1);
                    double stdev = std::sqrt(C(p + 1, p + 1));
                    y1arr->SetTuple1(p, mean - stdev);
                    y2arr->SetTuple1(p, stdev);
                    y3arr->SetTuple1(p, stdev);
                }

                auto table = vtkSmartPointer<vtkTable>::New();
                table->AddColumn(xarr);
                table->AddColumn(y1arr);
                table->AddColumn(y2arr);
                table->AddColumn(y3arr);

                auto colorseries = vtkSmartPointer<vtkColorSeries>::New();
                colorseries->SetColorScheme(vtkColorSeries::CUSTOM);
                colorseries->InsertColor(0, vtkColor3ub(255, 255, 255));
                colorseries->InsertColor(1, vtkColor3ub(255, 0, 0));
                colorseries->InsertColor(2, vtkColor3ub(255, 127, 127));

                auto plot = vtkSmartPointer<vtkPlotStacked>::New();
                plot->SetInputData(table);
                plot->SetInputArray(0, "x");
                plot->SetInputArray(1, "y1");
                plot->SetInputArray(2, "y2");
                plot->SetInputArray(3, "y3");
                plot->SetColorSeries(colorseries);

//                vtkSmartPointer<vtkChartXY> chart = vtkChartXY::SafeDownCast(m_contextViews[name - names.cbegin()][2 * c]->GetScene()->GetItem(0));
                vtkSmartPointer<vtkChartXY> chart = vtkChartXY::SafeDownCast(m_contextViews[name - names.cbegin()][c]->GetScene()->GetItem(0));
                chart->ClearPlots();
                chart->AddPlot(plot);

                auto pxarr = vtkSmartPointer<vtkFloatArray>::New();
                pxarr->SetNumberOfTuples(profilepoints);
                auto pyarr = vtkSmartPointer<vtkFloatArray>::New();
                pyarr->SetNumberOfTuples(profilepoints);

                if (m_volumes.find(*name) != m_volumes.end())
                {
                    ColumnVector profile = m_model->SampleProfile(vertex, profilepoints, profilespacing, *m_volumes[*name], m_transformation);

                    // TODO: Merge normalisation code with ShapeModel
                    double sxmin, sxmax, symin, symax, szmin, szmax;
                    m_model->GetShapeExtent(m_transformation, sxmin, sxmax, symin, symax, szmin, szmax);
                    double mean = m_model->GetMeanIntensity(*m_volumes[*name], m_normExclusionVolume.get(), sxmin, sxmax, symin, symax, szmin, szmax);

                    profile = m_model->NormaliseProfile(*name, profile, mean);

                    for (int p = 0; p < profilepoints; p++)
                    {
                        pxarr->SetTuple1(p, (p - profilepoints / 2.0) * profilespacing - displacement);
                        pxarr->SetName("x");
                        pyarr->SetTuple1(p, profile(p + 1));
                        pyarr->SetName("y");
                    }

                    auto ptable = vtkSmartPointer<vtkTable>::New();
                    ptable->AddColumn(pxarr);
                    ptable->AddColumn(pyarr);

                    auto pplot = vtkSmartPointer<vtkPlotLine>::New();
                    pplot->SetInputData(ptable, 0, 1);
                    pplot->SetColor(0, 0, 0, 255);

                    chart->AddPlot(pplot);
                }

                chart->SetTitle(*name + ": Component " + std::to_string(c) + ", theta("
                                + std::to_string(c) + ") = " + std::to_string(mixingcoefs[*name][c]));

//                m_contextViews[name - names.cbegin()][2 * c]->Render();
                m_contextViews[name - names.cbegin()][c]->Render();

//                auto mcxarr = vtkSmartPointer<vtkFloatArray>::New();
//                mcxarr->SetNumberOfTuples(profilepoints + steps);
//                mcxarr->SetName("x");
//                for (int p = 0; p < profilepoints + steps; p++)
//                    mcxarr->SetTuple1(p, (p - (profilepoints + steps) / 2.0) * profilespacing);

//                auto mctable = vtkSmartPointer<vtkTable>::New();
//                mctable->AddColumn(mcxarr);

//                auto meanssamples = pm->GetComponentMeansSamples()[*name][c];

//                // TODO: Remove columns if MCMC was not used
//                for (int i = 0; i < meanssamples.size(); i++)
//                {
//                    auto mcyarr = vtkSmartPointer<vtkFloatArray>::New();
//                    mcyarr->SetNumberOfTuples(profilepoints + steps);
//                    mcyarr->SetName((std::string("y") + std::to_string(i)).c_str());

//                    for (int p = 0; p < profilepoints + steps; p++)
//                        mcyarr->SetTuple1(p, meanssamples[i](p + 1));

//                    mctable->AddColumn(mcyarr);
//                }

//                vtkSmartPointer<vtkChartXY> mcchart = vtkChartXY::SafeDownCast(m_contextViews[name - names.cbegin()][2 * c + 1]->GetScene()->GetItem(0));
//                mcchart->ClearPlots();

//                for (int i = 0; i < meanssamples.size(); i++)
//                {
//                    auto mcline = mcchart->AddPlot(vtkChart::LINE);
//                    mcline->SetInputData(mctable);
//                    mcline->SetInputArray(0, "x");
//                    mcline->SetInputArray(1, (std::string("y") + std::to_string(i)).c_str());
//                }

//                mcchart->SetTitle("MCMC Samples of mean");

//                m_contextViews[name - names.cbegin()][2 * c + 1]->Render();
            }
        }
    }
    catch (std::exception &e)
    {
        QMessageBox::warning(this, "Error", (std::string("Caught exception while trying to plot profiles: ") + e.what()).c_str());
    }
    catch (Exception &e)
    {
        QMessageBox::warning(this, "Error", (std::string("Caught exception while trying to plot profiles: ") + e.what()).c_str());
    }
}

void PlotWindow::SetModel(boost::shared_ptr<const ShapeModel> model)
{
    // TODO: Error checking on dynamic casts!

    m_model = model;

    m_contextViews.clear();

    if (m_model)
    {
        boost::shared_ptr<const ProfileMixtures> vertexmodel =
                boost::dynamic_pointer_cast<const ProfileMixtures>(m_model->GetVertexModel(0));

        if (!vertexmodel)
        {
            QMessageBox::warning(this, "Error", "Dynamic cast of profile model to ProfileMixtures failed");
            return;
        }

        std::vector<std::string> names = m_model->GetModalityNames();
        int components = vertexmodel->GetNumberOfComponents();

    //    resize(2 * components * 300, names.size() * 300);
        resize(components * 300, names.size() * 300);

        for (int m = 0; m < names.size(); m++)
        {
            std::vector<vtkSmartPointer<vtkContextView> > cv;

    //        for (int c = 0; c < 2 * components; c++)
            for (int c = 0; c < components; c++)
            {
                auto chart = vtkSmartPointer<vtkChartXY>::New();
                chart->GetAxis(vtkAxis::BOTTOM)->SetTitle("Displacement (mm)");
                chart->GetAxis(vtkAxis::BOTTOM)->SetBehavior(1);
                float maxdisp = (m_model->GetProfilePoints() + vertexmodel->GetNumberOfSteps()) * m_model->GetProfileSpacing() / 2.0;
                chart->GetAxis(vtkAxis::BOTTOM)->SetRange(-maxdisp, maxdisp);
                chart->GetAxis(vtkAxis::LEFT)->SetTitle("Intensity");
                chart->GetAxis(vtkAxis::LEFT)->SetBehavior(0);

                auto qvtk = new QVTKWidget(this);
                auto contextview = vtkSmartPointer<vtkContextView>::New();
                contextview->SetRenderWindow(qvtk->GetRenderWindow());
                contextview->GetScene()->AddItem(chart);
                ui->grid->addWidget(qvtk, names.size() - 1 - m, c);
                cv.push_back(contextview);
            }

            m_contextViews.push_back(cv);
        }
    }
}

void PlotWindow::Clear()
{
    hide();

    QList<QVTKWidget *> widgets = findChildren<QVTKWidget *>();
    for (auto it = widgets.begin(); it != widgets.end(); ++it)
    {
        ui->grid->removeWidget(*it);
        delete *it;
    }

    m_model.reset();
    m_contextViews.clear();
    m_volumes.clear();
    m_normExclusionVolume.reset();
    m_transformation.reset();
}
