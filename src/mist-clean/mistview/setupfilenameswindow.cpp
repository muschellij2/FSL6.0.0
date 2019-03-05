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

#include "setupfilenameswindow.h"
#include "ui_setupfilenameswindow.h"
#include <QSettings>
#include <cstdlib>
#include <string>

// NB: QT will not show the dialog if it doesn't fit on the screen!

SetupFilenamesWindow::SetupFilenamesWindow(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SetupFilenamesWindow)
{
    ui->setupUi(this);
    setWindowFlags(Qt::WindowFlags(Qt::Dialog | Qt::CustomizeWindowHint | Qt::WindowTitleHint |
                                   Qt::WindowCloseButtonHint));

    std::string lastselected = LoadConfigurations();

    for (const auto &t : m_configurations)
    {
        ui->configurationComboBox->addItem(t.first.c_str());
        if (t.first == lastselected)
            ui->configurationComboBox->setCurrentIndex(ui->configurationComboBox->count() - 1);
    }
}

SetupFilenamesWindow::~SetupFilenamesWindow()
{
    delete ui;
}

std::string SetupFilenamesWindow::LoadConfigurations()
{
    m_configurations.clear();

    QSettings settings;
    int configurations = settings.beginReadArray("NewConfigurations");
    for (int i = 0; i < configurations; i++)
    {
        settings.setArrayIndex(i);
        Configuration config;
        config.modelFilename = settings.value("Model").toString().toStdString();
        config.transformFilename = settings.value("Transformation").toString().toStdString();
        config.displacementsBasename = settings.value("Displacements").toString().toStdString();
        config.normExclusionFilename = settings.value("NormalisationExclusionMask").toString().toStdString();
        config.outputFilename = settings.value("Output").toString().toStdString();
        config.extraShapeFilename = settings.value("ExtraShape").toString().toStdString();
        config.overlayFilename = settings.value("Overlay").toString().toStdString();

        int modalities = settings.beginReadArray("Modalities");
        for (int j = 0; j < modalities; j++)
        {
            settings.setArrayIndex(j);
            config.volumeFilenames[settings.value("ModalityName").toString().toStdString()] =
                    settings.value("ModalityVolume").toString().toStdString();
        }

        settings.endArray();
        m_configurations[settings.value("Name").toString().toStdString()] = config;
    }

    settings.endArray();

    int shapes = settings.beginReadArray("MultiShape");
    ui->multiShapeTable->setRowCount(shapes);

    for (int i = 0; i < shapes; i++)
    {
        settings.setArrayIndex(i);
        ui->multiShapeTable->setItem(i, 0, new QTableWidgetItem(settings.value("Filename").toString()));
        ui->multiShapeTable->setItem(i, 1, new QTableWidgetItem(settings.value("DisplacementsBasename").toString()));
        ui->multiShapeTable->setItem(i, 2, new QTableWidgetItem(std::to_string(settings.value("Red").toFloat()).c_str()));
        ui->multiShapeTable->setItem(i, 3, new QTableWidgetItem(std::to_string(settings.value("Green").toFloat()).c_str()));
        ui->multiShapeTable->setItem(i, 4, new QTableWidgetItem(std::to_string(settings.value("Blue").toFloat()).c_str()));
    }

    settings.endArray();

    return settings.value("LastSelected").toString().toStdString();
}

void SetupFilenamesWindow::SaveConfigurations(const std::string &selected)
{
    QSettings settings;

    settings.beginWriteArray("NewConfigurations");
    settings.remove("");
    int i = 0;
    for (const auto &t : m_configurations)
    {
        settings.setArrayIndex(i++);
        settings.setValue("Name", t.first.c_str());

        const Configuration &config = t.second;
        settings.setValue("Model", config.modelFilename.c_str());
        settings.setValue("Transformation", config.transformFilename.c_str());
        settings.setValue("Displacements", config.displacementsBasename.c_str());
        settings.setValue("NormalisationExclusionMask", config.normExclusionFilename.c_str());
        settings.setValue("Output", config.outputFilename.c_str());
        settings.setValue("ExtraShape", config.extraShapeFilename.c_str());
        settings.setValue("Overlay", config.overlayFilename.c_str());

        settings.beginWriteArray("Modalities");
        int j = 0;
        for (const auto &s : config.volumeFilenames)
        {
            settings.setArrayIndex(j++);
            settings.setValue("ModalityName", s.first.c_str());
            settings.setValue("ModalityVolume", s.second.c_str());
        }

        settings.endArray();
    }

    settings.endArray();

    settings.beginWriteArray("MultiShape");
    for (int i = 0; i < ui->multiShapeTable->rowCount(); i++)
    {
        settings.setArrayIndex(i);
        settings.setValue("Filename", ui->multiShapeTable->item(i, 0)->text());
        settings.setValue("DisplacementsBasename", ui->multiShapeTable->item(i, 1)->text());
        settings.setValue("Red", std::atof(ui->multiShapeTable->item(i, 2)->text().toStdString().c_str()));
        settings.setValue("Green", std::atof(ui->multiShapeTable->item(i, 3)->text().toStdString().c_str()));
        settings.setValue("Blue", std::atof(ui->multiShapeTable->item(i, 4)->text().toStdString().c_str()));
    }
    settings.endArray();

    settings.setValue("LastSelected", selected.c_str());
    settings.sync();
}

void SetupFilenamesWindow::ShowConfiguration(const std::string &name)
{
    if (!name.empty())
    {
        const Configuration &config = m_configurations.find(name)->second;

        ui->modelFilenameLineEdit->setText(config.modelFilename.c_str());
        ui->transformationFilenameLineEdit->setText(config.transformFilename.c_str());
        ui->displacementsBasenameLineEdit->setText(config.displacementsBasename.c_str());
        ui->normExclusionFilenameLineEdit->setText(config.normExclusionFilename.c_str());
        ui->outputFilenameLineEdit->setText(config.outputFilename.c_str());
        ui->extraShapeFilenameLineEdit->setText(config.extraShapeFilename.c_str());
        ui->overlayFilenameLineEdit->setText(config.overlayFilename.c_str());

        ui->volumeFilenameTable->setRowCount(config.volumeFilenames.size());
        int i = 0;
        for (const auto &t : config.volumeFilenames)
        {
            ui->volumeFilenameTable->setItem(i, 0, new QTableWidgetItem(t.first.c_str()));
            ui->volumeFilenameTable->setItem(i, 1, new QTableWidgetItem(t.second.c_str()));
            i++;
        }
    }
}

void SetupFilenamesWindow::EmptyConfiguration()
{
    ui->modelFilenameLineEdit->clear();
    ui->transformationFilenameLineEdit->clear();
    ui->displacementsBasenameLineEdit->clear();
    ui->normExclusionFilenameLineEdit->clear();
    ui->outputFilenameLineEdit->clear();
    ui->extraShapeFilenameLineEdit->clear();
    ui->overlayFilenameLineEdit->clear();

    ui->volumeFilenameTable->setRowCount(0);
}

void SetupFilenamesWindow::on_configurationComboBox_currentIndexChanged(const QString &text)
{
    ShowConfiguration(text.toStdString());
}

void SetupFilenamesWindow::on_deleteButton_clicked()
{
    auto match = m_configurations.find(ui->configurationComboBox->currentText().toStdString());
    if (match != m_configurations.end())
    {
        m_configurations.erase(match);
        ui->configurationComboBox->removeItem(ui->configurationComboBox->currentIndex());
    }
    else if (ui->configurationComboBox->count())
    {
        ui->configurationComboBox->setCurrentIndex(0);
    }
    else
    {
        ui->configurationComboBox->clearEditText();
        EmptyConfiguration();
    }
}

void SetupFilenamesWindow::accept()
{
    // This ensures that any open editors on the QTableWidget / QComboBox(?) are closed ...
    ui->buttonBox->setFocus();

    std::string selected = ui->configurationComboBox->currentText().toStdString();
    if (!selected.empty())
    {
        Configuration config;
        config.modelFilename = ui->modelFilenameLineEdit->text().toStdString();
        config.transformFilename = ui->transformationFilenameLineEdit->text().toStdString();
        config.displacementsBasename = ui->displacementsBasenameLineEdit->text().toStdString();
        config.normExclusionFilename = ui->normExclusionFilenameLineEdit->text().toStdString();
        config.outputFilename = ui->outputFilenameLineEdit->text().toStdString();
        config.extraShapeFilename = ui->extraShapeFilenameLineEdit->text().toStdString();
        config.overlayFilename = ui->overlayFilenameLineEdit->text().toStdString();

        for (int i = 0; i < ui->volumeFilenameTable->rowCount(); i++)
            config.volumeFilenames[ui->volumeFilenameTable->item(i, 0)->text().toStdString()] =
                    ui->volumeFilenameTable->item(i, 1)->text().toStdString();

        m_configurations[ui->configurationComboBox->currentText().toStdString()] = config;

        SaveConfigurations(selected);
    }
    else if (m_configurations.empty())
        SaveConfigurations(std::string());
    else
        SaveConfigurations(m_configurations.cbegin()->first);

    QDialog::accept();
}

int SetupFilenamesWindow::GetNumberOfVolumes() const
{
    return ui->volumeFilenameTable->rowCount();
}

void SetupFilenamesWindow::GetVolumeFilename(int i, std::string& name, std::string& filename)
{
    if (i < ui->volumeFilenameTable->rowCount())
    {
        name = ui->volumeFilenameTable->item(i, 0)->text().toStdString();
        filename = ui->volumeFilenameTable->item(i, 1)->text().toStdString();
    }
    else
    {
        name = std::string();
        filename = std::string();
    }
}

std::string SetupFilenamesWindow::GetModelFilename() const
{
    return ui->modelFilenameLineEdit->text().toStdString();
}

std::string SetupFilenamesWindow::GetTransformFilename() const
{
    return ui->transformationFilenameLineEdit->text().toStdString();
}

std::string SetupFilenamesWindow::GetDisplacementsBasename() const
{
    return ui->displacementsBasenameLineEdit->text().toStdString();
}

std::string SetupFilenamesWindow::GetNormExclusionFilename() const
{
    return ui->normExclusionFilenameLineEdit->text().toStdString();
}

std::string SetupFilenamesWindow::GetOutputFilename() const
{
    return ui->outputFilenameLineEdit->text().toStdString();
}

std::string SetupFilenamesWindow::GetExtraShapeFilename() const
{
    return ui->extraShapeFilenameLineEdit->text().toStdString();
}

std::string SetupFilenamesWindow::GetOverlayFilename() const
{
    return ui->overlayFilenameLineEdit->text().toStdString();
}

int SetupFilenamesWindow::GetNumberOfShapes() const
{
    return ui->multiShapeTable->rowCount();
}

void SetupFilenamesWindow::GetShape(int i, std::string &filename, std::string &displacementsbasename, float &red, float &green, float &blue) const
{
    if (i < ui->multiShapeTable->rowCount())
    {
        filename = ui->multiShapeTable->item(i, 0)->text().toStdString();
        displacementsbasename = ui->multiShapeTable->item(i, 1)->text().toStdString();
        red = std::atof(ui->multiShapeTable->item(i, 2)->text().toStdString().c_str());
        green = std::atof(ui->multiShapeTable->item(i, 3)->text().toStdString().c_str());
        blue = std::atof(ui->multiShapeTable->item(i, 4)->text().toStdString().c_str());
    }
    else
    {
        filename = std::string();
        displacementsbasename = std::string();
        red = 0.0;
        green = 0.0;
        blue = 0.0;
    }
}

void SetupFilenamesWindow::on_addButton_clicked()
{
    int row = ui->volumeFilenameTable->rowCount();
    ui->volumeFilenameTable->setRowCount(row + 1);
    ui->volumeFilenameTable->setItem(row, 0, new QTableWidgetItem());
    ui->volumeFilenameTable->setItem(row, 1, new QTableWidgetItem());
}

void SetupFilenamesWindow::on_removeButton_clicked()
{
    int index = ui->volumeFilenameTable->currentRow();
    if (index >= 0)
        ui->volumeFilenameTable->removeRow(index);
}

void SetupFilenamesWindow::on_multiShapeRemoveButton_clicked()
{
    int index = ui->multiShapeTable->currentRow();
    if (index >= 0)
        ui->multiShapeTable->removeRow(index);
}

void SetupFilenamesWindow::on_multiShapeAddButton_clicked()
{
    int row = ui->multiShapeTable->rowCount();
    ui->multiShapeTable->setRowCount(row + 1);
    ui->multiShapeTable->setItem(row, 0, new QTableWidgetItem());
    ui->multiShapeTable->setItem(row, 1, new QTableWidgetItem());
    ui->multiShapeTable->setItem(row, 2, new QTableWidgetItem());
    ui->multiShapeTable->setItem(row, 3, new QTableWidgetItem());
    ui->multiShapeTable->setItem(row, 4, new QTableWidgetItem());
}
