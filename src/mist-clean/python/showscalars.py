#   Copyright (C) 2016 University of Oxford 
#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   
#   LICENCE
#   
#   FMRIB Software Library, Release 5.0 (c) 2012, The University of
#   Oxford (the "Software")
#   
#   The Software remains the property of the University of Oxford ("the
#   University").
#   
#   The Software is distributed "AS IS" under this Licence solely for
#   non-commercial use in the hope that it will be useful, but in order
#   that the University as a charitable foundation protects its assets for
#   the benefit of its educational and research purposes, the University
#   makes clear that no condition is made or to be implied, nor is any
#   warranty given or to be implied, as to the accuracy of the Software,
#   or that it will be suitable for any particular purpose or for use
#   under any specific conditions. Furthermore, the University disclaims
#   all responsibility for the use which is made of the Software. It
#   further disclaims any liability for the outcomes arising from using
#   the Software.
#   
#   The Licensee agrees to indemnify the University and hold the
#   University harmless from and against any and all claims, damages and
#   liabilities asserted by third parties (including claims for
#   negligence) which arise directly or indirectly from the use of the
#   Software or the sale of any products based on the Software.
#   
#   No part of the Software may be reproduced, modified, transmitted or
#   transferred in any form or by any means, electronic or mechanical,
#   without the express permission of the University. The permission of
#   the University is not required if the said reproduction, modification,
#   transmission or transference is done without financial return, the
#   conditions of this Licence are imposed upon the receiver of the
#   product, and all original and amended source code is included in any
#   transmitted product. You may be held legally responsible for any
#   copyright infringement that is caused or encouraged by your failure to
#   abide by these terms and conditions.
#   
#   You are not permitted under this Licence to use this Software
#   commercially. Use for which any financial return is received shall be
#   defined as commercial use, and includes (1) integration of all or part
#   of the source code or the Software into a product for sale or license
#   by or on behalf of Licensee to third parties or (2) use of the
#   Software or any derivative of it for research with the final aim of
#   developing software products for sale or license to a third party or
#   (3) use of the Software or any derivative of it for research with the
#   final aim of developing non-software products for sale or license to a
#   third party, or (4) use of the Software to provide any service to an
#   external organisation for which payment is received. If you are
#   interested in using the Software commercially, please contact Oxford
#   University Innovation ("OUI"), the technology transfer company of the
#   University, to negotiate a licence. Contact details are:
#   Innovation@innovation.ox.ac.uk quoting reference DE/9564.
import nibabel.gifti as gii
import numpy as np
import os
import scipy.stats
import sys
import vtk
import vtk.util.numpy_support as vtknp

fsldir = os.getenv('FSLDIR')
sys.path.append(fsldir + '/python/mist')
import meshutils

def readmesh(meshfile):
    ext = os.path.splitext(meshfile)[1]
    if ext == '.mim':
        # Just use points and polys, not normals etc
        origpd = meshutils.loadmesh(meshfile)
        points, polys = meshutils.polydata_to_points_polys(origpd, True)
        pdfl = meshutils.points_polys_to_polydata(points, polys)
    
        # Internal meshes use radiological coordinates
        xfm = vtk.vtkTransform()
        xfm.Scale(-1.0, 1.0, 1.0)
        xfm.Update()
        xfmf = vtk.vtkTransformPolyDataFilter()
        xfmf.SetTransform(xfm)
        xfmf.SetInputData(pdfl)
        xfmf.Update()
        pd = xfmf.GetOutput()
    elif ext == '.gii':
        mesh = gii.giftiio.read(meshfile)
        pd = vtk.vtkPolyData()
        pts = vtk.vtkPoints()
        pts.SetData(vtknp.numpy_to_vtk(mesh.darrays[0].data, deep = 1))
        pd.SetPoints(pts)
        polysarr = np.c_[3 * np.ones(mesh.darrays[1].data.shape[0], dtype = 'int')[:, None],
                np.array(mesh.darrays[1].data, dtype = 'int')]
        polys = vtk.vtkCellArray()
        polys.SetCells(mesh.darrays[1].data.shape[0], vtknp.numpy_to_vtkIdTypeArray(polysarr, deep = 1))
        pd.SetPolys(polys)
    else:
        raise Exception('Unrecognised extension ({0})'.format(ext))

    return pd

def readscalars(scalarfile):
    if os.path.splitext(scalarfile)[1] == '.csv':
        scalars = np.loadtxt(scalarfile, delimiter = ',')
    else:
        scalars = np.loadtxt(scalarfile)

    return scalars

def readpvals(posfilename, negfilename = None):
    # This produces thresholded 1-p values (and p-1 values for the negative contrast)
    alpha = 0.05
    if negfilename is not None:
        alpha /= 2.0
    
    pos = readscalars(posfilename)
    scalars = np.zeros(len(pos))
    scalars[pos < alpha] = pos[pos < alpha]

    if negfilename is not None:
        neg = readscalars(negfilename)
        scalars[neg < alpha] = -neg[neg <= alpha]

    return scalars

def p_lut():
    lut = vtk.vtkLookupTable()
    lut.Allocate(501, 501)
    lut.SetTableRange(0.0, 0.05)
    lut.SetHueRange(0.0, 0.0)
    lut.SetSaturationRange(0.0, 0.0)
    lut.SetValueRange(0.5, 0.5)
    lut.Build()
    
    for i in range(0, 500):
        lut.SetTableValue(500 - i, 1.0, i / 500.0, 0.0)

    return lut

def p_lut_twosided():
    lut = vtk.vtkLookupTable()
    lut.Allocate(501, 501)
    lut.SetTableRange(-0.025, 0.025)
    lut.SetHueRange(0.0, 0.0)
    lut.SetSaturationRange(0.0, 0.0)
    lut.SetValueRange(0.5, 0.5)
    lut.Build()
    
    for i in range(0, 250):
        lut.SetTableValue(500 - i, 1.0, i / 250.0, 0.0)
        lut.SetTableValue(i, i / 500.0, i / 250.0, 1.0)

    return lut

def split_p_lut(lut):
    n = (lut.GetNumberOfTableValues() - 1) // 2
    poslut = vtk.vtkLookupTable()
    neglut = vtk.vtkLookupTable()
    poslut.Allocate(n, n)
    neglut.Allocate(n, n)
    poslut.SetTableRange(0.0, 0.025)
    neglut.SetTableRange(0.0, 0.025)
    
    for i in range(n):
        poslut.SetTableValue(i, lut.GetTableValue(n + 1 + i))
        neglut.SetTableValue(i, lut.GetTableValue(n - i))
        
    return poslut, neglut

def standard_lut():
    lut = vtk.vtkLookupTable()
    lut.Allocate(2000, 2000)
    lut.SetTableRange(-1.0, 1.0)
    lut.SetHueRange(0.66, 0.0)
    lut.SetSaturationRange(1.0, 1.0)
    lut.SetValueRange(1.0, 1.0)
    lut.Build()
    
    return lut

def getactor(polydata, scalars, smin, smax, lut = None):
    points = polydata.GetNumberOfPoints()
    
    if scalars.size != points:
        raise Exception('Incorrect number of points (expected {0})'.format(points))
    
    polydata.GetPointData().SetScalars(vtknp.numpy_to_vtk(scalars.copy(), deep = 1))
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    mapper.SetScalarVisibility(1)
    mapper.SetScalarRange(smin, smax)
    if lut:
        mapper.SetLookupTable(lut)
    else:
        mapper.SetLookupTable(standard_lut())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    return actor

def create_colourbar(lut, title = None):
    colourbar = vtk.vtkScalarBarActor()
    colourbar.SetBarRatio(0.2)
    colourbar.SetWidth(0.15)
    colourbar.SetHeight(0.8)
    colourbar.SetNumberOfLabels(2)
    colourbar.SetTextPositionToPrecedeScalarBar()
    colourbar.SetTextPad(5)
    colourbar.SetLabelFormat('%.2f')
    colourbar.GetLabelTextProperty().ItalicOff()
    colourbar.SetLookupTable(lut)
    colourbar.SetMaximumNumberOfColors(lut.GetNumberOfTableValues())
    colourbar.SetTitle(title if title else '')
    colourbar.GetTitleTextProperty().ItalicOff()

    return colourbar

def display(actors, width = 600, height = 600, viewnormal = None, viewup = None, colourbartitle = None, pvals = False, savepng = None, twosided = False):
    renderer = vtk.vtkRenderer()
    for actor in actors:
        renderer.AddActor(actor)
    
    if pvals:
        lut = actors[0].GetMapper().GetLookupTable()

        if twosided:
            poslut, neglut = split_p_lut(lut)
        else:
            poslut = lut

        poscolourbar = create_colourbar(poslut, colourbartitle)
        poscolourbar.SetLabelFormat('%.3f')
        poscolourbar.GetTitleTextProperty().SetOpacity(0.0)
        poscolourbar.SetPosition(0.8, 0.1)
        renderer.AddActor(poscolourbar)

        if twosided:
            negcolourbar = create_colourbar(neglut, colourbartitle)
            negcolourbar.SetLabelFormat('%.3f')
            negcolourbar.DrawTickLabelsOff()
            negcolourbar.SetPosition(0.85, 0.1)
            renderer.AddActor(negcolourbar)

    else:
        lut = actors[0].GetMapper().GetLookupTable()
        colourbar = create_colourbar(lut, colourbartitle)
        colourbar.SetPosition(0.8, 0.1)
        
        renderer.AddActor(colourbar)

    if viewnormal:
        if not viewup:
            raise Exception('When using viewnormal, viewup needs to be set as well')
        
        camera = renderer.GetActiveCamera()
        renderer.ResetCamera()
        camera.SetPosition(np.array(camera.GetFocalPoint()) - np.array(viewnormal))
        camera.SetViewUp(viewup)

    window = vtk.vtkRenderWindow()
    window.SetSize(width, height)
    window.AddRenderer(renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(window)
    
#    axes = vtk.vtkAxesActor()
#    axes.SetXAxisLabelText('R')
#    axes.SetYAxisLabelText('A')
#    axes.SetZAxisLabelText('S')
#    axeswidget = vtk.vtkOrientationMarkerWidget()
#    axeswidget.SetOrientationMarker(axes)
#    axeswidget.SetInteractor(interactor)
#    axeswidget.SetViewport(0.0, 0.0, 0.3, 0.3 * width / height)
#    axeswidget.EnabledOn()
   
    renderer.ResetCamera()
    renderer.GetActiveCamera().ParallelProjectionOff()
    renderer.GetActiveCamera().SetWindowCenter(0.1, 0.0)
    renderer.GetActiveCamera().Zoom(1.3)
    
    if savepng:
        window.Render()
        imfilt = vtk.vtkWindowToImageFilter()
        imfilt.SetInput(window)
        imfilt.ReadFrontBufferOff()
        imfilt.Update()
        pngwriter = vtk.vtkPNGWriter()
        pngwriter.SetInputData(imfilt.GetOutput())
        pngwriter.SetFileName(savepng)
        pngwriter.Write()
    else:
        interactor.Start()


