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
import io
import math
import nibabel.gifti as gii
import nibabel.nifti1 as nii
import numpy as np
import os
import subprocess
import sys
import vtk
import vtk.util.numpy_support as vtknp

fsldir = os.environ['FSLDIR']
sys.path.append(fsldir + '/python/mist')
import radiological

# There was a bug prior to VTK 6.3.0 in vtkImplicitPolyDataDistance that caused it to sometimes get the sign wrong
vtkversion = vtk.vtkVersion()
majorver = vtkversion.GetVTKMajorVersion()
minorver = vtkversion.GetVTKMinorVersion()

if majorver < 6 or (majorver == 6 and minorver < 3):
    raise Exception('Please use at least VTK version 6.3.0 (loaded version is ' + vtkversion.GetVTKVersion() + ')')

def unpackpolys(packed):
    unpacked = list()
    i = 0
    while i < packed.shape[0]:
        unpacked.append(packed[i + 1 : i + packed[i] + 1])
        i = i + packed[i] + 1
    
    return unpacked

def packpolys(unpacked):
    packed = np.ndarray([0], dtype = 'int64')
    for p in unpacked:
        packed = np.r_[packed, p.shape[0], p]

    return packed

def polydata_to_points_polys(polydata, get_polys = False):
    points = vtknp.vtk_to_numpy(polydata.GetPoints().GetData()).copy()
    
    if not get_polys:
        return points
    
    polys = unpackpolys(vtknp.vtk_to_numpy(polydata.GetPolys().GetData()).copy())

    return points, polys

def points_polys_to_polydata(points, polys):
    polydata = vtk.vtkPolyData()
    
    vtkpoints = vtk.vtkPoints()
    vtkpoints.SetData(vtknp.numpy_to_vtk(points.copy(), deep = 1))
    polydata.SetPoints(vtkpoints)
    
    vtkpolys = vtk.vtkCellArray()
    vtkpolys.SetCells(len(polys),
            vtknp.numpy_to_vtkIdTypeArray(packpolys(polys), deep = 1))
    polydata.SetPolys(vtkpolys)

    return polydata

def get_normals(polydata):
    normalfilter = vtk.vtkPolyDataNormals()
    normalfilter.SetInputData(polydata)
    normalfilter.ComputeCellNormalsOff()
    normalfilter.ComputePointNormalsOn()
    normalfilter.SplittingOff()
    normalfilter.ConsistencyOff()
    normalfilter.AutoOrientNormalsOn()
    normalfilter.Update()

    return vtknp.vtk_to_numpy(normalfilter.GetOutput().GetPointData().GetNormals()).copy()

def loadmesh(filename):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(filename)
    
    if reader.IsFilePolyData() == False:
        raise Exception('Not a VTK polydata file: ' + filename)
    
    reader.Update()
    polydata = reader.GetOutput()

    return polydata

def writemesh(filename, polydata):
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInputData(polydata)
    writer.Write()

def writegifti(filename, polydata, refniifn):
    refnii = radiological.load(refniifn)
    
    points, polys = polydata_to_points_polys(polydata, True)
    scaledpoints4 = np.hstack([points / refnii.header.get_zooms()[0 : 3], np.ones(points.shape[0])[:, None]])
    transformedpoints = refnii.affine.dot(scaledpoints4.transpose())[0 : 3, :].transpose()

    # Not setting xfm as nifti mat has already been applied - is this the right way to do it?
    giftipoints = gii.GiftiDataArray(data = transformedpoints, intent = 'NIFTI_INTENT_POINTSET',
                                     datatype = 'NIFTI_TYPE_FLOAT32', encoding = 'GIFTI_ENCODING_ASCII')
    giftipolys = gii.GiftiDataArray(data = polys, intent = 'NIFTI_INTENT_TRIANGLE', datatype = 'NIFTI_TYPE_INT32',
                                    encoding = 'GIFTI_ENCODING_ASCII')
 
    giftimesh = gii.GiftiImage(darrays = [giftipoints, giftipolys])
    giftimesh.to_filename(filename)

def meanmesh(meshes):
    # TODO: Check if all meshes have same polys
    _, polys = polydata_to_points_polys(meshes[0], True)

    allpoints = np.stack([polydata_to_points_polys(mesh) for mesh in meshes])
    meanpoints = allpoints.mean(axis = 0)

    return points_polys_to_polydata(meanpoints, polys)

def meshvolume(polydata):
    mp = vtk.vtkMassProperties()
    mp.SetInputData(polydata)
    
    return mp.GetVolume()

def register(inmesh, refmesh, allow_rotation, allow_scaling):
    if allow_scaling and not allow_rotation:
        raise Exception('Cannot compute scaling without rotating the input mesh')
    
    # See 'Procrustes analysis' on Wikipedia and Brian's thesis / Horn paper
    
    inpoints, inpolys = polydata_to_points_polys(inmesh, True)
    refpoints = polydata_to_points_polys(refmesh)

    if len(inpoints) != len(refpoints):
        raise Exception('The input and reference meshes do not have the same number of points')

    inmean = inpoints.mean(axis = 0)
    refmean = refpoints.mean(axis = 0)
    
    incentered = inpoints - inmean
    refcentered = refpoints - refmean

    outpoints = incentered

    if allow_rotation:
        M = incentered.transpose().dot(refcentered)
        U, S, Vt = np.linalg.svd(M)
        rotation = U.dot(Vt)

        outpoints = outpoints.dot(rotation)

        if allow_scaling:
            outpoints *= np.sqrt((refcentered ** 2).sum()) / np.sqrt((incentered ** 2).sum())
    
    outpoints += refmean

    return points_polys_to_polydata(outpoints, inpolys)


def point_distances(inmesh, refmesh):
    inpoints = polydata_to_points_polys(inmesh)
    refpoints = polydata_to_points_polys(refmesh)
    normals = get_normals(refmesh)

    dists = np.sum(normals * (inpoints - refpoints), axis = 1)
 
    return dists

def point_to_mesh_distances(inmesh, refmesh):
     distfunc = vtk.vtkImplicitPolyDataDistance()
     distfunc.SetInput(refmesh)
     
     points = polydata_to_points_polys(inmesh)
     dists = np.zeros(len(points))
     for i, p in enumerate(points):
         dists[i] = distfunc.FunctionValue(p[0], p[1], p[2])
  
     return dists

def non_overlapping_segmentation(meshfiles, reffile, resolve):
    refnii = radiological.load(reffile)
    meshes = [loadmesh(mf) for mf in meshfiles]

    allmasks = np.zeros(refnii.shape + (len(meshes), ))

    for i, mesh in enumerate(meshes):
        imd = vtk.vtkImageData()
        imd.SetExtent(0, refnii.shape[0] - 1, 0, refnii.shape[1] - 1, 0, refnii.shape[2] - 1)
        imd.SetSpacing(*refnii.get_header().get_zooms())

        filt = vtk.vtkSelectEnclosedPoints()
        filt.SetInputData(imd)
        filt.SetSurfaceData(mesh)
        filt.SetTolerance(0.00001)
        filt.Update()

        mask = np.reshape(vtknp.vtk_to_numpy(filt.GetOutput().GetPointData().GetArray('SelectedPoints')).copy(), refnii.shape, order = 'F')

        allmasks[..., i] = mask

    labels = (allmasks * np.arange(1, 1 + len(meshes))).sum(axis = 3)
    contested = np.transpose(np.nonzero(allmasks.sum(axis = 3) > 1))

    if resolve:
        ipdds = list()
        for mesh in meshes:
            ipdd = vtk.vtkImplicitPolyDataDistance()
            ipdd.SetInput(mesh)
            ipdds.append(ipdd)

        for voxel in contested:
            # NB: This also includes the meshes that do not include the point, but those meshes will have higher distances
            #     and will not influence the result
            dists = [ipdd.FunctionValue(*(voxel * refnii.get_header().get_zooms())) for ipdd in ipdds]
            labels[voxel[0], voxel[1], voxel[2]] = np.argmin(np.array(dists)) + 1
    else:
        for voxel in contested:
            labels[voxel[0], voxel[1], voxel[2]] = 0

    result = nii.Nifti1Image(labels, refnii.affine)
    result.set_qform(refnii.affine)

    return result

def affine_transform(polydata, affine):
    points, polys = polydata_to_points_polys(polydata, True)
    points4 = np.hstack([points, np.ones(points.shape[0])[:, None]])
    transformed_points = affine.dot(points4.transpose())[0 : 3, :].transpose()
    
    return points_polys_to_polydata(transformed_points, polys)

def warp(polydata, warpfn, inniifn, refniifn):
    inpoints, polys = polydata_to_points_polys(polydata, True)
    
    innii = radiological.load(inniifn)
    inscaled4 = np.hstack([inpoints / innii.header.get_zooms()[0 : 3], np.ones(inpoints.shape[0])[:, None]])
    inworld = innii.affine.dot(inscaled4.transpose())[0 : 3, :].transpose()

    # TODO: Factor out call to img2stdcoord - but note the difference between img2std and std2img
    i2s_in = io.BytesIO()
    np.savetxt(i2s_in, inworld)
    i2s_cmd = [fsldir + '/bin/img2stdcoord', '-mm', '-img', inniifn, '-std', refniifn, '-warp', warpfn]
    i2s_cmd.append('-')

    i2s = subprocess.Popen(i2s_cmd, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
    i2s_out, _ = i2s.communicate(i2s_in.getvalue())
    xfmdpoints = np.loadtxt(io.StringIO(i2s_out))

    refnii = radiological.load(refniifn)
    outscaled4 = np.linalg.inv(refnii.affine).dot(np.hstack([xfmdpoints, np.ones(xfmdpoints.shape[0])[:, None]]).transpose()).transpose()
    outpoints = outscaled4[:, 0 : 3] * refnii.header.get_zooms()[0 : 3]
    
    return points_polys_to_polydata(outpoints, polys)

