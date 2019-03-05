#!/usr/bin/env fslpython

import sys
import subprocess
import vtk
import vtk.util.numpy_support as vtknp
import nibabel.nifti1 as nii
import scipy.ndimage.morphology

print('NB: make_mesh.py assumes radiological orientation!')

infile = sys.argv[1]
outroot = sys.argv[2]
threshold = float(sys.argv[3])
distweight = float(sys.argv[4])
iterations = int(sys.argv[5])

img = nii.Nifti1Image.load(infile)
zooms = img.get_header().get_zooms()
vol = img.get_data()

if distweight == 0.0:
    print('Not using distance weighting')
else:
    print('Using distance weighting')

distmap = scipy.ndimage.morphology.distance_transform_edt(vol < 95)
vol *= distweight * distmap ** 2 / 100 + 1

imd = vtk.vtkImageData()
imd.SetDimensions(vol.shape[0], vol.shape[1], vol.shape[2])
imd.SetSpacing(zooms[0], zooms[1], zooms[2])
imd.GetPointData().SetScalars(vtknp.numpy_to_vtk(vol.flatten('F'), deep = 1))

mcf = vtk.vtkImageMarchingCubes()
mcf.SetInputData(imd)
mcf.SetValue(0, threshold)
mcf.Update()

writer = vtk.vtkPolyDataWriter()
writer.SetInputData(mcf.GetOutput())
writer.SetFileName(outroot + '_orig_mesh.vtk')
writer.SetFileTypeToASCII()
writer.Write()

ndt = scipy.ndimage.morphology.distance_transform_edt(vol > threshold)
ndtimg = nii.Nifti1Image(ndt, img.get_affine())
ndtimg.to_filename(outroot + '_distmap')

subprocess.call("run_mesh_utils --doDeformSurface -i {0} -m {1} -o {2} --w_im=-0.2 --w_tan=1.0 --w_tri=1.0 --w_norm=0.0 -d {3} -t 0.2".format(
        outroot + '_distmap', outroot + '_orig_mesh.vtk', outroot + '_deformed_mesh.vtk', iterations), shell = True)

