#!/usr/bin/env fslpython

import sys
import vtk
import numpy as np
import vtk.util.numpy_support as vtknp

if len(sys.argv) != 4:
    print("Usage: shrink_mesh.py [input mesh] [output mesh] [amount in mm]");
    sys.exit()

reader = vtk.vtkPolyDataReader()
reader.SetFileName(sys.argv[1])
reader.Update()

pd = reader.GetOutput()
filt = vtk.vtkPolyDataNormals()
filt.SetInputData(pd)
filt.SetComputeCellNormals(0)
filt.SetComputePointNormals(1)
filt.ConsistencyOff()
filt.SplittingOff()
filt.AutoOrientNormalsOn()
filt.Update()

newpoints = vtknp.vtk_to_numpy(pd.GetPoints().GetData()) \
        - float(sys.argv[3]) * vtknp.vtk_to_numpy(filt.GetOutput().GetPointData().GetNormals())
pd.GetPoints().SetData(vtknp.numpy_to_vtk(newpoints, deep = True))

writer = vtk.vtkPolyDataWriter()
writer.SetInputData(pd)
writer.SetFileName(sys.argv[2])
writer.SetFileTypeToASCII()
writer.Write()

