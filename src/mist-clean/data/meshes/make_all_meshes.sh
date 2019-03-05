#!/bin/bash

set -e -u

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm left_pallidum 6 1
../python/make_mesh.py left_pallidum.nii.gz left_pallidum 50 0 30

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm right_pallidum 17 1
../python/make_mesh.py right_pallidum.nii.gz right_pallidum 50 0 30

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm left_putamen 5 1
../python/make_mesh.py left_putamen.nii.gz left_putamen 50 0 30

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm right_putamen 16 1
../python/make_mesh.py right_putamen.nii.gz right_putamen 50 0 30

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm left_caudate 4 1
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm left_accumbens 10 1
fslmaths left_caudate -add left_accumbens left_caudate_accumbens
../python/make_mesh.py left_caudate_accumbens.nii.gz left_caudate_accumbens 50 1 100

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm right_caudate 15 1
fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm right_accumbens 20 1
fslmaths right_caudate -add right_accumbens right_caudate_accumbens
../python/make_mesh.py right_caudate_accumbens.nii.gz right_caudate_accumbens 50 1 100

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm left_thalamus 3 1
../python/make_mesh.py left_thalamus.nii.gz left_thalamus 50 0 100

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm right_thalamus 14 1
../python/make_mesh.py right_thalamus.nii.gz right_thalamus 50 0 100

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm left_hippocampus 8 1
../../python/make_mesh.py left_hippocampus.nii.gz left_hippocampus_large 40 1.0 100
../../python/shrink_mesh.py left_hippocampus_large_deformed_mesh.vtk left_hippocampus.vtk 1.0
cp left_hippocampus.vtk left_hippocampus.mim
../../bin/mist_mesh_utils gifti left_hippocampus.mim left_hippocampus.gii $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm right_hippocampus 18 1
../../python/make_mesh.py right_hippocampus.nii.gz right_hippocampus_large 40 1.0 100
../../python/shrink_mesh.py right_hippocampus_large_deformed_mesh.vtk right_hippocampus.vtk 1.0
cp right_hippocampus.vtk right_hippocampus.mim
../../bin/mist_mesh_utils gifti right_hippocampus.mim right_hippocampus.gii $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm left_amygdala 9 1
../../python/make_mesh.py left_amygdala.nii.gz left_amygdala 25 0 100
cp left_amygdala_deformed_mesh.vtk left_amygdala.mim
../../bin/mist_mesh_utils gifti left_amygdala.mim left_amygdala.gii $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz

fslroi $FSLDIR/data/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm right_amygdala 19 1
../../python/make_mesh.py right_amygdala.nii.gz right_amygdala 25 0 100
cp right_amygdala_deformed_mesh.vtk right_amygdala.mim
../../bin/mist_mesh_utils gifti right_amygdala.mim right_amygdala.gii $FSLDIR/data/standard/MNI152_T1_1mm.nii.gz

