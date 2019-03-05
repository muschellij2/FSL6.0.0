#!/usr/bin/env fslpython
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

import csv
import errno
import os
import os.path
import re
import shutil
import subprocess
import sys
import vtk

conffile = 'mist_filenames'
fsldir = os.getenv('FSLDIR')

fsldir = os.getenv('FSLDIR')
scriptdir = fsldir + '/python/mist'
sys.path.append(scriptdir)
import autosetup
import meshutils

def base_name(filename):
    for ext in ['.nii.gz', '.nii', '.img', '.hdr']:
        if filename.lower().endswith(ext):
            filename = filename[: -len(ext)]
            break

    return filename

def read_config():
    mtypes = dict()
    scans = dict()
    vsize = dict()
    alternate_extracted = None
    alternate_affine = None
    alternate_warp = None

    with open(conffile, 'rt') as fnf:
        for ln in csv.reader(fnf):
            if len(ln) > 0:
                name = ln[0]

                if name.lower() == 'alternate_extracted':
                    alternate_extracted = base_name(ln[1])
                elif name.lower() == 'alternate_affine':
                    alternate_affine = ln[1]
                elif name.lower() == 'alternate_warp':
                    alternate_warp = base_name(ln[1])
                else:
                    mtypes[name] = ln[1]
                    scans[name] = base_name(ln[2])
                    vsize[name] = float(ln[3])

    mtypes['__GMPVE__'] = 'GMPVE'
    scans['__GMPVE__'] = 'mist_t1_brain_pve_1'
    vsize['__GMPVE__'] = vsize['T1']

    return mtypes, scans, vsize, alternate_extracted, alternate_affine, alternate_warp

def get_number_of_points(meshfile):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(meshfile)
    reader.Update()

    return reader.GetOutput().GetNumberOfPoints()

def get_blocks(points, workers):
    blocks = [[] for i in range(workers)]

    for i in range(points):
        blocks[i % workers].append(i)

    return blocks

def read_directories(dirfile):
    with open(dirfile, 'r') as df:
        dirs = [d.rstrip() for d in df.readlines()]

    return dirs

def write_directories(outfile, dirs):
    with open(outfile, 'w') as df:
        df.writelines([d + '\n' for d in dirs])

def do_preproc(sdir):
    mtypes, scans, vsize, alternate_extracted, alternate_affine, alternate_warp = read_config()

    # NB: The image with *name* T1 is used for registration (not type, as that may not be unique)
    t1 = sdir + '/' + scans['T1']

    if alternate_extracted is not None:
        subprocess.check_call(['imcp', sdir + '/' + alternate_extracted, sdir + '/mist_t1_brain'])
    else:
        subprocess.check_call(['bet', t1, sdir + '/mist_t1_brain', '-R'])

    if alternate_affine is not None:
        shutil.copy(sdir + '/' + alternate_affine, sdir + '/mist_t1_brain_to_mni.mat')
        subprocess.check_call(['flirt', '-in', sdir + '/mist_t1_brain', '-ref', fsldir + '/data/standard/MNI152_T1_2mm_brain',
                               '-out', sdir + '/mist_t1_brain_to_mni_affine',
                               '-applyxfm', '-init', sdir + '/mist_t1_brain_to_mni.mat'])
    else:
        subprocess.check_call(['flirt', '-in', sdir + '/mist_t1_brain', '-ref', fsldir + '/data/standard/MNI152_T1_2mm_brain',
                               '-out', sdir + '/mist_t1_brain_to_mni_affine', '-omat', sdir + '/mist_t1_brain_to_mni.mat',
                               '-dof', '12'])

    if alternate_warp is not None:
        subprocess.check_call(['imcp', sdir + '/' + alternate_warp, sdir + '/mist_t1_to_mni_warp'])
        subprocess.check_call(['applywarp', '-i', t1, '-w', sdir + '/mist_t1_to_mni_warp',
                               '-r', fsldir + '/data/standard/MNI152_T1_2mm', '-o', sdir + '/mist_t1_to_mni_nonlin'])
    else:
        subprocess.check_call(['fnirt', '--in=' + t1, '--ref=' + fsldir + '/data/standard/MNI152_T1_2mm',
                               '--aff=' + sdir + '/mist_t1_brain_to_mni.mat', '--iout=' + sdir + '/mist_t1_to_mni_nonlin',
                               '--cout=' + sdir + '/mist_t1_to_mni_warp', '--config=' + fsldir + '/etc/flirtsch/T1_2_MNI152_2mm.cnf'])
    
    subprocess.check_call(['invwarp', '-w', sdir + '/mist_t1_to_mni_warp', '-r', fsldir + '/data/standard/MNI152_T1_2mm',
                           '-o', sdir + '/mist_t1_to_mni_warp_inv'])
    subprocess.check_call(['fast', '-R', '0.0', '-H', '0.0', sdir + '/mist_t1_brain'])

def do_autosetup(structures):
    mtypes, scans, vsize, _, _, _ = read_config()
    training_dirs = read_directories('mist_training_subjects')

    lcstructs = [s.lower() for s in structures]

    cmds = autosetup.generateall(training_dirs, mtypes, scans, 'mist_t1_brain_pve_0', 'mist_t1_to_mni_warp',
                                 'mist_t1_to_mni_warp_inv', vsize, 'mist_out', lcstructs if len(lcstructs) > 0 else None)
    
    with open('mist_out/training_commands', 'w') as f:
        f.writelines(ln + '\n' for ln in cmds)

def do_training(workers, workerid):
    with open('mist_out/training_commands', 'r') as f:
        for i, ln in enumerate(f):
            cmd = ln.rstrip()

            mesh = re.search(r"--shape=(\S*)", cmd).group(1)
            myblock = get_blocks(get_number_of_points(mesh), workers)[workerid]
            
            if len(myblock) > 0:
                workdir = "mist_out/work_" + str(i)
                mycmd = re.sub(r"--out=\S*", 
                                "--out=" + workdir + "/parallel --vertices=" + ",".join(str(v) for v in myblock), cmd)
            
                try:
                    os.mkdir(workdir)
                except OSError as e:
                    if e.errno != errno.EEXIST:
                        raise

                subprocess.check_call(mycmd, shell=True)

def do_merge():
    with open('mist_out/training_commands', 'r') as f:
        for i, ln in enumerate(f):
            cmd = ln.rstrip()
            
            workdir = "mist_out/work_" + str(i)
            mergecmd = re.sub(r" train ", " train --loadvertices=" + workdir + "/parallel ", cmd)

            subprocess.check_call(mergecmd, shell=True)

def do_fit(sdir):
    mtypes, scans, vsize, _, _, _ = read_config()

    meshes = []

    with open('mist_out/training_commands', 'r') as f:
        for ln in f:
            model = re.search(r"--out=(\S*)", ln).group(1)
            modelbase = re.match(r'model_(.*?)\.', os.path.basename(model)).group(1)
            modalities = re.search(r"--modalitynames=(\S*)", ln).group(1).split(',')

            subprocess.check_call(['convert_xfm', '-omat', sdir + '/mist_t1_brain_to_mni_inv.mat',
                                   '-inverse', sdir + '/mist_t1_brain_to_mni.mat'])
            
            subprocess.check_call(['mist', 'fit', '--model=' + model,
                                   '--modalitynames=' + ','.join(modalities),
                                   '--modalityimages=' + ','.join(sdir + '/' + scans[m] for m in modalities),
                                   '--warp=' + sdir + '/mist_t1_to_mni_warp',
                                   '--normexclusion=' + sdir + '/mist_t1_brain_pve_0',
                                   '--outregmat=' + sdir + '/mist_t1_brain_to_mni_inv.mat',
                                   '--outregref=' + fsldir + '/data/standard/MNI152_T1_1mm_brain',
                                   '--outbase=' + sdir + '/mist_' + modelbase])

            meshutils.writegifti(sdir + '/mist_' + modelbase + '_shape.gii',
                                 meshutils.loadmesh(sdir + '/mist_' + modelbase + '_shape.mim'),
                                 sdir + '/mist_t1_brain.nii.gz')

            meshutils.writegifti(sdir + '/mist_' + modelbase + '_shape_reg.gii',
                                 meshutils.loadmesh(sdir + '/mist_' + modelbase + '_shape_reg.mim'),
                                 fsldir + '/data/standard/MNI152_T1_1mm_brain.nii.gz')
            
            meshes.append(sdir + '/mist_' + modelbase + '_shape.mim')
    
    meshutils.non_overlapping_segmentation(meshes, sdir + '/mist_t1_brain.nii.gz', False).to_filename(sdir + '/mist_nonoverlapping.nii.gz')

def do_tables():
    all_dirs = read_directories('mist_subjects')

    with open('mist_out/training_commands', 'r') as f:
        for ln in f:
            model = re.search(r"--out=(\S*)", ln).group(1)
            modelbase = re.match(r'model_(.*?)\.', os.path.basename(model)).group(1)

            reference = meshutils.loadmesh(fsldir + '/data/mist/meshes/' + modelbase + '.mim')
            
            rowsnative = list()
            rowsmni = list()

            for sdir in all_dirs:
                # 6 DOF registration for native meshes
                nativemesh = meshutils.loadmesh(sdir + '/mist_' + modelbase + '_shape.mim')
                rowsnative.append(meshutils.point_distances(meshutils.register(nativemesh, reference, True, False), reference))
                
                # 3 DOF registration for MNI meshes
                mnimesh = meshutils.loadmesh(sdir + '/mist_' + modelbase + '_shape_reg.mim')
                rowsmni.append(meshutils.point_distances(meshutils.register(mnimesh, reference, False, False), reference))
            
            with open('mist_out/' + modelbase + '_distances_native.csv', 'w') as f:
                cw = csv.writer(f)
                cw.writerows(rowsnative)
            
            with open('mist_out/' + modelbase + '_distances_mni.csv', 'w') as f:
                cw = csv.writer(f)
                cw.writerows(rowsmni)
                
if __name__ == '__main__':
    verb = sys.argv[1];

    if verb == 'preproc':
        do_preproc(sys.argv[2])
    elif verb == 'autosetup':
        do_autosetup(sys.argv[2 :])
    elif verb == 'train':
        do_training(int(sys.argv[2]), int(sys.argv[3]))
    elif verb == 'merge':
        do_merge()
    elif verb == 'fit':
        do_fit(sys.argv[2])
    elif verb == 'maketables':
        do_tables()
    else:
        raise Exception('Invalid verb')

