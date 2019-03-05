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
import copy
import nibabel.nifti1 as nii
import numpy as np
import os
import subprocess

fsldir = os.environ['FSLDIR']

mist_bin = fsldir + '/bin/mist'
mask_dir = fsldir + '/data/mist/masks'
mesh_dir = fsldir + '/data/mist/meshes'

mask_putamen = mask_dir + '/putamen_thr75' 
mask_pallidum = mask_dir + '/pallidum_thr75' 
mask_caudate_accumbens = mask_dir + '/caudate_accumbens_thr75' 
mask_thalamus = mask_dir + '/thalamus_thr75' 
mask_hippocampus = mask_dir + '/hippocampus_thr75' 
mask_amygdala = mask_dir + '/amygdala_thr75' 
mask_wm = mask_dir + '/wm_thr75' 
mask_csf = mask_dir + '/csf_thr75' 
mask_red_nucleus = mask_dir + '/red_nuclei_eroded'
mask_substantia_nigra = mask_dir + '/substantia_nigra_eroded'
mask_subthalamic_nucleus = mask_dir + '/subthalamic_nuclei_eroded'


def get_mask_medians(images, masks, warpinv, dirs):
    # Get mean across subjects of median intensities within masks

    numbers = subprocess.check_output([
            mist_bin,
            'maskmedians',
            '--masks={0}'.format(','.join([masks[k] for k in sorted(masks)])),
            '--images={0}'.format(','.join([images[k] for k in sorted(images)])),
            '--warp={0}'.format(warpinv)]
            + dirs).split()

    ni = iter(numbers)

    result = {mk : {ik : 0.0 for ik in images} for mk in masks} 

    for d in dirs:
        for mk in sorted(masks):
            for ik in sorted(images):
                result[mk][ik] += float(next(ni)) / len(dirs)

    try:
        _ = next(ni)
    except StopIteration:
        return result

    raise Exception('MIST returned too many numbers')


def drop_modalities(scans, modalitytypes, droptypes):
    newscans = copy.deepcopy(scans)

    for mname, mtype in modalitytypes.items():
        if mtype in droptypes:
            newscans.pop(mname)

    return newscans

def get_all_medians(scans, warpinv, dirs):
    masks = {'putamen': mask_putamen,
             'pallidum': mask_pallidum,
             'caudate_accumbens': mask_caudate_accumbens,
             'thalamus': mask_thalamus,
             'hippocampus': mask_hippocampus,
             'amygdala': mask_amygdala,
             'csf': mask_csf,
             'wm': mask_wm,
             'red_nucleus': mask_red_nucleus,
             'substantia_nigra': mask_substantia_nigra,
             'subthalamic_nucleus': mask_subthalamic_nucleus}

    return get_mask_medians(scans, masks, warpinv, dirs)

def setup_putamen(dirs, modalitytypes, scans, warpinv, scalarvoxelsizes, medians):
    minvoxelsize = min(scalarvoxelsizes.values())

    meansspecs = []
    covsspecs = []

    for name in sorted(scans):
        resfrac = scalarvoxelsizes[name] / minvoxelsize
        mt = modalitytypes[name]

        if mt == 'T1' or mt == 'T2':
            # First two components are for pallidum / WM
            # Third component handles vessels (now based on CSF to make autosetup simpler)
            meansspecs.append('"exp(4.0;{0};{1}),exp(1.0;{0};{1}),exp(1.0;{0};{2})"'.format(
                medians['putamen'][name],
                (medians['pallidum'][name] + medians['wm'][name]) / 2.0,
                (medians['putamen'][name] + medians['csf'][name]) / 2.0))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.1 * resfrac * medians['putamen'][name]) ** 2))
        elif mt == 'FA':
            # First two componentents are for the case where there is white matter outside (lateral boundary)
            # Third component is for when there is no WM (medial boundary)
            meansspecs.append('"exp(4.0;0.1;0.7),exp(1.0;0.1;0.4),flat(0.1)"')
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.1 * resfrac) ** 2))
        elif mt == 'QSM':
            meansspecs.append('"exp(2.0;{0};{1}),step({0};{2}),flat({0})"'.format(
                medians['putamen'][name],
                medians['pallidum'][name],
                medians['wm'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.3 * resfrac * medians['putamen'][name]) ** 2))
        else:
            raise Exception('Modality type ' + mt + ' is not supported for putamen')

    return '--components=3 --priormeans={0} --priorcovcoefs={1}'.format(
            ','.join(meansspecs), ','.join(covsspecs))

def setup_pallidum(dirs, modalitytypes, scans, warpinv, scalarvoxelsizes, medians):
    minvoxelsize = min(scalarvoxelsizes.values())

    meansspecs = []
    covsspecs = []

    for name in sorted(scans):
        resfrac = scalarvoxelsizes[name] / minvoxelsize
        mt = modalitytypes[name]

        if mt == 'T1':
            # First component is for boundary with putamen, second is for boundary with WM
            meansspecs.append('"exp(3.0;{0};{1}),exp(3.0;{0};{2})"'.format(
                medians['pallidum'][name],
                medians['putamen'][name],
                medians['pallidum'][name] + 3.0 * (medians['wm'][name] - medians['pallidum'][name])))
            covsspecs.append('"flat({0}),flat({0})"'.format((0.1 * resfrac * medians['pallidum'][name]) ** 2))
        elif mt == 'T2':
            # First component is for boundary with putamen, second is for boundary with WM
            meansspecs.append('"exp(3.0;{0};{1}),exp(3.0;{0};{2})"'.format(
                medians['pallidum'][name],
                medians['putamen'][name],
                medians['wm'][name]))
            covsspecs.append('"flat({0}),flat({0})"'.format((0.1 * resfrac * medians['pallidum'][name]) ** 2))
        elif mt == 'FA':
            # First componentent is for the case where there is white matter outside (medial boundary)
            # Second component is for when there is no WM (lateral boundary)
            meansspecs.append('"step(0.2;0.6),flat(0.2)"')
            covsspecs.append('"flat({0}),flat({0})"'.format((0.1 * resfrac) ** 2))
        elif mt == 'QSM':
            meansspecs.append('"step({0};{1}),step({0};{2})"'.format(
                medians['pallidum'][name],
                medians['putamen'][name],
                medians['wm'][name]))
            covsspecs.append('"flat({0}),flat({0})"'.format((0.3 * resfrac * medians['pallidum'][name]) ** 2))
        else:
            raise Exception('Modality type ' + mt + ' is not supported for pallidum')

    return '--components=2 --priormeans={0} --priorcovcoefs={1}'.format(
            ','.join(meansspecs), ','.join(covsspecs))

def setup_caudate_accumbens(dirs, modalitytypes, scans, warpinv, scalarvoxelsizes, medians):
    minvoxelsize = min(scalarvoxelsizes.values())

    meansspecs = []
    covsspecs = []

    for name in sorted(scans):
        resfrac = scalarvoxelsizes[name] / minvoxelsize
        mt = modalitytypes[name]

        if mt == 'T1':
            meansspecs.append('"exp(1.0;{0};{1}),exp(4.0;{0};{2}),exp2(1.0;{0};{2};{3}),exp(1.0;{0};{3}),exp(4.0;{0};{4})"'.format(
                medians['caudate_accumbens'][name],
                medians['caudate_accumbens'][name] - 0.33 * (medians['caudate_accumbens'][name] - medians['csf'][name]),
                medians['caudate_accumbens'][name] - 1.5 * (medians['caudate_accumbens'][name] - medians['csf'][name]),
                medians['caudate_accumbens'][name] - 0.33 * (medians['caudate_accumbens'][name] - medians['wm'][name]),
                medians['caudate_accumbens'][name] - 1.0 * (medians['caudate_accumbens'][name] - medians['wm'][name])))
            covsspecs.append('"flat({0}),flat({0}),flat({0}),flat({0}),flat({0})"'.format((0.1 * resfrac * medians['caudate_accumbens'][name]) ** 2))
        elif mt == 'T2':
            meansspecs.append('"step({0};{1}),exp(4.0;{0};{1}),exp(1.0;{0};{1}),step({0};{2}),exp(1.0;{0};{2})"'.format(
                medians['caudate_accumbens'][name], medians['csf'][name], medians['wm'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0}),flat({0}),flat({0})"'.format((0.1 * resfrac * medians['caudate_accumbens'][name]) ** 2))
        else:
            raise Exception('Modality type ' + mt + ' is not supported for caudate nucleus + nucleus accumbens')

    return '--components=5 --priormeans={0} --priorcovcoefs={1}'.format(
            ','.join(meansspecs), ','.join(covsspecs))


def setup_thalamus(dirs, modalitytypes, scans, warpinv, scalarvoxelsizes, medians):
    minvoxelsize = min(scalarvoxelsizes.values())

    meansspecs = []
    covsspecs = []

    for name in sorted(scans):
        resfrac = scalarvoxelsizes[name] / minvoxelsize
        mt = modalitytypes[name]

        if mt == 'T1' or mt == 'T2':
            meansspecs.append('"exp(1.0;{0};{1}),exp(6.0;{0};{1}),exp(1.0;{0};{2})"'.format(
                medians['thalamus'][name],
                medians['csf'][name],
                medians['wm'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.1 * resfrac * medians['thalamus'][name]) ** 2))
        elif mt == 'FA':
            meansspecs.append('"exp(1.0;0.2;0.5),exp(6.0;0.2;0.5),flat(0.2)"')
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.1 * resfrac) ** 2))
        else:
            raise Exception('Modality type ' + mt + ' is not supported for pallidum')

    return '--components=3 --priormeans={0} --priorcovcoefs={1}'.format(
            ','.join(meansspecs), ','.join(covsspecs))

def setup_hippocampus(dirs, modalitytypes, scans, warpinv, scalarvoxelsizes, medians):
    minvoxelsize = min(scalarvoxelsizes.values())

    meansspecs = []
    covsspecs = []

    for name in sorted(scans):
        resfrac = scalarvoxelsizes[name] / minvoxelsize
        mt = modalitytypes[name]

        if mt == 'T1':
            meansspecs.append('"exp(1.0;{0};{6}),exp2(1.0;{0};{6};{3}),exp(1.0;{0};{5}),doubleexp(1.0;2.0;{0};{2};{5}),doubleexp(2.0;4.0;{0};{2};{5}),doubleexp2(1.0;2.0;{0};{2};{5};{2})"'.format(
                medians['hippocampus'][name],
                medians['wm'][name],
                (medians['wm'][name] + medians['hippocampus'][name]) / 2.0,
                (medians['wm'][name] + 2.0 * medians['hippocampus'][name]) / 3.0,
                (medians['csf'][name] + medians['hippocampus'][name]) / 2.0,
                (medians['csf'][name] + 3.0 * medians['hippocampus'][name]) / 4.0,
                (3.0 * medians['wm'][name] + medians['hippocampus'][name]) / 4.0))
            covsspecs.append('"flat({0}),flat({0}),flat({0}),flat({0}),flat({0}),flat({0})"'.format((0.03 * resfrac * medians['hippocampus'][name]) ** 2))
        elif mt == 'T2' or mt == 'WMS':
            meansspecs.append('"exp(1.0;{0};{6}),exp(3.0;{0};{6}),exp(1.0;{0};{5}),doubleexp(1.0;1.5;{0};{2};{5}),doubleexp(2.0;4.0;{0};{2};{5}),flat({0})"'.format(
                medians['hippocampus'][name],
                medians['wm'][name],
                (medians['wm'][name] + medians['hippocampus'][name]) / 2.0,
                (medians['wm'][name] + 2.0 * medians['hippocampus'][name]) / 3.0,
                (medians['csf'][name] + medians['hippocampus'][name]) / 2.0,
                (medians['csf'][name] + 3.0 * medians['hippocampus'][name]) / 4.0,
                (3.0 * medians['wm'][name] + medians['hippocampus'][name]) / 4.0))
            covsspecs.append('"flat({0}),flat({0}),flat({0}),flat({0}),flat({0}),flat({0})"'.format((0.03 * resfrac * medians['hippocampus'][name]) ** 2))
        elif mt == 'GMPVE':
            meansspecs.append('"step(1.0;0.5),flat(1.0),flat(1.0),flat(1.0),flat(1.0),flat(1.0)"')
            covsspecs.append('"flat(0.001),flat(0.001),flat(0.001),flat(0.001),flat(0.001),flat(0.001)"')
        else:
            raise Exception('Modality type ' + mt + ' is not supported for hippocampus')

    return '--components=6 --priormeans={0} --priorcovcoefs={1}'.format(
            ','.join(meansspecs), ','.join(covsspecs))

def setup_amygdala(dirs, modalitytypes, scans, warpinv, scalarvoxelsizes, medians):
    minvoxelsize = min(scalarvoxelsizes.values())

    meansspecs = []
    covsspecs = []

    for name in sorted(scans):
        resfrac = scalarvoxelsizes[name] / minvoxelsize
        mt = modalitytypes[name]

        if mt == 'T1':
            meansspecs.append('"exp(3.0;{0};{1}),exp(1.0;{0};{1}),exp(2.0;{0};{2})"'.format(
                medians['hippocampus'][name],
                (medians['wm'][name] + medians['hippocampus'][name]) / 2.0,
                (medians['csf'][name] + medians['hippocampus'][name]) / 2.0))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.03 * resfrac * medians['hippocampus'][name]) ** 2))
        elif mt == 'T2' or mt == 'WMS':
            meansspecs.append('"flat({0}),exp(2.0;{0};{1}),exp(2.0;{0};{2})"'.format(
                medians['hippocampus'][name],
                (medians['wm'][name] + medians['hippocampus'][name]) / 2.0,
                (medians['csf'][name] + medians['hippocampus'][name]) / 2.0))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.03 * resfrac * medians['hippocampus'][name]) ** 2))
        elif mt == 'GMPVE':
            meansspecs.append('"exp(2.0;1.0;0.5),flat(1.0),flat(1.0)"')
            covsspecs.append('"flat(0.001),flat(0.001),flat(0.001)"')
        else:
            raise Exception('Modality type ' + mt + ' is not supported for amygdala')

    return '--components=3 --priormeans={0} --priorcovcoefs={1}'.format(
            ','.join(meansspecs), ','.join(covsspecs))

def setup_red_nucleus(dirs, modalitytypes, scans, warpinv, scalarvoxelsizes, medians):
    minvoxelsize = min(scalarvoxelsizes.values())

    meansspecs = []
    covsspecs = []

    for name in sorted(scans):
        resfrac = scalarvoxelsizes[name] / minvoxelsize
        mt = modalitytypes[name]

        if mt == 'T1':
            meansspecs.append('"flat({0}),flat({0})"'.format(medians['red_nucleus'][name]))
            covsspecs.append('"flat({0}),flat({0})"'.format((0.1 * resfrac * medians['red_nucleus'][name]) ** 2))
        elif mt == 'T2':
            meansspecs.append('"exp(1.0;{0};{1}),exp(3.0;{0};{1})"'.format(medians['red_nucleus'][name], 1.33 * medians['red_nucleus'][name]))
            covsspecs.append('"flat({0}),flat({0})"'.format((0.1 * resfrac * medians['red_nucleus'][name]) ** 2))
        elif mt == 'WMS':
            meansspecs.append('"exp(1.0;{0};{1}),exp(3.0;{0};{1})"'.format(medians['red_nucleus'][name], 1.33 * medians['red_nucleus'][name]))
            covsspecs.append('"flat({0}),flat({0})"'.format((0.5 * resfrac * medians['red_nucleus'][name]) ** 2))
        elif mt == 'QSM':
            meansspecs.append('"exp(1.0;{0};{1}),exp(3.0;{0};{1})"'.format(medians['red_nucleus'][name], medians['wm'][name]))
            covsspecs.append('"flat({0}),flat({0})"'.format((0.3 * resfrac * medians['red_nucleus'][name]) ** 2))
        else:
            raise Exception('Modality type ' + mt + ' is not supported for red nucleus')

    return '--components=2 --priormeans={0} --priorcovcoefs={1}'.format(
            ','.join(meansspecs), ','.join(covsspecs))


def setup_substantia_nigra(dirs, modalitytypes, scans, warpinv, scalarvoxelsizes, medians):
    minvoxelsize = min(scalarvoxelsizes.values())

    meansspecs = []
    covsspecs = []

    for name in sorted(scans):
        resfrac = scalarvoxelsizes[name] / minvoxelsize
        mt = modalitytypes[name]

        if mt == 'T1':
            meansspecs.append('"exp(2.0;{0};{1}),exp(5.0;{0};{1}),exp(3.0;{0};{2})"'.format(medians['substantia_nigra'][name], 1.33 * medians['substantia_nigra'][name], 0.5 * medians['substantia_nigra'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.5 * resfrac * medians['substantia_nigra'][name]) ** 2))
        elif mt == 'T2':
            meansspecs.append('"exp(1.0;{0};{1}),exp(3.0;{0};{1}),exp(1.0;{0};{2})"'.format(medians['substantia_nigra'][name], 1.33 * medians['substantia_nigra'][name], 0.67 * medians['substantia_nigra'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.1 * resfrac * medians['substantia_nigra'][name]) ** 2))
        elif mt == 'WMS':
            meansspecs.append('"exp(1.0;{0};{1}),exp(4.0;{0};{1}),exp2(2.0;{0};{1};{2})"'.format(medians['substantia_nigra'][name], 0.5 * medians['substantia_nigra'][name], 1.33 * medians['substantia_nigra'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.5 * resfrac * medians['substantia_nigra'][name]) ** 2))
        elif mt == 'QSM':
            meansspecs.append('"exp(1.0;{0};{1}),exp(1.0;{0};{1}),exp(2.0;{0};{1})"'.format(medians['substantia_nigra'][name], medians['wm'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.3 * resfrac * medians['substantia_nigra'][name]) ** 2))
        else:
            raise Exception('Modality type ' + mt + ' is not supported for substantia nigra')

    return '--components=3 --priormeans={0} --priorcovcoefs={1}'.format(
            ','.join(meansspecs), ','.join(covsspecs))


def setup_subthalamic_nucleus(dirs, modalitytypes, scans, warpinv, scalarvoxelsizes, medians):
    minvoxelsize = min(scalarvoxelsizes.values())

    meansspecs = []
    covsspecs = []

    for name in sorted(scans):
        resfrac = scalarvoxelsizes[name] / minvoxelsize
        mt = modalitytypes[name]

        if mt == 'T2':
            meansspecs.append('"exp(1.0;{0};{1}),exp(3.0;{0};{1}),exp(1.0;{0};{2})"'.format(medians['subthalamic_nucleus'][name], 1.2 * medians['subthalamic_nucleus'][name], 0.5 * medians['subthalamic_nucleus'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.1 * resfrac * medians['subthalamic_nucleus'][name]) ** 2))
        elif mt == 'WMS':
            meansspecs.append('"exp(1.0;{0};{1}),exp(2.0;{0};{1}),exp(4.0;{0};{1})"'.format(medians['subthalamic_nucleus'][name], 0.67 * medians['subthalamic_nucleus'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.5 * resfrac * medians['subthalamic_nucleus'][name]) ** 2))
        elif mt == 'QSM':
            meansspecs.append('"step({0};{1}),exp(3.0;{0};{1}),exp(1.0;{0};{1})"'.format(medians['subthalamic_nucleus'][name], medians['wm'][name]))
            covsspecs.append('"flat({0}),flat({0}),flat({0})"'.format((0.3 * resfrac * medians['subthalamic_nucleus'][name]) ** 2))
        else:
            raise Exception('Modality type ' + mt + ' is not supported for subthalamic nucleus')

    return '--components=3 --priormeans={0} --priorcovcoefs={1}'.format(
            ','.join(meansspecs), ','.join(covsspecs))


def setup_common(dirs, modalitytypes, scans, csfmask, warp, scalarvoxelsizes, searchmm, deltastdev, mrfweight):
    names = sorted(scans)
    spacing = 0.5 * min(scalarvoxelsizes.values())
    # Make sure steps is even
    steps = 2 * int(searchmm / spacing)
    smooth = [0.5] * len(names)
    normlist = []
    for name in sorted(names):
        mt = modalitytypes[name]

        if mt == 'FA':
            normlist.append('none')
        elif mt == 'QSM':
            normlist.append('add')
        else:
            normlist.append('mul')
    
    return '--ftol=1e-10 --modalitynames={0} {1} --modalityimages={2} --normexclusion={3} --warp={4} --smoothness={5} --profilealpha=-3 --profilen0=3 --deltastdev={6} --profilespacing={7} --profilepoints={8} --steps={8} --usemrf --mrfweight={9} --mrfmeanfrac=0.0'.format(
            ','.join(names),
            '--normalisation=' + ','.join(normlist),
            ','.join([scans[name] for name in names]),
            csfmask,
            warp,
            ','.join([str(s) for s in smooth]),
            deltastdev / spacing,
            spacing,
            steps,
            mrfweight)


def generateall(dirs, modalitytypes, scans, csfmask, warp, warpinv, scalarvoxelsizes, outdir, structures = None):
    medians = get_all_medians(scans, warpinv, dirs)

    scans_putamen_pallidum = drop_modalities(scans, modalitytypes, ['GMPVE'])

    common_putamen_pallidum = mist_bin + ' --loglevel=info train ' + \
            setup_common(dirs, modalitytypes, scans_putamen_pallidum, csfmask, warp, scalarvoxelsizes, 3.0, 2.0, 10.0)

    cmds = []

    if structures is None or 'putamen' in structures:
        putamen_priors = setup_putamen(dirs, modalitytypes, scans_putamen_pallidum, warpinv, scalarvoxelsizes, medians)

        cmds.append(common_putamen_pallidum + ' --shape=' + mesh_dir + '/left_putamen.mim' + \
                ' --out=' + outdir + '/model_left_putamen.txt ' + putamen_priors + ' ' + ' '.join(dirs))

        cmds.append(common_putamen_pallidum + ' --shape=' + mesh_dir + '/right_putamen.mim' + \
                ' --out=' + outdir + '/model_right_putamen.txt ' + putamen_priors + ' ' + ' '.join(dirs))

    if structures is None or 'pallidum' in structures:
        pallidum_priors = setup_pallidum(dirs, modalitytypes, scans_putamen_pallidum, warpinv, scalarvoxelsizes, medians)
        
        cmds.append(common_putamen_pallidum + ' --shape=' + mesh_dir + '/left_pallidum.mim' + \
                ' --out=' + outdir + '/model_left_pallidum.txt ' + pallidum_priors + ' ' + ' '.join(dirs))

        cmds.append(common_putamen_pallidum + ' --shape=' + mesh_dir + '/right_pallidum.mim' + \
                ' --out=' + outdir + '/model_right_pallidum.txt ' + pallidum_priors + ' ' + ' '.join(dirs))

    
    if structures is None or 'caudate_accumbens' in structures:
        scans_caudate_accumbens = drop_modalities(scans, modalitytypes, ['FA', 'QSM', 'GMPVE'])

        common_caudate_accumbens = mist_bin + ' --loglevel=info train ' + \
                setup_common(dirs, modalitytypes, scans_caudate_accumbens, csfmask, warp, scalarvoxelsizes, 3.0, 2.0, 10.0)
        
        caudate_accumbens_priors = setup_caudate_accumbens(dirs, modalitytypes, scans_caudate_accumbens, warpinv, scalarvoxelsizes, medians)

        cmds.append(common_caudate_accumbens + ' --shape=' + mesh_dir + '/left_caudate_accumbens.mim' + \
                ' --out=' + outdir + '/model_left_caudate_accumbens.txt ' + caudate_accumbens_priors + ' ' + ' '.join(dirs))

        cmds.append(common_caudate_accumbens + ' --shape=' + mesh_dir + '/right_caudate_accumbens.mim' + \
                ' --out=' + outdir + '/model_right_caudate_accumbens.txt ' + caudate_accumbens_priors + ' ' + ' '.join(dirs))

    if structures is None or 'thalamus' in structures:
        scans_thalamus = drop_modalities(scans, modalitytypes, ['QSM', 'GMPVE'])
        
        common_thalamus = mist_bin + ' --loglevel=info train ' + \
                setup_common(dirs, modalitytypes, scans_thalamus, csfmask, warp, scalarvoxelsizes, 3.0, 2.0, 10.0)
        
        thalamus_priors = setup_thalamus(dirs, modalitytypes, scans_thalamus, warpinv, scalarvoxelsizes, medians)

        cmds.append(common_thalamus + ' --shape=' + mesh_dir + '/left_thalamus.mim' + \
                ' --out=' + outdir + '/model_left_thalamus.txt ' + thalamus_priors + ' ' + ' '.join(dirs))

        cmds.append(common_thalamus + ' --shape=' + mesh_dir + '/right_thalamus.mim' + \
                ' --out=' + outdir + '/model_right_thalamus.txt ' + thalamus_priors + ' ' + ' '.join(dirs))

    scans_hippocampus_amygdala = drop_modalities(scans, modalitytypes, ['FA', 'QSM'])

    if structures is None or 'hippocampus' in structures:
        common_hippocampus = mist_bin + ' --loglevel=info train ' + \
                setup_common(dirs, modalitytypes, scans_hippocampus_amygdala, csfmask, warp, scalarvoxelsizes, 3.0, 10.0, 10.0)

        hippocampus_priors = setup_hippocampus(dirs, modalitytypes, scans_hippocampus_amygdala, warpinv, scalarvoxelsizes, medians)

        cmds.append(common_hippocampus + ' --shape=' + mesh_dir + '/left_hippocampus.mim' + \
                ' --out=' + outdir + '/model_left_hippocampus.txt ' + hippocampus_priors + ' ' + ' '.join(dirs))

        cmds.append(common_hippocampus + ' --shape=' + mesh_dir + '/right_hippocampus.mim' + \
                ' --out=' + outdir + '/model_right_hippocampus.txt ' + hippocampus_priors + ' ' + ' '.join(dirs))

    if structures is None or 'amygdala' in structures:
        common_amygdala = mist_bin + ' --loglevel=info train ' + \
                setup_common(dirs, modalitytypes, scans_hippocampus_amygdala, csfmask, warp, scalarvoxelsizes, 3.0, 10.0, 5.0)

        amygdala_priors = setup_amygdala(dirs, modalitytypes, scans_hippocampus_amygdala, warpinv, scalarvoxelsizes, medians)
       
        cmds.append(common_amygdala + ' --shape=' + mesh_dir + '/left_amygdala.mim' + \
                ' --out=' + outdir + '/model_left_amygdala.txt ' + amygdala_priors + ' ' + ' '.join(dirs))
       
        cmds.append(common_amygdala + ' --shape=' + mesh_dir + '/right_amygdala.mim' + \
                ' --out=' + outdir + '/model_right_amygdala.txt ' + amygdala_priors + ' ' + ' '.join(dirs))
    
    scans_brainstem = drop_modalities(scans, modalitytypes, ['FA', 'T1', 'GMPVE'])
    
    common_rn_sn = mist_bin + ' --loglevel=info train ' + \
            setup_common(dirs, modalitytypes, scans_brainstem, csfmask, warp, scalarvoxelsizes, 2.0, 2.0, 1.0)
    
    if structures is None or 'red_nucleus' in structures:
        red_nucleus_priors = setup_red_nucleus(dirs, modalitytypes, scans_brainstem, warpinv, scalarvoxelsizes, medians)

        cmds.append(common_rn_sn + ' --mrfweight=10.0 --mrfmeanfrac=0.0 --shape=' + mesh_dir + \
                '/left_red_nucleus.mim --out=' + outdir + '/model_left_red_nucleus.txt ' + \
                red_nucleus_priors + ' ' + ' '.join(dirs))

        cmds.append(common_rn_sn + ' --mrfweight=10.0 --mrfmeanfrac=0.0 --shape=' + mesh_dir + \
                '/right_red_nucleus.mim --out=' + outdir + '/model_right_red_nucleus.txt ' + \
                red_nucleus_priors + ' ' + ' '.join(dirs))

    if structures is None or 'substantia_nigra' in structures:
        substantia_nigra_priors = setup_substantia_nigra(dirs, modalitytypes, scans_brainstem, warpinv, scalarvoxelsizes, medians)

        cmds.append(common_rn_sn + ' --mrfweight=10.0 --mrfmeanfrac=0.0 --shape=' + mesh_dir + \
                '/left_substantia_nigra.mim --out=' + outdir + '/model_left_substantia_nigra.txt ' + \
                substantia_nigra_priors + ' ' + ' '.join(dirs))

        cmds.append(common_rn_sn + ' --mrfweight=10.0 --mrfmeanfrac=0.0 --shape=' + mesh_dir + \
                '/right_substantia_nigra.mim --out=' + outdir + '/model_right_substantia_nigra.txt ' + \
                substantia_nigra_priors + ' ' + ' '.join(dirs))

    if structures is None or 'subthalamic_nucleus' in structures:
        common_stn = mist_bin + ' --loglevel=info train ' + \
                setup_common(dirs, modalitytypes, scans_brainstem, csfmask, warp, scalarvoxelsizes, 2.0, 2.0, 1.0)
        
        subthalamic_nucleus_priors = setup_subthalamic_nucleus(dirs, modalitytypes, scans_brainstem, warpinv, scalarvoxelsizes, medians)

        cmds.append(common_stn + ' --mrfweight=100.0 --mrfmeanfrac=0.0 --shape=' + mesh_dir + \
                '/left_subthalamic_nucleus.mim --out=' + outdir + '/model_left_subthalamic_nucleus.txt ' + \
                subthalamic_nucleus_priors + ' ' + ' '.join(dirs))

        cmds.append(common_stn + ' --mrfweight=100.0 --mrfmeanfrac=0.0 --shape=' + mesh_dir + \
                '/right_subthalamic_nucleus.mim --out=' + outdir + '/model_right_subthalamic_nucleus.txt ' + \
                subthalamic_nucleus_priors + ' ' + ' '.join(dirs))

    return cmds

