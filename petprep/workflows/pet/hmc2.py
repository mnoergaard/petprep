#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 16:39:38 2022

@author: martinnorgaard
"""

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl, freesurfer as fs

from ...config import DEFAULT_MEMORY_MIN_GB


def init_pet_hmc_wf(mem_gb, omp_nthreads, metadata, name='pet_hmc_wf'):
    """
    Build a workflow to estimate head-motion parameters.

    This workflow estimates the motion parameters to perform
    :abbr:`HMC (head motion correction)` over the input
    :abbr:`PET (Positron Emission Tomography)` image.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from petprep.workflows.pet import init_pet_hmc_wf
            wf = init_pet_hmc_wf(
                mem_gb=3,
                omp_nthreads=1)

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of PET file in GB
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    metadata : :obj:`dict`
        BIDS metadata for PET file
    name : :obj:`str`
        Name of workflow (default: ``pet_hmc_wf``)

    Inputs
    ------
    pet_file
        PET NIfTI file

    Outputs
    -------
    pet_mc_file
        ITKTransform file aligning each volume to ``ref_image``
    hmc_confounds
        MCFLIRT motion parameters, normalized to SPM format (X, Y, Z, Rx, Ry, Rz)
    translation
        Framewise displacement as measured by ``fsl_motion_outliers`` [Jenkinson2002]_.
    rotation
        Framewise displacement as measured by ``fsl_motion_outliers`` [Jenkinson2002]_.
    movement
        Framewise displacement as measured by ``fsl_motion_outliers`` [Jenkinson2002]_.

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    workflow = Workflow(name=name)
    workflow.__desc__ = """\
        Head-motion parameters with respect to the PET reference
        (transformation matrices, and six corresponding rotation and translation
         parameters) are estimated before any spatiotemporal filtering using
        `mcflirt` [FSL {fsl_ver}, @mcflirt].
        """.format(fsl_ver=fsl.Info().version() or '<ver>')

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['pet_file']),
        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['pet_mc_file', 'hmc_confounds', 'translation', 'rotation', 'movement']),
        name='outputnode')
    
    mid_frames = metadata['FrameTimesStart'] + metadata['FrameDuration']/2

    inputnode.min_frame = next(x for x, val in enumerate(mid_frames)
                                  if val > 120)   

    split_pet = pe.Node(interface = fs.MRIConvert(split = True), 
                     name = "split_pet")
    
    smooth_frame = pe.MapNode(interface=fsl.Smooth(fwhm=10), 
                           name="smooth_frame", 
                           iterfield=['in_file'])
    
    thres_frame = pe.MapNode(interface = fsl.maths.Threshold(thresh = 20, use_robust_range = True),
                          name = "thres_frame", 
                          iterfield = ['in_file'])
    
    estimate_motion = pe.Node(interface = fs.RobustTemplate(auto_detect_sensitivity = True,
                                            intensity_scaling = True,
                                            average_metric = 'mean',
                                            args = '--cras'),
                           name="estimate_motion", iterfield=['in_files'])
    
    correct_motion = pe.MapNode(interface = fs.ApplyVolTransform(), 
                             name = "correct_motion", 
                             iterfield = ['source_file', 'reg_file', 'transformed_file'])
    
    concat_frames = pe.Node(interface = fs.Concatenate(concatenated_file = 'mc.nii.gz'), 
                         name = "concat_frames")
    
    lta2xform = pe.MapNode(interface = fs.utils.LTAConvert(), 
                        name = "lta2xform", 
                        iterfield = ['in_lta', 'out_fsl'])
    
    est_trans_rot = pe.MapNode(interface = fsl.AvScale(all_param = True), 
                            name = "est_trans_rot", 
                            iterfield = ['mat_file'])
    
    upd_frame_list = pe.Node(pe.Function(input_names = ['in_file','min_frame'],
                                   output_names = ['upd_list_frames'],
                                   function = update_list_frames),
                          name = "upd_frame_list")
    
    upd_transform_list = pe.Node(pe.Function(input_names = ['in_file','min_frame'],
                                       output_names = ['upd_list_transforms'],
                                       function = update_list_transforms),
                              name = "upd_transform_list")
    
    hmc_movement_output = pe.Node(pe.Function(input_names = ['translations', 'rot_angles', 'rotation_translation_matrix','in_file'],
                                           output_names = ['hmc_confounds'],
                                           function = combine_hmc_outputs),
                               name = "hmc_movement_output")
    
    plot_motion = pe.Node(pe.Function(input_names = ['in_file'],
                                           function = plot_motion_outputs),
                               name = "plot_motion")
    
    # Connect workflow - init_pet_hmc_wf
    workflow.config['execution']['remove_unnecessary_outputs'] = 'false'
    workflow.connect([
        (inputnode, split_pet, [('pet', 'in_file')]),
        (split_pet,smooth_frame,[('out_file', 'in_file')]),
        (smooth_frame,thres_frame,[('smoothed_file', 'in_file')]),
        (thres_frame,upd_frame_list,[('out_file', 'in_file')]),
        (inputnode,upd_frame_list,[('min_frame', 'min_frame')]),
        (upd_frame_list,estimate_motion,[('upd_list_frames', 'in_files')]),
        (thres_frame,upd_transform_list,[('out_file', 'in_file')]),
        (inputnode,upd_transform_list,[('min_frame', 'min_frame')]),
        (upd_transform_list,estimate_motion,[('upd_list_transforms', 'transform_outputs')]),
        (split_pet,correct_motion,[('out_file', 'source_file')]),
        (estimate_motion,correct_motion,[('transform_outputs', 'reg_file')]),
        (estimate_motion,correct_motion,[('out_file', 'target_file')]),
        (split_pet,correct_motion,[(('out_file', add_mc_ext), 'transformed_file')]),
        (correct_motion,concat_frames,[('transformed_file', 'in_files')]),
        (estimate_motion,lta2xform,[('transform_outputs', 'in_lta')]),
        (estimate_motion,lta2xform,[(('transform_outputs', lta2mat), 'out_fsl')]),
        (lta2xform,est_trans_rot,[('out_fsl', 'mat_file')]),
        (est_trans_rot,hmc_movement_output,[('translations', 'translations'),('rot_angles', 'rot_angles'),('rotation_translation_matrix','rotation_translation_matrix')]),
        (upd_frame_list,hmc_movement_output,[('upd_list_frames', 'in_file')]),
        (hmc_movement_output,plot_motion,[('hmc_confounds','in_file')]),
        (concat_frames, outputnode, [('transformed_file', 'pet_mc_file')]),
        (hmc_movement_output, outputnode, [('hmc_confounds', 'hmc_confounds')]),
        (plot_motion, outputnode, [('translation', 'translation')]),
        (plot_motion, outputnode, [('rotation', 'rotation')]),
        (plot_motion, outputnode, [('movement', 'movement')])
                         ])
    return workflow

# HELPER FUNCTIONS
def update_list_frames(in_file, min_frame):
    
    """
    
    Arguments
    ---------
    """
    
    new_list = [in_file[min_frame]] * min_frame + in_file[min_frame:]
    return new_list

def update_list_transforms(in_file, min_frame):
    
    """
    
    Arguments
    ---------
    """
    
    new_list = [in_file[min_frame]] * min_frame + in_file[min_frame:]
    lta_list = [ext.replace('nii.gz','lta') for ext in new_list]  
    return lta_list

def add_mc_ext(in_file):
    
    """
    
    Arguments
    ---------
    """
    
    if len(in_file) > 1:
        mc_list = [ext.replace('.nii.gz','_mc.nii.gz') for ext in in_file] # list comphrehension
    else:
        mc_list = in_file.replace('.nii.gz','_mc.nii.gz')
    return mc_list

def lta2mat(in_file):
    
    """
    
    Arguments
    ---------
    """
    
    mat_list = [ext.replace('.lta','.mat') for ext in in_file]
    return mat_list 

def combine_hmc_outputs(translations, rot_angles, rotation_translation_matrix, in_file):
    
    """
    
    Arguments
    ---------
    """
    
    import os
    import pandas as pd
    import numpy as np
    import nibabel as nib
    
    new_pth = os.getcwd()
    
    movement = []
    for idx, trans in enumerate(translations):
        
        img = nib.load(in_file[idx])
        vox_ind = np.asarray(np.nonzero(img.get_fdata()))
        pos_bef = np.concatenate((vox_ind,np.ones((1,len(vox_ind[0,:])))))
        pos_aft = rotation_translation_matrix[idx] @ pos_bef
        diff_pos = pos_bef-pos_aft
        
        max_x = abs(max(diff_pos[0,:], key=abs))
        max_y = abs(max(diff_pos[1,:], key=abs))
        max_z = abs(max(diff_pos[2,:], key=abs))
        overall = np.sqrt(diff_pos[0,:] ** 2 + diff_pos[1,:] ** 2 + diff_pos[2,:] ** 2)
        max_tot = np.max(overall)
        median_tot = np.median(overall)
        
        movement.append(np.concatenate((translations[idx],rot_angles[idx], [max_x, max_y, max_z, max_tot, median_tot])))
        
    confounds = pd.DataFrame(movement, columns=['trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'max_x', 'max_y', 
                                                'max_z', 'max_tot', 'median_tot'])
    confounds.to_csv(os.path.join(new_pth,'hmc_confounds.tsv'), sep='\t')
    # np.savetxt(os.path.join(new_pth,'hmc_confounds.tsv'), movement, fmt='%10.5f', delimiter='\t', header='trans_x    trans_y    trans_z    rot_x    rot_y    rot_z')
    
    return os.path.join(new_pth,'hmc_confounds.tsv')

def plot_motion_outputs(in_file):
    
    """
    
    Arguments
    ---------
    """
    
    import os
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    
    confounds = pd.read_csv(in_file, sep='\t')
    
    new_pth = os.getcwd()
    
    n_frames = len(confounds.index)
    
    plt.figure(figsize=(11,5))
    plt.plot(np.arange(0,n_frames), confounds['trans_x'], "-r", label='trans_x')
    plt.plot(np.arange(0,n_frames), confounds['trans_y'], "-g", label='trans_y')
    plt.plot(np.arange(0,n_frames), confounds['trans_z'], "-b", label='trans_z')
    plt.legend(loc="upper left")
    plt.ylabel('Translation [mm]')
    plt.xlabel('frame #')
    plt.grid(visible=True)
    plt.savefig(os.path.join(new_pth,'translation.png'), format='png')
    plt.close()
    
    plt.figure(figsize=(11,5))
    plt.plot(np.arange(0,n_frames), confounds['rot_x'], "-r", label='rot_x')
    plt.plot(np.arange(0,n_frames), confounds['rot_y'], "-g", label='rot_y')
    plt.plot(np.arange(0,n_frames), confounds['rot_z'], "-b", label='rot_z')
    plt.legend(loc="upper left")
    plt.ylabel('Rotation [degrees]')
    plt.xlabel('frame #')
    plt.grid(visible=True)
    plt.savefig(os.path.join(new_pth,'rotation.png'), format='png')
    plt.close()
    
    plt.figure(figsize=(11,5))
    plt.plot(np.arange(0,n_frames), confounds['max_x'], "--r", label='max_x')
    plt.plot(np.arange(0,n_frames), confounds['max_y'], "--g", label='max_y')
    plt.plot(np.arange(0,n_frames), confounds['max_z'], "--b", label='max_z')
    plt.plot(np.arange(0,n_frames), confounds['max_tot'], "-k", label='max_total')
    plt.plot(np.arange(0,n_frames), confounds['median_tot'], "-m", label='median_tot')
    plt.legend(loc="upper left")
    plt.ylabel('Movement [mm]')
    plt.xlabel('frame #')
    plt.grid(visible=True)
    plt.savefig(os.path.join(new_pth,'movement.png'), format='png')
    plt.close()
    
    translation = os.path.join(new_pth,'translation.png')
    rotation = os.path.join(new_pth,'rotation.png')
    movement = os.path.join(new_pth,'movement.png')
    
    return translation, rotation, movement
