#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 17:46:59 2022

@author: martinnorgaard
"""

"""
Head-Motion Estimation and Correction (HMC) of PET images
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_pet_hmc_wf

"""

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl, freesurfer as fs

def init_pet_pvc_wf(mem_gb, omp_nthreads, name='pet_pvc_wf'):
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
    name : :obj:`str`
        Name of workflow (default: ``pet_hmc_wf``)

    Inputs
    ------
    pet_file
        PET NIfTI file

    Outputs
    -------
    xforms
        ITKTransform file aligning each volume to ``ref_image``
    movpar_file
        MCFLIRT motion parameters, normalized to SPM format (X, Y, Z, Rx, Ry, Rz)
    rms_file
        Framewise displacement as measured by ``fsl_motion_outliers`` [Jenkinson2002]_.

    """
    
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['pet_file', 
                                      '']),
        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['xforms', 'movpar_file', 'rmsd_file']),
        name='outputnode')
    
    