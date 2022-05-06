#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  5 17:56:52 2022

@author: martinnorgaard
"""

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, freesurfer as fs

from ...config import DEFAULT_MEMORY_MIN_GB

def init_gtmseg_wf(mem_gb, omp_nthreads, name='gtmseg_wf'):
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
        Name of workflow (default: ``gtmseg_wf``)

    Inputs
    ------
    subject_id
        ID of subject to analyze
    fs_subjects_dir
        FreeSurfer subject directory
        

    Outputs
    -------
    gtmseg_file
        file of anatomical segmentation for the geometric transfer matrix (GTM)

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    
    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['subject_id', 'fs_subjects_dir']),
        name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['gtmseg_file']),
        name='outputnode')
    
    gtmseg = pe.Node(fs.petsurfer.GTMSeg(),
                     name="gtmseg",
                     mem_gb=DEFAULT_MEMORY_MIN_GB)
    
    workflow.connect([
        (inputnode, gtmseg, [('subject_id', 'subject_id'),
                             ('fs_subjects_dir', 'subjects_dir')]),
        (gtmseg, outputnode, [('out_file', 'gtmseg_file')])
        ])

    return workflow
    