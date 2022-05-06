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
    
    