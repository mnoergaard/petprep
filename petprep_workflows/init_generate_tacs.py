#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 14:30:45 2022



@author: martinnorgaard
"""

import nipype.interfaces.freesurfer as fs
import os
import pandas as pd
import nibabel as nib

def init_seg2tacs(pet_t1, seg_file, metadata, name='seg2tacs_wf'):

    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    workflow = Workflow(name=name)
    
    ss = fs.SegStats()
    ss.inputs.segmentation_file = seg_file
    ss.inputs.in_file = pet_t1
    ss.inputs.color_table_file = os.path.join(os.environ['FREESURFER_HOME'],'FreeSurferColorLUT.txt')
    ss.inputs.avgwf_file = True
    ss.inputs.avgwf_txt_file = True
    ss.inputs.exclude_id = 0
    out = ss.run()
    
    return workflow


def seg2table(summary_file, avgwf_txt_file):
    
    """
    
    A function to compute a time weighted average over
        the time frames for a given pet volume

        Parameters
        ----------
        summary_file : str 
            summary file path (str) from mri_segstats of PET data and segmentation
        avgwf_txt_file : str 
            avgwf_txt_file file path (str) of TACS from PET

        Returns
        -------
        out_file : str 
            output file path (str) of TACS file

    """
    
    new_pth = os.getcwd()
    
    summary = pd.read_csv(summary_file, comment = '#', delim_whitespace = True,
                         names=['Index', 
                                'SegId', 
                                'NVoxels', 
                                'Volume_mm3', 
                                'StructName', 
                                'Mean', 
                                'StdDev', 
                                'Min', 
                                'Max', 
                                'Range']).dropna(axis = 0)
    
    tacs = pd.read_csv(avgwf_txt_file, 
                       delim_whitespace = True, 
                       names = summary['StructName'].tolist())
    
    tacs.to_csv(os.path.join(new_pth, avgwf_txt_file.replace('_avgwf.txt', '_tacs.tsv')), sep='\t')
    
    return os.path.join(new_pth, avgwf_txt_file.replace('_avgwf.txt', '_tacs.tsv'))
