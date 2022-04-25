#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 11:45:37 2022

Set of nipype classes to run FreeSurfer subsegmentation of the brainstem, thalamus, hippocampus and amygdala.

@author: martinnorgaard
"""

from nipype.interfaces.base import (
                                    CommandLine, 
                                    CommandLineInputSpec, 
                                    File,
                                    traits,
                                    Directory,
                                    TraitedSpec
                                    )
import os

class SegmentHippocampusAmygdalaInputSpec(CommandLineInputSpec):
    subject_id = traits.Str(mandatory=True, argstr="sub-%s", desc = "Subject ID", position=1)
    fs_subjects_dir = traits.Str(mandatory=True, argstr="%s", desc="FreeeSurfer subjects directory", position=2)
    
class SegmentHippocampusAmygdalaOutputSpec(TraitedSpec):
    left_hippoamyg_seg_file = File(desc = "Output volume", exists = True)
    right_hippoamyg_seg_file = File(desc = "Output volume", exists = True)
    left_hippoamyg_vol_file = File(desc = "Output volume", exists = True)
    right_hippoamyg_vol_file = File(desc = "Output volume", exists = True)

class SegmentHippocampusAmygdala(CommandLine):
    _cmd = 'segmentHA_T1.sh'
    input_spec = SegmentHippocampusAmygdalaInputSpec
    output_spec = SegmentHippocampusAmygdalaOutputSpec
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['left_hippoamyg_seg_file'] = os.path.join(self.inputs.fs_subjects_dir, f"sub-{self.inputs.subject_id}", 'mri','lh.hippoAmygLabels-T1.v21.FSvoxelSpace.mgz')
        outputs['right_hippoamyg_seg_file'] = os.path.join(self.inputs.fs_subjects_dir, f"sub-{self.inputs.subject_id}", 'mri','rh.hippoAmygLabels-T1.v21.FSvoxelSpace.mgz')
        outputs['left_hippoamyg_vol_file'] = os.path.join(self.inputs.fs_subjects_dir, f"sub-{self.inputs.subject_id}", 'mri','lh.amygNucVolumes-T1.v21')
        outputs['right_hippoamyg_vol_file'] = os.path.join(self.inputs.fs_subjects_dir, f"sub-{self.inputs.subject_id}", 'mri','rh.amygNucVolumes-T1.v21')
        return outputs


class SegmentThalamusInputSpec(CommandLineInputSpec):
    subject_id = traits.Str(mandatory=True, argstr="sub-%s", desc = "Subject ID", position=1)
    fs_subjects_dir = traits.Str(mandatory=True, argstr="%s", desc="FreeeSurfer subjects directory", position=2)
    
class SegmentThalamusOutputSpec(TraitedSpec):
    thalamus_seg_file = File(desc = "Output volume", exists = True)
    thalamus_vol_file = File(desc = "Output volume", exists = True)

class SegmentThalamus(CommandLine):
    _cmd = 'segmentThalamicNuclei.sh'
    input_spec = SegmentThalamusInputSpec
    output_spec = SegmentThalamusOutputSpec
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['thalamus_seg_file'] = os.path.join(self.inputs.fs_subjects_dir, f"sub-{self.inputs.subject_id}", 'mri','ThalamicNuclei.v12.T1.FSvoxelSpace.mgz')
        outputs['thalamus_vol_file'] = os.path.join(self.inputs.fs_subjects_dir, f"sub-{self.inputs.subject_id}", 'mri','ThalamicNuclei.v12.T1.volumes')
        return outputs
    
class SegmentBrainstemInputSpec(CommandLineInputSpec):
    subject_id = traits.Str(mandatory=True, argstr="sub-%s", desc = "Subject ID", position=1)
    fs_subjects_dir = traits.Str(mandatory=True, argstr="%s", desc="FreeeSurfer subjects directory", position=2)
    
class SegmentBrainstemOutputSpec(TraitedSpec):
    brainstem_seg_file = File(desc = "Output volume", exists = True)

class SegmentBrainstem(CommandLine):
    _cmd = 'segmentBS.sh'
    input_spec = SegmentBrainstemInputSpec
    output_spec = SegmentBrainstemOutputSpec
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['brainstem_seg_file'] = os.path.join(self.inputs.fs_subjects_dir, f"sub-{self.inputs.subject_id}", 'mri','brainstemSsLabels.v12.FSvoxelSpace.mgz')
        outputs['brainstem_vol_file'] = os.path.join(self.inputs.fs_subjects_dir, f"sub-{self.inputs.subject_id}", 'mri','brainstemSsVolumes.v12.txt')
        return outputs
    
    