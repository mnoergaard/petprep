# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright 2021 The NiPreps Developers <nipreps@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# We support and encourage derived works from this project, please read
# about our expectations at
#
#     https://www.nipreps.org/community/licensing/
#
"""
Orchestrating the PET-preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_pet_preproc_wf
.. autofunction:: init_pet_derivatives_wf

"""
from ... import config

import os

import nibabel as nb
from nipype.interfaces.fsl import Split as FSLSplit
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from niworkflows.utils.connections import pop_file, listify

from ...interfaces import DerivativesDataSink
from ...interfaces.reports import PETSummary

# PET workflows
from .hmc import init_pet_hmc_wf
from .registration import init_pet_t1_trans_wf, init_pet_reg_wf
from .resampling import (
    init_pet_surf_wf,
    init_pet_std_trans_wf,
    init_pet_preproc_trans_wf,
)
from .outputs import init_pet_derivatives_wf


def init_pet_preproc_wf(pet_file):
    """
    This workflow controls the PET preprocessing stages of *PETPrep*.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from petprep.workflows.tests import mock_config
            from petprep import config
            from petprep.workflows.pet.base import init_pet_preproc_wf
            with mock_config():
                pet_file = config.execution.bids_dir / "sub-01" / "pet" \
                    / "sub-01_ses-baseline_pet.nii.gz"
                wf = init_pet_preproc_wf(str(pet_file))

    Parameters
    ----------
    pet_file
        Path to NIfTI file

    Inputs
    ------
    pet_file
        PET NIfTI file
    t1w_preproc
        Bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    t1w_asec
        Segmentation of structural image, done with FreeSurfer.
    t1w_aparc
        Parcellation of structural image, done with FreeSurfer.
    t1w_tpms
        List of tissue probability maps in T1w space
    template
        List of templates to target
    anat2std_xfm
        List of transform files, collated with templates
    std2anat_xfm
        List of inverse transform files, collated with templates
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    t1w2fsnative_xfm
        LTA-style affine matrix translating from T1w to FreeSurfer-conformed subject space
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed subject space to T1w

    Outputs
    -------
    pet_t1
        PET data resampled to T1w space
    pet_mask_t1
        PET data mask in T1w space
    pet_std
        PET data resampled to template space
    pet_mask_std
        PET data mask in template space
    confounds
        TSV of confounds
    surfaces
        PET data resampled to FreeSurfer surfaces
    pet_cifti
        PET CIFTI image
    cifti_variant
        combination of target spaces for `pet_cifti`

    See Also
    --------

    * :py:func:`~niworkflows.func.util.init_pet_reference_wf`
    * :py:func:`~petprep.workflows.pet.hmc.init_pet_hmc_wf`
    * :py:func:`~petprep.workflows.pet.registration.init_pet_t1_trans_wf`
    * :py:func:`~petprep.workflows.pet.registration.init_pet_reg_wf`
    * :py:func:`~petprep.workflows.pet.confounds.init_pet_confs_wf`
    * :py:func:`~petprep.workflows.pet.resampling.init_pet_std_trans_wf`
    * :py:func:`~petprep.workflows.pet.resampling.init_pet_preproc_trans_wf`
    * :py:func:`~petprep.workflows.pet.resampling.init_pet_surf_wf`

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nibabel import ApplyMask
    from niworkflows.interfaces.utility import KeySelect, DictMerge

    if nb.load(
        pet_file[0] if isinstance(pet_file, (list, tuple)) else pet_file
    ).shape[3:] <= (5 - config.execution.sloppy,):
        config.loggers.workflow.warning(
            f"Too short PET series (<= 5 timepoints). Skipping processing of <{pet_file}>."
        )
        return

    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}
    pet_tlen = 10

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    freesurfer = config.workflow.run_reconall
    spaces = config.workflow.spaces
    petprep_dir = str(config.execution.petprep_dir)

    # Extract BIDS entities and metadata from PET file(s)
    entities = extract_entities(pet_file)
    layout = config.execution.layout

    # Extract metadata
    all_metadata = [layout.get_metadata(fname) for fname in listify(pet_file)]

    # Take first file as reference
    ref_file = pop_file(pet_file)
    metadata = all_metadata[0]
    # get original image orientation
    ref_orientation = get_img_orientation(ref_file)

    if os.path.isfile(ref_file):
        bold_tlen, mem_gb = _create_mem_gb(ref_file)

    wf_name = _get_wf_name(ref_file)
    config.loggers.workflow.debug(
        "Creating bold processing workflow for <%s> (%.2f GB / %d TRs). "
        "Memory resampled/largemem=%.2f/%.2f GB.",
        ref_file,
        mem_gb["filesize"],
        bold_tlen,
        mem_gb["resampled"],
        mem_gb["largemem"],
    )

    # Build workflow
    workflow = Workflow(name=wf_name)
    workflow.__postdesc__ = """\
All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. head-motion
transform matrices, susceptibility distortion correction when available,
and co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `antsApplyTransforms` (ANTs),
configured with Lanczos interpolation to minimize the smoothing
effects of other kernels [@lanczos].
Non-gridded (surface) resamplings were performed using `mri_vol2surf`
(FreeSurfer).
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "pet_file",
                "subjects_dir",
                "subject_id",
                "t1w_preproc",
                "t1w_mask",
                "t1w_dseg",
                "t1w_tpms",
                "t1w_aseg",
                "t1w_aparc",
                "anat2std_xfm",
                "std2anat_xfm",
                "template",
                "t1w2fsnative_xfm",
                "fsnative2t1w_xfm",
            ]
        ),
        name="inputnode",
    )
    inputnode.inputs.pet_file = pet_file

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "pet_t1",
                "pet_t1_ref",
                "pet2anat_xfm",
                "anat2pet_xfm",
                "pet_mask_t1",
                "pet_aseg_t1",
                "pet_aparc_t1",
                "pet_std",
                "pet_std_ref",
                "pet_mask_std",
                "pet_aseg_std",
                "pet_aparc_std",
                "pet_native",
                "pet_native_ref",
                "pet_mask_native",
                "pet_echos_native",
                "pet_cifti",
                "cifti_variant",
                "cifti_metadata",
                "cifti_density",
                "surfaces",
                "confounds",
                "nonaggr_denoised_file",
                "confounds_metadata",
            ]
        ),
        name="outputnode",
    )

    # Generate a brain-masked conversion of the t1w
    t1w_brain = pe.Node(ApplyMask(), name="t1w_brain")

    # PET source: track original PET file(s)
    pet_source = pe.Node(niu.Select(inlist=pet_file), name="pet_source")


    summary = pe.Node(
        PETSummary(
            registration=("FSL", "FreeSurfer")[freesurfer],
            registration_dof=config.workflow.pet2t1w_dof,
            registration_init=config.workflow.pet2t1w_init,
            orientation=ref_orientation,
        ),
        name="summary",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )
    summary.inputs.dummy_scans = config.workflow.dummy_scans

    pet_derivatives_wf = init_pet_derivatives_wf(
        bids_root=layout.root,
        cifti_output=config.workflow.cifti_output,
        freesurfer=freesurfer,
        all_metadata=all_metadata,
        output_dir=petprep_dir,
        spaces=spaces,
    )
    pet_derivatives_wf.inputs.inputnode.all_source_files = pet_file

    # fmt:off
    workflow.connect([
        (outputnode, pet_derivatives_wf, [
            ("pet_t1", "inputnode.pet_t1"),
            ("pet_t1_ref", "inputnode.pet_t1_ref"),
            ("pet2anat_xfm", "inputnode.pet2anat_xfm"),
            ("pet2bold_xfm", "inputnode.pet2bold_xfm"),
            ("pet_aseg_t1", "inputnode.pet_aseg_t1"),
            ("pet_aparc_t1", "inputnode.pet_aparc_t1"),
            ("pet_mask_t1", "inputnode.pet_mask_t1"),
            ("pet_native", "inputnode.pet_native"),
            ("pet_native_ref", "inputnode.pet_native_ref"),
            ("pet_mask_native", "inputnode.pet_mask_native"),
            ("surfaces", "inputnode.surf_files"),
            ("pet_cifti", "inputnode.pet_cifti"),
            ("cifti_variant", "inputnode.cifti_variant"),
            ("cifti_metadata", "inputnode.cifti_metadata"),
            ("cifti_density", "inputnode.cifti_density")
        ]),
    ])
    # fmt:on

    # Generate a tentative boldref
    initial_petref_wf = init_pet_reference_wf(
        name="initial_petref_wf",
        omp_nthreads=omp_nthreads,
        pet_file=pet_file
    )
    initial_petref_wf.inputs.inputnode.dummy_scans = config.workflow.dummy_scans

    # Select validated PET files (orientations checked or corrected)
    select_pet = pe.Node(niu.Select(), name="select_pet")

    # HMC on the PET
    pet_hmc_wf = init_pet_hmc_wf(
        name="pet_hmc_wf", mem_gb=mem_gb["filesize"], omp_nthreads=omp_nthreads
    )

    # calculate BOLD registration to T1w
    pet_reg_wf = init_pet_reg_wf(
        pet2t1w_dof=config.workflow.pet2t1w_dof,
        pet2t1w_init=config.workflow.pet2t1w_init,
        freesurfer=freesurfer,
        mem_gb=mem_gb["resampled"],
        name="pet_reg_wf",
        omp_nthreads=omp_nthreads,
        sloppy=config.execution.sloppy,
        use_coreg=config.workflow.use_coreg,
        use_compression=False,
    )

    # apply PET registration to T1w
    pet_t1_trans_wf = init_pet_t1_trans_wf(
        name="pet_t1_trans_wf",
        freesurfer=freesurfer,
        mem_gb=mem_gb["resampled"],
        omp_nthreads=omp_nthreads,
        use_compression=False,
    )
    pet_t1_trans_wf.inputs.inputnode.fieldwarp = "identity"

    pet_final = pe.Node(
        niu.IdentityInterface(fields=["bold", "boldref", "mask", "bold_echos"]),
        name="bold_final"
    )


    # MAIN WORKFLOW STRUCTURE #######################################################
    # fmt:off
    workflow.connect([
        # Prepare masked T1w image
        (inputnode, t1w_brain, [("t1w_preproc", "in_file"),
                                ("t1w_mask", "in_mask")]),
        # Select validated bold files per-echo
        (initial_petref_wf, select_pet, [("outputnode.all_pet_files", "inlist")]),
        # HMC
        (initial_petref_wf, pet_hmc_wf, [
            ("outputnode.raw_ref_image", "inputnode.raw_ref_image"),
            ("outputnode.pet_file", "inputnode.pet_file"),
        ]),
        # PET-T1w registration workflow
        (inputnode, pet_reg_wf, [
            ("t1w_dseg", "inputnode.t1w_dseg"),
            # Undefined if --fs-no-reconall, but this is safe
            ("subjects_dir", "inputnode.subjects_dir"),
            ("subject_id", "inputnode.subject_id"),
            ("fsnative2t1w_xfm", "inputnode.fsnative2t1w_xfm"),
        ]),
        (pet_final, pet_reg_wf, [
            ("petref", "inputnode.ref_pet_brain")]),
        (t1w_brain, pet_reg_wf, [("out_file", "inputnode.t1w_brain")]),
        (inputnode, pet_t1_trans_wf, [
            ("pet_file", "inputnode.name_source"),
            ("t1w_mask", "inputnode.t1w_mask"),
            ("t1w_aseg", "inputnode.t1w_aseg"),
            ("t1w_aparc", "inputnode.t1w_aparc"),
        ]),
        (t1w_brain, pet_t1_trans_wf, [("out_file", "inputnode.t1w_brain")]),
        (pet_reg_wf, outputnode, [
            ("outputnode.itk_pet_to_t1", "pet2anat_xfm"),
            ("outputnode.itk_t1_to_pet", "anat2pet_xfm"),
        ]),
        (pet_reg_wf, pet_t1_trans_wf, [
            ("outputnode.itk_pet_to_t1", "inputnode.itk_pet_to_t1"),
        ]),
        (pet_final, pet_t1_trans_wf, [
            ("mask", "inputnode.ref_pet_mask"),
            ("petref", "inputnode.ref_pet_brain"),
        ]),
        (pet_t1_trans_wf, outputnode, [
            ("outputnode.pet_t1", "pet_t1"),
            ("outputnode.pet_t1_ref", "pet_t1_ref"),
            ("outputnode.pet_aseg_t1", "pet_aseg_t1"),
            ("outputnode.pet_aparc_t1", "pet_aparc_t1"),
        ]),
        # Connect bold_confounds_wf
        (inputnode, bold_confounds_wf, [
            ("t1w_tpms", "inputnode.t1w_tpms"),
            ("t1w_mask", "inputnode.t1w_mask"),
        ]),
        (bold_hmc_wf, bold_confounds_wf, [
            ("outputnode.movpar_file", "inputnode.movpar_file"),
            ("outputnode.rmsd_file", "inputnode.rmsd_file"),
        ]),
        (bold_reg_wf, bold_confounds_wf, [
            ("outputnode.itk_t1_to_bold", "inputnode.t1_bold_xform")
        ]),
        (initial_boldref_wf, bold_confounds_wf, [
            ("outputnode.skip_vols", "inputnode.skip_vols"),
        ]),
        (initial_boldref_wf, final_boldref_wf, [
            ("outputnode.skip_vols", "inputnode.dummy_scans"),
        ]),
        (final_boldref_wf, bold_final, [
            ("outputnode.ref_image", "boldref"),
            ("outputnode.bold_mask", "mask"),
        ]),
        (bold_final, bold_confounds_wf, [
            ("bold", "inputnode.bold"),
            ("mask", "inputnode.bold_mask"),
        ]),
        (bold_confounds_wf, outputnode, [
            ("outputnode.confounds_file", "confounds"),
            ("outputnode.confounds_metadata", "confounds_metadata"),
            ("outputnode.acompcor_masks", "acompcor_masks"),
            ("outputnode.tcompcor_mask", "tcompcor_mask"),
        ]),
        # Native-space BOLD files (if calculated)
        (bold_final, outputnode, [
            ("bold", "bold_native"),
            ("boldref", "bold_native_ref"),
            ("mask", "bold_mask_native"),
            ("bold_echos", "bold_echos_native"),
        ]),
        # Summary
        (initial_boldref_wf, summary, [("outputnode.algo_dummy_scans", "algo_dummy_scans")]),
        (bold_reg_wf, summary, [("outputnode.fallback", "fallback")]),
        (outputnode, summary, [("confounds", "confounds_file")]),
        # Select echo indices for original/validated BOLD files
        (echo_index, bold_source, [("echoidx", "index")]),
        (echo_index, select_bold, [("echoidx", "index")]),
    ])
    # fmt:on

    # for standard EPI data, pass along correct file
    if not multiecho:
        # fmt:off
        workflow.connect([
            (inputnode, func_derivatives_wf, [("bold_file", "inputnode.source_file")]),
            (bold_split, bold_t1_trans_wf, [("out_files", "inputnode.bold_split")]),
            (bold_hmc_wf, bold_t1_trans_wf, [("outputnode.xforms", "inputnode.hmc_xforms")]),
        ])
        # fmt:on
    else:  # for meepi, use optimal combination
        # fmt:off
        workflow.connect([
            # update name source for optimal combination
            (inputnode, func_derivatives_wf, [
                (("bold_file", combine_meepi_source), "inputnode.source_file"),
            ]),
            (join_echos, bold_t2s_wf, [("bold_files", "inputnode.bold_file")]),
            (join_echos, bold_final, [("bold_files", "bold_echos")]),
            (bold_t2s_wf, split_opt_comb, [("outputnode.bold", "in_file")]),
            (split_opt_comb, bold_t1_trans_wf, [("out_files", "inputnode.bold_split")]),
            (bold_t2s_wf, bold_final, [("outputnode.bold", "bold")]),
        ])
        # fmt:on

        # Already applied in bold_bold_trans_wf, which inputs to bold_t2s_wf
        bold_t1_trans_wf.inputs.inputnode.hmc_xforms = "identity"

    # Map final BOLD mask into T1w space (if required)
    nonstd_spaces = set(spaces.get_nonstandard())
    if nonstd_spaces.intersection(("T1w", "anat")):
        from niworkflows.interfaces.fixes import (
            FixHeaderApplyTransforms as ApplyTransforms,
        )

        boldmask_to_t1w = pe.Node(
            ApplyTransforms(interpolation="MultiLabel"),
            name="boldmask_to_t1w",
            mem_gb=0.1,
        )
        # fmt:off
        workflow.connect([
            (bold_reg_wf, boldmask_to_t1w, [("outputnode.itk_bold_to_t1", "transforms")]),
            (bold_t1_trans_wf, boldmask_to_t1w, [("outputnode.bold_mask_t1", "reference_image")]),
            (bold_final, boldmask_to_t1w, [("mask", "input_image")]),
            (boldmask_to_t1w, outputnode, [("output_image", "bold_mask_t1")]),
        ])
        # fmt:on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        # Only use uncompressed output if AROMA is to be run
        bold_std_trans_wf = init_bold_std_trans_wf(
            freesurfer=freesurfer,
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            name="bold_std_trans_wf",
            use_compression=not config.execution.low_mem,
        )
        bold_std_trans_wf.inputs.inputnode.fieldwarp = "identity"

        # fmt:off
        workflow.connect([
            (inputnode, bold_std_trans_wf, [
                ("template", "inputnode.templates"),
                ("anat2std_xfm", "inputnode.anat2std_xfm"),
                ("bold_file", "inputnode.name_source"),
                ("t1w_aseg", "inputnode.bold_aseg"),
                ("t1w_aparc", "inputnode.bold_aparc"),
            ]),
            (bold_final, bold_std_trans_wf, [
                ("mask", "inputnode.bold_mask"),
            ]),
            (bold_reg_wf, bold_std_trans_wf, [
                ("outputnode.itk_bold_to_t1", "inputnode.itk_bold_to_t1"),
            ]),
            (bold_std_trans_wf, outputnode, [
                ("outputnode.bold_std", "bold_std"),
                ("outputnode.bold_std_ref", "bold_std_ref"),
                ("outputnode.bold_mask_std", "bold_mask_std"),
            ]),
        ])
        # fmt:on

        if freesurfer:
            # fmt:off
            workflow.connect([
                (bold_std_trans_wf, func_derivatives_wf, [
                    ("outputnode.bold_aseg_std", "inputnode.bold_aseg_std"),
                    ("outputnode.bold_aparc_std", "inputnode.bold_aparc_std"),
                ]),
                (bold_std_trans_wf, outputnode, [
                    ("outputnode.bold_aseg_std", "bold_aseg_std"),
                    ("outputnode.bold_aparc_std", "bold_aparc_std"),
                ]),
            ])
            # fmt:on

        if not multiecho:
            # fmt:off
            workflow.connect([
                (bold_split, bold_std_trans_wf, [("out_files", "inputnode.bold_split")]),
                (bold_hmc_wf, bold_std_trans_wf, [
                    ("outputnode.xforms", "inputnode.hmc_xforms"),
                ]),
            ])
            # fmt:on
        else:
            # fmt:off
            workflow.connect([
                (split_opt_comb, bold_std_trans_wf, [("out_files", "inputnode.bold_split")])
            ])
            # fmt:on

            # Already applied in bold_bold_trans_wf, which inputs to bold_t2s_wf
            bold_std_trans_wf.inputs.inputnode.hmc_xforms = "identity"

        # fmt:off
        # func_derivatives_wf internally parametrizes over snapshotted spaces.
        workflow.connect([
            (bold_std_trans_wf, func_derivatives_wf, [
                ("outputnode.template", "inputnode.template"),
                ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ("outputnode.bold_std_ref", "inputnode.bold_std_ref"),
                ("outputnode.bold_std", "inputnode.bold_std"),
                ("outputnode.bold_mask_std", "inputnode.bold_mask_std"),
            ]),
        ])
        # fmt:on

        if config.workflow.use_aroma:  # ICA-AROMA workflow
            from .confounds import init_ica_aroma_wf

            ica_aroma_wf = init_ica_aroma_wf(
                mem_gb=mem_gb["resampled"],
                metadata=metadata,
                omp_nthreads=omp_nthreads,
                err_on_aroma_warn=config.workflow.aroma_err_on_warn,
                aroma_melodic_dim=config.workflow.aroma_melodic_dim,
                name="ica_aroma_wf",
            )

            join = pe.Node(
                niu.Function(output_names=["out_file"], function=_to_join),
                name="aroma_confounds",
            )

            mrg_conf_metadata = pe.Node(
                niu.Merge(2),
                name="merge_confound_metadata",
                run_without_submitting=True,
            )
            mrg_conf_metadata2 = pe.Node(
                DictMerge(),
                name="merge_confound_metadata2",
                run_without_submitting=True,
            )
            # fmt:off
            workflow.disconnect([
                (bold_confounds_wf, outputnode, [
                    ("outputnode.confounds_file", "confounds"),
                ]),
                (bold_confounds_wf, outputnode, [
                    ("outputnode.confounds_metadata", "confounds_metadata"),
                ]),
            ])
            workflow.connect([
                (inputnode, ica_aroma_wf, [("bold_file", "inputnode.name_source")]),
                (bold_hmc_wf, ica_aroma_wf, [
                    ("outputnode.movpar_file", "inputnode.movpar_file"),
                ]),
                (initial_boldref_wf, ica_aroma_wf, [
                    ("outputnode.skip_vols", "inputnode.skip_vols"),
                ]),
                (bold_confounds_wf, join, [("outputnode.confounds_file", "in_file")]),
                (bold_confounds_wf, mrg_conf_metadata, [
                    ("outputnode.confounds_metadata", "in1"),
                ]),
                (ica_aroma_wf, join, [("outputnode.aroma_confounds", "join_file")]),
                (ica_aroma_wf, mrg_conf_metadata, [("outputnode.aroma_metadata", "in2")]),
                (mrg_conf_metadata, mrg_conf_metadata2, [("out", "in_dicts")]),
                (ica_aroma_wf, outputnode, [
                    ("outputnode.aroma_noise_ics", "aroma_noise_ics"),
                    ("outputnode.melodic_mix", "melodic_mix"),
                    ("outputnode.nonaggr_denoised_file", "nonaggr_denoised_file"),
                ]),
                (join, outputnode, [("out_file", "confounds")]),
                (mrg_conf_metadata2, outputnode, [("out_dict", "confounds_metadata")]),
                (bold_std_trans_wf, ica_aroma_wf, [
                    ("outputnode.bold_std", "inputnode.bold_std"),
                    ("outputnode.bold_mask_std", "inputnode.bold_mask_std"),
                    ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ]),
            ])
            # fmt:on

    # SURFACES ##################################################################################
    # Freesurfer
    freesurfer_spaces = spaces.get_fs_spaces()
    if freesurfer and freesurfer_spaces:
        config.loggers.workflow.debug("Creating BOLD surface-sampling workflow.")
        bold_surf_wf = init_bold_surf_wf(
            mem_gb=mem_gb["resampled"],
            surface_spaces=freesurfer_spaces,
            medial_surface_nan=config.workflow.medial_surface_nan,
            name="bold_surf_wf",
        )
        # fmt:off
        workflow.connect([
            (inputnode, bold_surf_wf, [
                ("subjects_dir", "inputnode.subjects_dir"),
                ("subject_id", "inputnode.subject_id"),
                ("t1w2fsnative_xfm", "inputnode.t1w2fsnative_xfm"),
            ]),
            (bold_t1_trans_wf, bold_surf_wf, [("outputnode.bold_t1", "inputnode.source_file")]),
            (bold_surf_wf, outputnode, [("outputnode.surfaces", "surfaces")]),
            (bold_surf_wf, func_derivatives_wf, [("outputnode.target", "inputnode.surf_refs")]),
        ])
        # fmt:on

        # CIFTI output
        if config.workflow.cifti_output:
            from .resampling import init_bold_grayords_wf

            bold_grayords_wf = init_bold_grayords_wf(
                grayord_density=config.workflow.cifti_output,
                mem_gb=mem_gb["resampled"],
                repetition_time=metadata["RepetitionTime"],
            )

            # fmt:off
            workflow.connect([
                (inputnode, bold_grayords_wf, [("subjects_dir", "inputnode.subjects_dir")]),
                (bold_std_trans_wf, bold_grayords_wf, [
                    ("outputnode.bold_std", "inputnode.bold_std"),
                    ("outputnode.spatial_reference", "inputnode.spatial_reference"),
                ]),
                (bold_surf_wf, bold_grayords_wf, [
                    ("outputnode.surfaces", "inputnode.surf_files"),
                    ("outputnode.target", "inputnode.surf_refs"),
                ]),
                (bold_grayords_wf, outputnode, [
                    ("outputnode.cifti_bold", "bold_cifti"),
                    ("outputnode.cifti_variant", "cifti_variant"),
                    ("outputnode.cifti_metadata", "cifti_metadata"),
                    ("outputnode.cifti_density", "cifti_density"),
                ]),
            ])
            # fmt:on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        carpetplot_wf = init_carpetplot_wf(
            mem_gb=mem_gb["resampled"],
            metadata=metadata,
            cifti_output=config.workflow.cifti_output,
            name="carpetplot_wf",
        )

        if config.workflow.cifti_output:
            workflow.connect(
                bold_grayords_wf,
                "outputnode.cifti_bold",
                carpetplot_wf,
                "inputnode.cifti_bold",
            )
        else:
            # Xform to "MNI152NLin2009cAsym" is always computed.
            carpetplot_select_std = pe.Node(
                KeySelect(fields=["std2anat_xfm"], key="MNI152NLin2009cAsym"),
                name="carpetplot_select_std",
                run_without_submitting=True,
            )

            # fmt:off
            workflow.connect([
                (inputnode, carpetplot_select_std, [("std2anat_xfm", "std2anat_xfm"),
                                                    ("template", "keys")]),
                (carpetplot_select_std, carpetplot_wf, [
                    ("std2anat_xfm", "inputnode.std2anat_xfm"),
                ]),
                (bold_final, carpetplot_wf, [
                    ("bold", "inputnode.bold"),
                    ("mask", "inputnode.bold_mask"),
                ]),
                (bold_reg_wf, carpetplot_wf, [
                    ("outputnode.itk_t1_to_bold", "inputnode.t1_bold_xform"),
                ]),
            ])
            # fmt:on

        # fmt:off
        workflow.connect([
            (bold_confounds_wf, carpetplot_wf, [
                ("outputnode.confounds_file", "inputnode.confounds_file"),
            ])
        ])
        # fmt:on

    # REPORTING ############################################################
    ds_report_summary = pe.Node(
        DerivativesDataSink(
            desc="summary", datatype="figures", dismiss_entities=("echo",)
        ),
        name="ds_report_summary",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    ds_report_validation = pe.Node(
        DerivativesDataSink(
            desc="validation", datatype="figures", dismiss_entities=("echo",)
        ),
        name="ds_report_validation",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (summary, ds_report_summary, [("out_report", "in_file")]),
        (initial_boldref_wf, ds_report_validation, [("outputnode.validation_report", "in_file")]),
    ])
    # fmt:on

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split(".")[-1].startswith("ds_report"):
            workflow.get_node(node).inputs.base_directory = fmriprep_dir
            workflow.get_node(node).inputs.source_file = ref_file

    if not has_fieldmap:
        # Finalize workflow without SDC connections
        summary.inputs.distortion_correction = "None"

        # Resample in native space in just one shot
        bold_bold_trans_wf = init_bold_preproc_trans_wf(
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            use_compression=not config.execution.low_mem,
            use_fieldwarp=False,
            name="bold_bold_trans_wf",
        )
        bold_bold_trans_wf.inputs.inputnode.fieldwarp = "identity"

        # fmt:off
        workflow.connect([
            # Connect bold_bold_trans_wf
            (bold_source, bold_bold_trans_wf, [("out", "inputnode.name_source")]),
            (bold_split, bold_bold_trans_wf, [("out_files", "inputnode.bold_file")]),
            (bold_hmc_wf, bold_bold_trans_wf, [
                ("outputnode.xforms", "inputnode.hmc_xforms"),
            ]),
        ])

        workflow.connect([
            (bold_bold_trans_wf, bold_final, [("outputnode.bold", "bold")]),
            (bold_bold_trans_wf, final_boldref_wf, [
                ("outputnode.bold", "inputnode.bold_file"),
            ]),
        ] if not multiecho else [
            (initial_boldref_wf, bold_t2s_wf, [
                ("outputnode.bold_mask", "inputnode.bold_mask"),
            ]),
            (bold_bold_trans_wf, join_echos, [
                ("outputnode.bold", "bold_files"),
            ]),
            (join_echos, final_boldref_wf, [
                ("bold_files", "inputnode.bold_file"),
            ]),
        ])
        # fmt:on
        return workflow

    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from niworkflows.interfaces.utility import KeySelect
    from sdcflows.workflows.apply.registration import init_coeff2epi_wf
    from sdcflows.workflows.apply.correction import init_unwarp_wf

    coeff2epi_wf = init_coeff2epi_wf(
        debug="fieldmaps" in config.execution.debug,
        omp_nthreads=config.nipype.omp_nthreads,
        write_coeff=True,
    )
    unwarp_wf = init_unwarp_wf(
        debug="fieldmaps" in config.execution.debug,
        omp_nthreads=config.nipype.omp_nthreads,
    )
    unwarp_wf.inputs.inputnode.metadata = metadata

    output_select = pe.Node(
        KeySelect(fields=["fmap", "fmap_ref", "fmap_coeff", "fmap_mask", "sdc_method"]),
        name="output_select",
        run_without_submitting=True,
    )
    output_select.inputs.key = estimator_key[0]
    if len(estimator_key) > 1:
        config.loggers.workflow.warning(
            f"Several fieldmaps <{', '.join(estimator_key)}> are "
            f"'IntendedFor' <{bold_file}>, using {estimator_key[0]}"
        )

    sdc_report = pe.Node(
        SimpleBeforeAfter(
            before_label="Distorted",
            after_label="Corrected",
            dismiss_affine=True,
        ),
        name="sdc_report",
        mem_gb=0.1,
    )

    ds_report_sdc = pe.Node(
        DerivativesDataSink(
            base_directory=fmriprep_dir,
            desc="sdc",
            suffix="bold",
            datatype="figures",
            dismiss_entities=("echo",),
        ),
        name="ds_report_sdc",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, output_select, [("fmap", "fmap"),
                                    ("fmap_ref", "fmap_ref"),
                                    ("fmap_coeff", "fmap_coeff"),
                                    ("fmap_mask", "fmap_mask"),
                                    ("sdc_method", "sdc_method"),
                                    ("fmap_id", "keys")]),
        (output_select, coeff2epi_wf, [
            ("fmap_ref", "inputnode.fmap_ref"),
            ("fmap_coeff", "inputnode.fmap_coeff"),
            ("fmap_mask", "inputnode.fmap_mask")]),
        (output_select, summary, [("sdc_method", "distortion_correction")]),
        (initial_boldref_wf, coeff2epi_wf, [
            ("outputnode.ref_image", "inputnode.target_ref"),
            ("outputnode.bold_mask", "inputnode.target_mask")]),
        (coeff2epi_wf, unwarp_wf, [
            ("outputnode.fmap_coeff", "inputnode.fmap_coeff")]),
        (bold_hmc_wf, unwarp_wf, [
            ("outputnode.xforms", "inputnode.hmc_xforms")]),
        (initial_boldref_wf, sdc_report, [
            ("outputnode.ref_image", "before")]),
        (bold_split, unwarp_wf, [
            ("out_files", "inputnode.distorted")]),
        (final_boldref_wf, sdc_report, [
            ("outputnode.ref_image", "after"),
            ("outputnode.bold_mask", "wm_seg")]),
        (inputnode, ds_report_sdc, [("bold_file", "source_file")]),
        (sdc_report, ds_report_sdc, [("out_report", "in_file")]),

    ])
    # fmt:on

    if not multiecho:
        # fmt:off
        workflow.connect([
            (unwarp_wf, bold_final, [("outputnode.corrected", "bold")]),
            # remaining workflow connections
            (unwarp_wf, final_boldref_wf, [
                ("outputnode.corrected", "inputnode.bold_file"),
            ]),
            (unwarp_wf, bold_t1_trans_wf, [
                # TEMPORARY: For the moment we can't use frame-wise fieldmaps
                (("outputnode.fieldwarp", pop_file), "inputnode.fieldwarp"),
            ]),
            (unwarp_wf, bold_std_trans_wf, [
                # TEMPORARY: For the moment we can't use frame-wise fieldmaps
                (("outputnode.fieldwarp", pop_file), "inputnode.fieldwarp"),
            ]),
        ])
        # fmt:on
        return workflow

    # Finalize connections if ME-EPI
    join_sdc_echos = pe.JoinNode(
        niu.IdentityInterface(
            fields=[
                "fieldmap",
                "fieldwarp",
                "corrected",
                "corrected_ref",
                "corrected_mask",
            ]
        ),
        joinsource="echo_index",
        joinfield=[
            "fieldmap",
            "fieldwarp",
            "corrected",
            "corrected_ref",
            "corrected_mask",
        ],
        name="join_sdc_echos",
    )

    def _dpop(list_of_lists):
        return list_of_lists[0][0]

    # fmt:off
    workflow.connect([
        (unwarp_wf, join_echos, [
            ("outputnode.corrected", "bold_files"),
        ]),
        (unwarp_wf, join_sdc_echos, [
            ("outputnode.fieldmap", "fieldmap"),
            ("outputnode.fieldwarp", "fieldwarp"),
            ("outputnode.corrected", "corrected"),
            ("outputnode.corrected_ref", "corrected_ref"),
            ("outputnode.corrected_mask", "corrected_mask"),
        ]),
        # remaining workflow connections
        (join_sdc_echos, final_boldref_wf, [
            ("corrected", "inputnode.bold_file"),
        ]),
        (join_sdc_echos, bold_t2s_wf, [
            (("corrected_mask", pop_file), "inputnode.bold_mask"),
        ]),
        (join_sdc_echos, bold_t1_trans_wf, [
            # TEMPORARY: For the moment we can't use frame-wise fieldmaps
            (("fieldwarp", _dpop), "inputnode.fieldwarp"),
        ]),
        (join_sdc_echos, bold_std_trans_wf, [
            # TEMPORARY: For the moment we can't use frame-wise fieldmaps
            (("fieldwarp", _dpop), "inputnode.fieldwarp"),
        ]),
    ])
    # fmt:on

    return workflow


def _create_mem_gb(bold_fname):
    bold_size_gb = os.path.getsize(bold_fname) / (1024 ** 3)
    bold_tlen = nb.load(bold_fname).shape[-1]
    mem_gb = {
        "filesize": bold_size_gb,
        "resampled": bold_size_gb * 4,
        "largemem": bold_size_gb * (max(bold_tlen / 100, 1.0) + 4),
    }

    return bold_tlen, mem_gb


def _get_wf_name(bold_fname):
    """
    Derive the workflow name for supplied BOLD file.

    >>> _get_wf_name("/completely/made/up/path/sub-01_task-nback_bold.nii.gz")
    'func_preproc_task_nback_wf'
    >>> _get_wf_name("/completely/made/up/path/sub-01_task-nback_run-01_echo-1_bold.nii.gz")
    'func_preproc_task_nback_run_01_echo_1_wf'

    """
    from nipype.utils.filemanip import split_filename

    fname = split_filename(bold_fname)[1]
    fname_nosub = "_".join(fname.split("_")[1:])
    name = "func_preproc_" + fname_nosub.replace(".", "_").replace(" ", "").replace(
        "-", "_"
    ).replace("_bold", "_wf")

    return name


def _to_join(in_file, join_file):
    """Join two tsv files if the join_file is not ``None``."""
    from niworkflows.interfaces.utility import JoinTSVColumns

    if join_file is None:
        return in_file
    res = JoinTSVColumns(in_file=in_file, join_file=join_file).run()
    return res.outputs.out_file


def extract_entities(file_list):
    """
    Return a dictionary of common entities given a list of files.

    Examples
    --------
    >>> extract_entities("sub-01/anat/sub-01_T1w.nii.gz")
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(["sub-01/anat/sub-01_T1w.nii.gz"] * 2)
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(["sub-01/anat/sub-01_run-1_T1w.nii.gz",
    ...                   "sub-01/anat/sub-01_run-2_T1w.nii.gz"])
    {'subject': '01', 'run': [1, 2], 'suffix': 'T1w', 'datatype': 'anat',
     'extension': '.nii.gz'}

    """
    from collections import defaultdict
    from bids.layout import parse_file_entities

    entities = defaultdict(list)
    for e, v in [
        ev_pair
        for f in listify(file_list)
        for ev_pair in parse_file_entities(f).items()
    ]:
        entities[e].append(v)

    def _unique(inlist):
        inlist = sorted(set(inlist))
        if len(inlist) == 1:
            return inlist[0]
        return inlist

    return {k: _unique(v) for k, v in entities.items()}


def get_img_orientation(imgf):
    """Return the image orientation as a string"""
    img = nb.load(imgf)
    return "".join(nb.aff2axcodes(img.affine))
