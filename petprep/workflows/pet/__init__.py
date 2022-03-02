# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""

Pre-processing fMRI - BOLD signal workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: fmriprep.workflows.bold.base
.. automodule:: fmriprep.workflows.bold.hmc
.. automodule:: fmriprep.workflows.bold.registration
.. automodule:: fmriprep.workflows.bold.resampling
.. automodule:: fmriprep.workflows.bold.confounds


"""

from .base import init_pet_preproc_wf
from .hmc import init_pet_hmc_wf
from .registration import (
    init_pet_t1_trans_wf,
    init_pet_reg_wf,
)
from .resampling import (
    init_pet_std_trans_wf,
    init_pet_surf_wf,
    init_pet_preproc_trans_wf,
)

from .confounds import (
    init_pet_confs_wf
)

__all__ = [
    'init_pet_confs_wf',
    'init_pet_hmc_wf',
    'init_pet_std_trans_wf',
    'init_pet_preproc_trans_wf',
    'init_pet_reg_wf',
    'init_pet_surf_wf',
    'init_pet_t1_trans_wf',
    'init_pet_preproc_wf'
]
