# Import decorated batch alignment functions first
from .imod import batch_fiducials as _batch_imod_fiducials
from .imod import batch_patch_tracking as _batch_imod_patch_tracking
from .aretomo import batch_aretomo as _batch_aretomo
from ._job_utils import write_aligned_tilt_series_star_file as _write_global_output

# Then import the cli, it will be decorated with all programs
from ._cli import cli
