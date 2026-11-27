# Import decorated command line interface functions first
from .imod import fiducials_cli
from .imod import patch_tracking_cli
try:
    from .aretomo import aretomo_cli
except ImportError:
    # AreTomo support is optional (it's a closed-source binary tied to old
    # CUDA versions) - fall back gracefully so the IMOD-based backends
    # (fiducials/patch tracking) still work without it.
    aretomo_cli = None
from ._job_utils import write_global_output

# Then import the cli, it will be decorated with all programs (subcommands)
from ._cli import cli
