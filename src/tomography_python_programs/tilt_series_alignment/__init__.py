# Import decorated command line interface functions first
from .imod import fiducials_cli
from .imod import patch_tracking_cli
from .aretomo import aretomo_cli
from ._job_utils import write_global_output

# Then import the cli, it will be decorated with all programs (subcommands)
from ._cli import cli
