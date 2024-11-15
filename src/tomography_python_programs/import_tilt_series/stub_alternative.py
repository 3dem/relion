"""Stub alternative CLI command for import_tilt_series job."""

from ._cli import cli


@cli.command(name='alternative', hidden=True)
def stub():
    """Stub for an alternative to SerialEM tomo metadata import_tilt_series."""
    pass
