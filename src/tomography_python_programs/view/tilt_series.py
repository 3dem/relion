from pathlib import Path

import typer

from .._qt.components.tilt_series_browser import TiltSeriesBrowserWidget
from .._metadata_models.relion.tilt_series_set import RlnTiltSeriesSet
from ..view._cli import cli

COMMAND_NAME = 'tilt-series'


@cli.command(name=COMMAND_NAME, no_args_is_help=True)
def view_tilt_series(
    tilt_series_star_file: Path = typer.Option(...),
    cache_size: int = typer.Option(5, help='number of cached tilt-series')
):
    import napari  # local to speedup CLI
    viewer = napari.Viewer(ndisplay=3)
    relion_metadata = RlnTiltSeriesSet.from_star_file(tilt_series_star_file)
    for tilt_series in relion_metadata.tilt_series:
        tilt_series.data = tilt_series.data.sort_values(
            by='rlnTomoNominalStageTiltAngle')
    gui_model = relion_metadata.as_gui_model()
    dock_widget = TiltSeriesBrowserWidget(viewer, gui_model, cache_size=cache_size)
    viewer.window.add_dock_widget(
        dock_widget, name='RELION tilt-series viewer', area='left'
    )
    viewer.text_overlay.text = """
    keyboard shortcuts
    - '[' and ']' for previous/next tilt-series
    """
    viewer.text_overlay.visible = True
    napari.run()
