from pathlib import Path

import napari
import typer

from tomography_python_programs.gui.components.tilt_series_browser import TiltSeriesBrowserWidget
from ..metadata_model import RelionTiltSeriesSet

cli = typer.Typer(add_completion=False)


@cli.command(name='relion_tomo_view_tilt_series')
def tilt_series_viewer(
        tilt_series_star_file: Path = typer.Option(...),
        cache_size: int = typer.Option(5, help='number of cached tilt-series')
):
    viewer = napari.Viewer(ndisplay=2)
    relion_metadata = RelionTiltSeriesSet.from_star_file(tilt_series_star_file)
    for tilt_series in relion_metadata.tilt_series:
        tilt_series.data = tilt_series.data.sort_values(by='rlnTomoNominalStageTiltAngle')
    gui_model = RelionTiltSeriesSet.as_gui_model()
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
