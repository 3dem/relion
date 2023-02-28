from pathlib import Path

import napari
import typer

from .._gui.components.tomogram_browser import TomogramBrowserWidget
from .._metadata_models.relion.tilt_series_set import TiltSeriesSet

cli = typer.Typer(add_completion=False)


@cli.command(name='relion_tomo_view_tomograms')
def tomogram_viewer(
        tilt_series_star_file: Path = typer.Option(...),
        cache_size: int = typer.Option(3, help='number of cached tomograms')
):
    viewer = napari.Viewer(ndisplay=3)
    relion_metadata = TiltSeriesSet.from_star_file(tilt_series_star_file)
    gui_model = relion_metadata.as_gui_model()
    dock_widget = TomogramBrowserWidget(viewer, gui_model, cache_size=cache_size)
    viewer.window.add_dock_widget(
        dock_widget, name='RELION tomogram viewer', area='left'
    )
    viewer.axes.visible = True
    viewer.axes.labels = False
    viewer.camera.angles = (-15, 30, 150)
    viewer.camera.zoom *= 0.7
    viewer.text_overlay.text = """
    keyboard shortcuts
    - '[' and ']' for previous/next tomogram
    - 'x'/'y'/'z'/'o' to reorient the plane (around the cursor)
    
    mouse controls
    - 'Shift' + click + drag to move the plane along its normal
    """
    viewer.text_overlay.visible = True
    napari.run()
