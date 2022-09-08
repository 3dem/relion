from pathlib import Path

import napari
import typer
from napari.settings import get_settings
from napari.utils.notifications import NotificationSeverity

from .relion_tilt_image_excluder import RelionTiltImageExcluderWidget
from ..metadata_model import RelionTiltSeriesSet
from ..utils.relion import relion_pipeline_job

cli = typer.Typer(add_completion=False)


@cli.command(name='relion_tomo_exclude_tilt_images')
@relion_pipeline_job
def exclude_tilt_images_cli(
        tilt_series_star_file: Path = typer.Option(...),
        output_directory: Path = typer.Option(...),
        cache_size: int = typer.Option(5, help='number of cached tilt-series')
):
    viewer = napari.Viewer()
    settings = get_settings()
    settings.application.gui_notification_level = NotificationSeverity.ERROR
    relion_metadata = RelionTiltSeriesSet.from_star_file(tilt_series_star_file)
    dock_widget = RelionTiltImageExcluderWidget(
        viewer=viewer,
        relion_metadata=relion_metadata,
        cache_size=cache_size,
        output_global_metadata_file=output_directory / 'selected_tilt_series.star'
    )
    viewer.window.add_dock_widget(dock_widget, name='RELION tilt-image excluder')
    viewer.text_overlay.text = """
    keyboard shortcuts
    - '[' and ']' for previous/next tilt-series
    - '↑' and '↓' for previous/next image
    - 'PageUp' and 'PageDown' for first/last image
    - 'Enter' to select/deselect an image
    """
    viewer.text_overlay.visible = True
    napari.run()
