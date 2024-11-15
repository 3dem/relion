from pathlib import Path

import napari

from .._qt.tilt_image_excluder import TiltImageExcluderWidget
from .._metadata_models.relion.tilt_series_set import RlnTiltSeriesSet


class RelionTiltImageExcluderWidget(TiltImageExcluderWidget):
    def __init__(
            self,
            viewer: napari.Viewer,
            relion_metadata: RlnTiltSeriesSet,
            cache_size: int,
            output_global_metadata_file: Path,
            *args,
            **kwargs
    ):
        # sort relion metadata by tilt angle first
        relion_metadata = relion_metadata.copy(deep=True)
        for tilt_series in relion_metadata.tilt_series:
            tilt_series.data = tilt_series.data.sort_values(by='rlnTomoNominalStageTiltAngle')

        super().__init__(
            viewer=viewer,
            tilt_series=relion_metadata.as_gui_model(),
            cache_size=cache_size,
            save_button_label='save tilt-series STAR file',
            *args,
            **kwargs
        )
        self.relion_metadata = relion_metadata
        self.output_filename = output_global_metadata_file

    def save_output(self):
        relion_metadata = self.relion_metadata.copy(deep=True)
        for tilt_series in relion_metadata.tilt_series:
            selected_image_idx = self.selected_tilt_images[tilt_series.name]
            tilt_series.data = tilt_series.data[selected_image_idx]
        relion_metadata.write_star_file(self.output_filename)
        self._post_successful_save()
