from tomography_python_programs.metadata_model import RelionTiltSeriesSet, RelionTiltSeries


def test_read_single_tilt_series_metadata(tilt_series_star_file):
    """Simple parsing test."""
    tilt_series = RelionTiltSeries.from_star_file(tilt_series_star_file)
    assert tilt_series.data.shape == (41, 29)


def test_read_global_metadata(global_relion_metadata_file):
    """Simple parsing test."""
    relion_metadata = RelionTiltSeriesSet.from_star_file(global_relion_metadata_file)
    assert isinstance(relion_metadata, RelionTiltSeriesSet)
    assert len(relion_metadata.tilt_series) == 5
    assert len(relion_metadata.global_data) == 5
    assert relion_metadata.global_data.shape == (5, 13)


def test_global_metadata_correctly_ordered(relion_metadata):
    """Ordering of tilt-series should be consistent across metadata"""
    names_a = relion_metadata.global_data['rlnTomoName']
    names_b = [tilt_series.name for tilt_series in relion_metadata.tilt_series]
    assert all([a == b for a, b in zip(names_a, names_b)])


def test_gui_data_model_from_relion_metadata(relion_metadata):
    """Test correct generation of simple data model compatible with GUI components."""
    gui_model = relion_metadata.as_gui_model()
    assert isinstance(gui_model, dict)
    assert len(gui_model) == 5
    assert all([
        a == b
        for a, b
        in zip(relion_metadata.global_data['rlnTomoName'], gui_model.keys())
    ]
    )


def test_gui_data_model_includes_tomogram_file(relion_metadata):
    """Test that gui model generated from relion metadata references tomograms."""
    gui_model = relion_metadata.as_gui_model()
    tilt_series = gui_model[list(gui_model.keys())[0]]
    assert tilt_series.tomogram_file is not None
