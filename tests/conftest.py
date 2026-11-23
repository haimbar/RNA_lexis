"""Shared pytest fixtures."""

import pytest


@pytest.fixture(autouse=True)
def _never_open_files_in_a_viewer(monkeypatch):
    """The plot functions open every saved figure in the system viewer
    (Preview, the default PDF reader, ...).  The tests save into temporary
    directories that are deleted straight away, so without this the viewer is
    left showing "The file ... couldn't be opened because there is no such
    file".  Replace the openers with no-ops for every test."""
    import rna_lexis.plots as plots
    monkeypatch.setattr(plots, "open_file_with_default_software", lambda *a, **k: None)
    monkeypatch.setattr(plots, "open_pdf", lambda *a, **k: None)
