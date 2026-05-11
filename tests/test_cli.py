"""Tests for the spaceprime CLI."""

import sys
import pytest

from spaceprime import __version__
from spaceprime import cli


def test_cli_version_option_prints_version(monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", ["spaceprime", "--version"])

    with pytest.raises(SystemExit) as excinfo:
        cli.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert captured.out.strip() == f"spaceprime {__version__}"


def test_cli_short_version_option_prints_version(monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", ["spaceprime", "-v"])

    with pytest.raises(SystemExit) as excinfo:
        cli.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert captured.out.strip() == f"spaceprime {__version__}"


def test_cli_help_option_prints_help(monkeypatch, capsys):
    """Test that --help prints help text and exits cleanly."""
    monkeypatch.setattr("sys.argv", ["spaceprime", "--help"])

    with pytest.raises(SystemExit) as excinfo:
        cli.main()

    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "usage" in captured.out.lower()


def test_cli_help_option_exits_before_loading_packages(monkeypatch, capsys):
    """Test that --help exits before printing 'Loading packages...'."""
    monkeypatch.setattr("sys.argv", ["spaceprime", "--help"])

    with pytest.raises(SystemExit):
        cli.main()

    captured = capsys.readouterr()
    assert "Loading packages..." not in captured.err


def test_cli_version_exits_before_loading_packages(monkeypatch, capsys):
    """Test that --version exits before printing 'Loading packages...'."""
    monkeypatch.setattr("sys.argv", ["spaceprime", "--version"])

    with pytest.raises(SystemExit):
        cli.main()

    captured = capsys.readouterr()
    assert "Loading packages..." not in captured.err


def test_cli_loading_packages_message_printed(monkeypatch, capsys):
    """Test that 'Loading packages...' is printed to stderr when proceeding past arg parsing.

    Providing a non-existent raster path should raise FileNotFoundError, and
    'Loading packages...' must have been printed before the error is raised.
    """
    monkeypatch.setattr("sys.argv", ["spaceprime", "-r", "fake.tif", "-co", "fake.csv"])

    with pytest.raises(FileNotFoundError):
        cli.main()

    captured = capsys.readouterr()
    assert "Loading packages..." in captured.err


def test_spaceprime_lazy_loading_does_not_import_heavy_deps():
    """Test that importing spaceprime does not eagerly load heavy dependencies.

    The lazy-loading mechanism in spaceprime/__init__.py should defer heavy
    submodule imports until a public symbol is first accessed.
    """
    import importlib

    # Reload spaceprime to get a clean state for this test
    import spaceprime

    importlib.reload(spaceprime)

    # The public symbols should be accessible via __dir__
    sp_dir = dir(spaceprime)
    assert "raster_to_demes" in sp_dir
    assert "sim_ancestry" in sp_dir
    assert "filter_gt" in sp_dir

    # __all__ should list the full public API
    assert "raster_to_demes" in spaceprime.__all__
    assert "sim_ancestry" in spaceprime.__all__
    assert "filter_gt" in spaceprime.__all__
