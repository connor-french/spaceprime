"""Tests for the spaceprime CLI."""

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
    assert "Loading packages..." not in captured.out


def test_cli_version_exits_before_loading_packages(monkeypatch, capsys):
    """Test that --version exits before printing 'Loading packages...'."""
    monkeypatch.setattr("sys.argv", ["spaceprime", "--version"])

    with pytest.raises(SystemExit):
        cli.main()

    captured = capsys.readouterr()
    assert "Loading packages..." not in captured.out


def test_cli_loading_packages_message_printed(monkeypatch, capsys):
    """Test that 'Loading packages...' is printed when proceeding past arg parsing."""
    monkeypatch.setattr("sys.argv", ["spaceprime", "-r", "fake.tif", "-co", "fake.csv"])

    # An error is expected before a simulation runs (raster file not found or
    # a validation error), but the loading message must appear first.
    with pytest.raises((FileNotFoundError, TypeError, ValueError)):
        cli.main()

    captured = capsys.readouterr()
    assert "Loading packages..." in captured.out
