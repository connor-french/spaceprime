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
