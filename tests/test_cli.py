#!/usr/bin/env python

"""Tests for `cli` module."""
import pytest
from spaceprime.cli import sci_notation_int


def test_sci_notation_int_plain_int():
    assert sci_notation_int("1000") == 1000
    assert isinstance(sci_notation_int("1000"), int)


def test_sci_notation_int_scientific_notation():
    assert sci_notation_int("1e6") == 1000000
    assert isinstance(sci_notation_int("1e6"), int)


def test_sci_notation_int_scientific_notation_small_exponent():
    assert sci_notation_int("1e3") == 1000
    assert sci_notation_int("2.5e4") == 25000


def test_sci_notation_int_raises_on_non_numeric():
    with pytest.raises((ValueError, TypeError)):
        sci_notation_int("abc")
