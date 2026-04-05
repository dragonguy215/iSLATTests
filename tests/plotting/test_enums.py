# -*- coding: utf-8 -*-
"""Tests for GapMode and XScaling enums."""

import pytest

from iSLAT.Modules.Plotting import GapMode, XScaling


class TestGapMode:
    def test_values(self):
        assert GapMode.CONNECT.value == "connect"
        assert GapMode.SKIP.value == "skip"

    def test_from_string(self):
        assert GapMode("connect") is GapMode.CONNECT
        assert GapMode("skip") is GapMode.SKIP

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            GapMode("invalid")


class TestXScaling:
    def test_values(self):
        assert XScaling.WAVELENGTH.value == "wavelength"
        assert XScaling.DATA_DENSITY.value == "data_density"

    def test_from_string(self):
        assert XScaling("wavelength") is XScaling.WAVELENGTH
        assert XScaling("data_density") is XScaling.DATA_DENSITY

    def test_invalid_raises(self):
        with pytest.raises(ValueError):
            XScaling("bad")
