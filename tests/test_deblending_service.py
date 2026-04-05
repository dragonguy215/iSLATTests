# -*- coding: utf-8 -*-
"""
Unit tests for iSLAT.Modules.DataProcessing.DeblendingService.

Tests component extraction, formatting, and display.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from iSLAT.Modules.DataProcessing.DeblendingService import DeblendingService


# ============================================================================
# Construction
# ============================================================================

class TestDeblendingServiceInit:
    """Tests for DeblendingService initialization."""

    def test_can_instantiate(self):
        svc = DeblendingService()
        assert svc is not None


# ============================================================================
# extract_deblended_components
# ============================================================================

class TestExtractDeblendedComponents:
    """Tests for extracting deblended component info from fit results."""

    @staticmethod
    def _mock_molecule_line(lam, lev_up='0|10', lev_low='0|9',
                             a_stein=0.01, e_up=4586.0, e_low=3379.0,
                             g_up=21, g_low=19):
        """Create a mock MoleculeLine."""
        ml = MagicMock()
        ml.lam = lam
        ml.get_dict.return_value = {
            'lam': lam, 'lev_up': lev_up, 'lev_low': lev_low,
            'a_stein': a_stein, 'e_up': e_up, 'e_low': e_low,
            'g_up': g_up, 'g_low': g_low,
        }
        return ml

    def test_single_component(self):
        svc = DeblendingService()
        line_params = {
            'component_0': {
                'center': 15.001,
                'center_stderr': 0.0001,
                'fwhm': 12.5,
                'fwhm_stderr': 0.5,
                'area': 1e-15,
                'area_stderr': 1e-17,
            }
        }
        ml = self._mock_molecule_line(15.0)
        line_info = [(ml, 5e-16, 0.3)]

        components = svc.extract_deblended_components(line_params, line_info, 'H2O')

        assert len(components) == 1
        comp = components[0]
        assert comp['index'] == 0
        assert comp['center'] == 15.001
        assert comp['molecule_name'] == 'H2O'
        assert comp['lam'] == 15.0
        assert comp['lev_up'] == '0|10'
        assert 'doppler' in comp

    def test_two_components(self):
        svc = DeblendingService()
        line_params = {
            'component_0': {
                'center': 15.001, 'center_stderr': 0.0001,
                'fwhm': 12.5, 'fwhm_stderr': 0.5,
                'area': 1e-15, 'area_stderr': 1e-17,
            },
            'component_1': {
                'center': 15.05, 'center_stderr': 0.0002,
                'fwhm': 10.0, 'fwhm_stderr': 0.8,
                'area': 5e-16, 'area_stderr': 2e-17,
            },
        }
        ml1 = self._mock_molecule_line(15.0)
        ml2 = self._mock_molecule_line(15.04, lev_up='0|8', lev_low='0|7')
        line_info = [(ml1, 5e-16, 0.3), (ml2, 3e-16, 0.1)]

        components = svc.extract_deblended_components(line_params, line_info, 'H2O')

        assert len(components) == 2
        assert components[0]['index'] == 0
        assert components[1]['index'] == 1
        assert components[1]['center'] == 15.05

    def test_no_components(self):
        svc = DeblendingService()
        line_params = {'center': 15.0}  # No component_N keys
        components = svc.extract_deblended_components(line_params, [], 'H2O')
        assert components == []

    def test_more_components_than_line_info(self):
        """Extra components should still be extracted (without molecular data)."""
        svc = DeblendingService()
        line_params = {
            'component_0': {
                'center': 15.0, 'center_stderr': 0.0001,
                'fwhm': 12.0, 'fwhm_stderr': 0.5,
                'area': 1e-15, 'area_stderr': 1e-17,
            },
            'component_1': {
                'center': 15.1, 'center_stderr': 0.0002,
                'fwhm': 10.0, 'fwhm_stderr': 0.8,
                'area': 5e-16, 'area_stderr': 2e-17,
            },
        }
        # Only one line_info entry
        ml = self._mock_molecule_line(15.0)
        line_info = [(ml, 5e-16, 0.3)]

        components = svc.extract_deblended_components(line_params, line_info, 'CO')
        assert len(components) == 2
        # First has molecular data, second doesn't
        assert components[0].get('lam') == 15.0
        assert 'lam' not in components[1] or components[1].get('lam', 0.0) == 0.0


# ============================================================================
# format_component_for_save
# ============================================================================

class TestFormatComponentForSave:
    """Tests for formatting component dict for file output."""

    def test_basic_format(self):
        svc = DeblendingService()
        comp = {
            'molecule_name': 'H2O',
            'lev_up': '0|10', 'lev_low': '0|9',
            'lam': 15.0, 'tau': 0.3, 'intens': 5e-16,
            'a_stein': 0.01, 'e_up': 4586.0, 'e_low': 3379.0,
            'g_up': 21, 'g_low': 19,
            'center': 15.001, 'center_stderr': 0.0001,
            'fwhm': 12.5, 'fwhm_stderr': 0.5,
            'area': 1e-15, 'area_stderr': 1e-17,
            'doppler': 20.0,
        }
        result = svc.format_component_for_save(comp)

        assert result['species'] == 'H2O'
        assert result['Flux_fit'] == 1e-15
        assert result['Err_fit'] == 1e-17
        assert result['FWHM_fit'] == 12.5
        assert result['Centr_fit'] == 15.001
        assert result['Doppler'] == 20.0

    def test_missing_molecule_name_uses_unknown(self):
        svc = DeblendingService()
        comp = {
            'center': 15.0, 'center_stderr': 0.0,
            'fwhm': 10.0, 'fwhm_stderr': 0.0,
            'area': 1e-15, 'area_stderr': 0.0,
        }
        result = svc.format_component_for_save(comp)
        assert result['species'] == 'Unknown'


# ============================================================================
# format_component_display
# ============================================================================

class TestFormatComponentDisplay:
    """Tests for display formatting."""

    def test_basic_display(self):
        svc = DeblendingService()
        comp = {
            'center': 15.00123,
            'center_stderr': 0.00012,
            'fwhm': 12.5,
            'fwhm_stderr': 0.8,
            'area': 1.234e-15,
            'area_stderr': 5.678e-17,
        }
        messages = svc.format_component_display(comp)

        assert len(messages) == 3
        assert 'Centroid' in messages[0]
        assert 'FWHM' in messages[1]
        assert 'Area' in messages[2]

    def test_none_stderr_handled(self):
        svc = DeblendingService()
        comp = {
            'center': 15.0,
            'center_stderr': None,
            'fwhm': 10.0,
            'fwhm_stderr': None,
            'area': 1e-15,
            'area_stderr': None,
        }
        messages = svc.format_component_display(comp)
        assert len(messages) == 3
        # Should not crash — "0" fallback for None
        assert '0' in messages[0]


# ============================================================================
# format_single_gaussian_display
# ============================================================================

class TestFormatSingleGaussianDisplay:
    """Tests for single Gaussian display formatting."""

    def test_basic_single_gaussian(self):
        svc = DeblendingService()
        params = {
            'center': 14.99,
            'center_stderr': 0.001,
            'fwhm': 11.0,
            'fwhm_stderr': 0.5,
            'area': 2e-15,
            'area_stderr': 3e-17,
        }
        messages = svc.format_single_gaussian_display(params)
        assert len(messages) == 3
        assert 'Centroid' in messages[0]

    def test_no_center_key(self):
        svc = DeblendingService()
        params = {'amplitude': 1.0}
        messages = svc.format_single_gaussian_display(params)
        assert len(messages) == 1
        assert 'Could not extract' in messages[0]
