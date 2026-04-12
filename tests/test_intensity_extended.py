# -*- coding: utf-8 -*-
"""Extended unit tests for the Intensity class — coverage gaps.

Covers:
- calc_intensity_batch_fast() with pre-computed context
- prepare_batch_context() and prepare_overlap_structure()
- build_table() — full_range and wavelength_range modes
- get_lines_in_range_with_intensity()
- get_table_in_range()
- bulk_parameter_update_vectorized()
- register_strategy() / _get_strategy()
- curve_growth_no_overlap method
- _repr_html_() notebook display
"""

import pytest
import numpy as np

import iSLAT.Constants as c


class _IntensityTestBase:
    """Shared helpers for Intensity tests."""

    def _make_line_list(self, n_lines=10, lam_start=10.0, lam_step=1.0,
                        molecule_id='TestMol'):
        from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
        rng = np.random.default_rng(42)
        lines_data = []
        for i in range(n_lines):
            lam = lam_start + i * lam_step
            freq = c.SPEED_OF_LIGHT_MICRONS / lam
            lines_data.append({
                'nr': i + 1,
                'lev_up': f'0|{2 * i + 2}',
                'lev_low': f'0|{2 * i + 1}',
                'lam': lam,
                'freq': freq,
                'a_stein': rng.uniform(1e-3, 5e-2),
                'e_up': 1000.0 + i * 500,
                'e_low': 500.0 + i * 500,
                'g_up': 2 * i + 3,
                'g_low': 2 * i + 1,
            })
        mll = MoleculeLineList(molecule_id=molecule_id, lines_data=lines_data)
        mll.partition_function = mll._partition_type(
            t=np.array([100, 200, 300, 500, 1000, 2000], dtype=np.float64),
            q=np.array([10, 50, 150, 500, 2000, 8000], dtype=np.float64),
        )
        mll._molar_mass = 18.015
        return mll

    def _make_intensity(self, mll=None, t_kin=500.0, n_mol=1e18, dv=1.0):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        if mll is None:
            mll = self._make_line_list()
        inten = Intensity(mll)
        inten.calc_intensity(t_kin=t_kin, n_mol=n_mol, dv=dv)
        return inten


class TestCalcIntensityBatchFast(_IntensityTestBase):
    """Tests for the fast batch API used in fitting."""

    def test_prepare_batch_context_keys(self):
        """prepare_batch_context should return the expected key set."""
        inten = self._make_intensity()
        ctx = inten.prepare_batch_context()
        expected_keys = {'freq', 'e_up', 'e_low', 'e_delta', 'line_scalar',
                         'freq_ratio', 'partition_t', 'partition_q',
                         'x_quad', 'w_quad', 'exp_neg_x2', 'sqrt_pi_sum',
                         'overlap', 'n_lines',
                         'e_low_f32', 'e_delta_f32', 'line_scalar_f32',
                         'freq_ratio_f32', 'freq_f32',
                         'x_quad_f32', 'w_quad_f32', 'exp_neg_x2_f32',
                         'sqrt_pi_sum_f32', 'two_sqrt_ln2_f32', 'tau_thin_f32',
                         'line_indices'}
        assert expected_keys.issubset(ctx.keys())

    def test_prepare_batch_context_n_lines(self):
        mll = self._make_line_list(n_lines=7)
        inten = self._make_intensity(mll)
        ctx = inten.prepare_batch_context()
        assert ctx['n_lines'] == 7
        assert ctx['freq'].shape == (7,)

    def test_prepare_batch_context_line_indices(self):
        """Providing line_indices should slice all arrays."""
        inten = self._make_intensity()
        idx = np.array([0, 2, 5])
        ctx = inten.prepare_batch_context(line_indices=idx)
        assert ctx['n_lines'] == 3
        assert ctx['freq'].shape == (3,)
        assert ctx['line_indices'] is None  # baked in

    def test_prepare_overlap_structure_all_isolated(self):
        """With well-separated lines, all should be isolated."""
        mll = self._make_line_list(n_lines=5, lam_step=2.0)
        inten = self._make_intensity(mll)
        ctx = inten.prepare_batch_context()
        inten.prepare_overlap_structure(ctx, dv=1.0)
        assert ctx['overlap']['all_isolated'] is True

    def test_prepare_overlap_structure_single_line(self):
        mll = self._make_line_list(n_lines=1)
        inten = self._make_intensity(mll)
        ctx = inten.prepare_batch_context()
        inten.prepare_overlap_structure(ctx, dv=1.0)
        assert ctx['overlap']['all_isolated'] is True

    def test_batch_fast_matches_regular_batch(self):
        """calc_intensity_batch_fast should give similar results to calc_intensity_batch."""
        mll = self._make_line_list(n_lines=6, lam_step=1.5)
        inten = self._make_intensity(mll)
        ctx = inten.prepare_batch_context()
        inten.prepare_overlap_structure(ctx, dv=1.0)

        t_arr = np.array([400.0, 600.0, 800.0])
        n_arr = np.array([1e17, 1e18, 1e19])
        dv_arr = np.array([1.0, 1.0, 1.0])

        result_fast = inten.calc_intensity_batch_fast(t_arr, n_arr, dv_arr, ctx)
        result_reg = inten.calc_intensity_batch(t_arr, n_arr, dv_arr)

        # float32 vs float64 means tolerance is looser
        np.testing.assert_allclose(result_fast, result_reg, rtol=1e-2, atol=1e-30)

    def test_batch_fast_shape(self):
        mll = self._make_line_list(n_lines=8)
        inten = self._make_intensity(mll)
        ctx = inten.prepare_batch_context()
        inten.prepare_overlap_structure(ctx, dv=1.0)

        t_arr = np.array([300.0, 500.0])
        n_arr = np.array([1e17, 1e18])
        dv_arr = np.array([1.0, 1.0])

        result = inten.calc_intensity_batch_fast(t_arr, n_arr, dv_arr, ctx)
        assert result.shape == (2, 8)

    def test_batch_fast_nonnegative(self):
        inten = self._make_intensity()
        ctx = inten.prepare_batch_context()
        inten.prepare_overlap_structure(ctx, dv=1.0)

        result = inten.calc_intensity_batch_fast(
            np.array([500.0]), np.array([1e18]), np.array([1.0]), ctx
        )
        assert np.all(result >= 0)


class TestBuildTable(_IntensityTestBase):
    """Tests for build_table() — DataFrame generation."""

    def test_build_table_default(self):
        """Default build_table returns active-range lines with intensity."""
        inten = self._make_intensity()
        df = inten.build_table()
        assert 'lam' in df.columns
        assert 'tau' in df.columns
        assert 'intens' in df.columns
        assert 'a_stein' in df.columns
        assert len(df) == 10

    def test_build_table_wavelength_range_filter(self):
        """wavelength_range kwarg should filter the table."""
        inten = self._make_intensity()
        df = inten.build_table(wavelength_range=(12.0, 16.0))
        # Only lines with lam in [12, 16]
        assert all(df['lam'] >= 12.0)
        assert all(df['lam'] <= 16.0)
        assert len(df) < 10

    def test_build_table_no_intensity(self):
        """If calc_intensity not called, tau/intens should be NaN."""
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_line_list()
        inten = Intensity(mll)
        df = inten.build_table()
        assert all(np.isnan(df['tau']))
        assert all(np.isnan(df['intens']))


class TestGetLinesInRange(_IntensityTestBase):
    """Tests for get_lines_in_range_with_intensity()."""

    def test_returns_tuples(self):
        inten = self._make_intensity()
        results = inten.get_lines_in_range_with_intensity(12.0, 16.0)
        assert isinstance(results, list)
        for item in results:
            assert len(item) == 3  # (line, intensity, tau)

    def test_empty_before_calc(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        mll = self._make_line_list()
        inten = Intensity(mll)
        results = inten.get_lines_in_range_with_intensity(10.0, 20.0)
        assert results == []

    def test_all_lines_in_full_range(self):
        inten = self._make_intensity()
        results = inten.get_lines_in_range_with_intensity(9.0, 21.0)
        assert len(results) == 10

    def test_filtered_by_range(self):
        inten = self._make_intensity()
        results = inten.get_lines_in_range_with_intensity(12.0, 14.5)
        lams = [r[0].lam for r in results]
        for l in lams:
            assert 12.0 <= l <= 14.5

    def test_intensity_and_tau_positive(self):
        inten = self._make_intensity()
        results = inten.get_lines_in_range_with_intensity(10.0, 20.0)
        for _, i_val, tau_val in results:
            assert i_val >= 0
            assert tau_val >= 0


class TestGetTableInRange(_IntensityTestBase):
    """Tests for get_table_in_range() convenience method."""

    def test_returns_dataframe(self):
        inten = self._make_intensity()
        df = inten.get_table_in_range(12.0, 16.0)
        assert hasattr(df, 'columns')
        assert 'lam' in df.columns

    def test_filters_correctly(self):
        inten = self._make_intensity()
        df = inten.get_table_in_range(13.0, 15.0)
        assert all(df['lam'] >= 13.0)
        assert all(df['lam'] <= 15.0)


class TestBulkParameterUpdate(_IntensityTestBase):
    """Tests for bulk_parameter_update_vectorized()."""

    def test_basic_batch(self):
        inten = self._make_intensity()
        combos = [
            {'t_kin': 400.0, 'n_mol': 1e17, 'dv': 1.0},
            {'t_kin': 600.0, 'n_mol': 1e18, 'dv': 2.0},
            {'t_kin': 800.0, 'n_mol': 1e19, 'dv': 1.5},
        ]
        result = inten.bulk_parameter_update_vectorized(combos)
        assert result.shape == (3, 10)
        assert np.all(result >= 0)

    def test_empty_combos(self):
        inten = self._make_intensity()
        result = inten.bulk_parameter_update_vectorized([])
        assert len(result) == 0

    def test_single_combo_matches_calc_intensity(self):
        inten = self._make_intensity()
        combos = [{'t_kin': 500.0, 'n_mol': 1e18, 'dv': 1.0}]
        result = inten.bulk_parameter_update_vectorized(combos)
        # Compare to regular calc_intensity
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0)
        np.testing.assert_allclose(result[0], inten.intensity, rtol=1e-10)


class TestStrategyPattern(_IntensityTestBase):
    """Tests for the strategy registry and dispatch."""

    def test_default_strategies_exist(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        Intensity._ensure_strategy_registry()
        assert 'curve_growth' in Intensity._strategy_registry
        assert 'radex' in Intensity._strategy_registry
        assert 'curve_growth_no_overlap' in Intensity._strategy_registry

    def test_get_strategy_unknown_raises(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity
        with pytest.raises(ValueError, match="Unknown intensity method"):
            Intensity._get_strategy('nonexistent_method')

    def test_register_custom_strategy(self):
        from iSLAT.Modules.DataTypes.Intensity import Intensity

        class DummyStrategy:
            def calculate(self, intensity_obj, center_tau, dv_vals,
                         bb_vals, freq_ratio, sqrt_ln2_inv, frequencies):
                return np.zeros_like(center_tau)

        Intensity.register_strategy('dummy_test', DummyStrategy())
        strategy = Intensity._get_strategy('dummy_test')
        assert isinstance(strategy, DummyStrategy)
        # Clean up
        del Intensity._strategy_registry['dummy_test']

    def test_curve_growth_no_overlap_method(self):
        """curve_growth_no_overlap should produce valid results."""
        inten = self._make_intensity()
        inten.invalidate_cache()
        inten.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0,
                            method='curve_growth_no_overlap')
        assert inten.intensity is not None
        assert np.all(inten.intensity >= 0)

    def test_curve_growth_vs_no_overlap_similar_for_isolated(self):
        """With well-separated lines, overlap and no-overlap should agree closely."""
        mll = self._make_line_list(n_lines=5, lam_step=2.0)

        from iSLAT.Modules.DataTypes.Intensity import Intensity
        inten1 = Intensity(mll)
        inten1.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0,
                             method='curve_growth')
        inten2 = Intensity(mll)
        inten2.calc_intensity(t_kin=500.0, n_mol=1e18, dv=1.0,
                             method='curve_growth_no_overlap')
        np.testing.assert_allclose(inten1.intensity, inten2.intensity, rtol=1e-2)


class TestIntensityReprHtml(_IntensityTestBase):
    """Test the _repr_html_ method if it exists."""

    def test_repr_basic(self):
        """__repr__ should contain molecule name."""
        inten = self._make_intensity()
        r = repr(inten)
        assert 'TestMol' in r
        assert 'Intensity' in r

    def test_repr_with_calculated_values(self):
        """__repr__ should show calculated parameters."""
        inten = self._make_intensity(t_kin=700.0)
        r = repr(inten)
        assert '700' in r
