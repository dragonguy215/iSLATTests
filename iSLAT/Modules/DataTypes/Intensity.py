# -*- coding: utf-8 -*-

"""
The class Intensity calculates the intensities

* The same algorithm as in the Fortran 90 code are used. This is explained in the appendix of Banzatti et al. 2012.
        frequencies : np.ndarray
            Line frequencies in Hz with shape (n_lines,)
        sqrt_ln2_inv : float
            Normalization constant

        Returns read access to the fields is granted through properties

- 01/06/2020: SB, initial version
- 11/01/2025 Johnny McCaskill, redesigned for performance and to enable overlapping line treatment

"""

import numpy as np
import time
from typing import Optional, Union, Literal, TYPE_CHECKING, Any
import iSLAT.Constants as c

# Performance logging
from iSLAT.Modules.Debug.PerformanceLogger import perf_log, log_timing, PerformanceSection

# Lazy imports for performance
_scipy_imported = False
_pd_imported = False

def _get_scipy():
    """Lazy import of scipy components"""
    global _scipy_imported
    if not _scipy_imported:
        global fixed_quad
        from scipy.integrate import fixed_quad
        _scipy_imported = True
    return fixed_quad

def _get_pandas():
    """Lazy import of pandas"""
    global _pd_imported, pd
    if not _pd_imported:
        try:
            import pandas as pd
            _pd_imported = True
        except ImportError:
            pd = None
            _pd_imported = True
    return pd

if TYPE_CHECKING:
    import pandas
    from .MoleculeLineList import MoleculeLineList
else:
    # Lazy import for runtime
    def _get_molecule_line_list():
        from .MoleculeLineList import MoleculeLineList
        return MoleculeLineList

__all__ = ["Intensity"]

class Intensity:
    __slots__ = ('_molecule', '_intensity', '_tau', '_t_kin', '_n_mol', '_dv', '_cache_valid', '_sorted_idx', '_sorted_freq',
                 '_cached_line_scalar', '_cached_e_delta', '_cached_dv_star',
                 '_cached_freq_ratio', '_wavelength_range')
    
    def __init__(self, molecule_line_list: 'MoleculeLineList',
                 wavelength_range: Optional[tuple] = None) -> None:
        """Initialize an intensity class which calculates the intensities for a given molecule and provided
        physical parameters.

        Parameters
        ----------
        molecule_line_list: MoleculeLineList
            Molecular line list data to calculate the intensity
        wavelength_range: tuple of (float, float), optional
            Wavelength range ``(lam_min, lam_max)`` in microns.  When set,
            only lines within this range are used for intensity calculations.
            The range is propagated to the underlying MoleculeLineList.
        """
        if TYPE_CHECKING:
            self._molecule: "MoleculeLineList" = molecule_line_list
        else:
            # Runtime: accept any compatible object
            self._molecule = molecule_line_list
        
        self._wavelength_range = wavelength_range
        if wavelength_range is not None:
            self._molecule.wavelength_range = wavelength_range
            
        self._intensity: Optional[np.ndarray] = None
        self._tau: Optional[np.ndarray] = None
        self._t_kin: Optional[float] = None
        self._n_mol: Optional[float] = None
        self._dv: Optional[float] = None
        self._cache_valid: bool = False
        self._sorted_idx: Optional[np.ndarray] = None
        self._sorted_freq: Optional[np.ndarray] = None
        # Pre-computed line constants (molecule-dependent, not parameter-dependent)
        self._cached_line_scalar: Optional[np.ndarray] = None
        self._cached_e_delta: Optional[np.ndarray] = None
        self._cached_dv_star: Optional[np.ndarray] = None  # per-neighbor dv thresholds for overlap
        self._cached_freq_ratio: Optional[np.ndarray] = None  # freq / c, shape (1, n_lines)

    # Pre-computed constants for blackbody calculation (class-level for efficiency)
    _BB_COEFF_RJ = 2.0 * c.BOLTZMANN_CONSTANT / (c.SPEED_OF_LIGHT_CGS ** 2)
    _BB_COEFF_PLANCK = 2.0 * c.PLANCK_CONSTANT / (c.SPEED_OF_LIGHT_CGS ** 2)
    _BB_X_FACTOR = c.PLANCK_CONSTANT / c.BOLTZMANN_CONSTANT
    
    # Pre-computed constants for intensity calculation (avoid recomputing every call)
    _SQRT_LN2_INV = 1.0 / (2.0 * np.sqrt(np.log(2.0)))  # ~ 0.6005
    _C3_OVER_8PI = (c.SPEED_OF_LIGHT_CGS ** 3) / (8.0 * np.pi)  # c^3/(8π)
    _INV_FGAUSS_1E5 = 1.0 / (1e5 * c.FGAUSS_PREFACTOR)  # 1/(10^5 * FGAUSS_PREFACTOR)
    _TWO_SQRT_LN2 = 2.0 * np.sqrt(np.log(2.0))  # ~ 1.6651, for converting dv to x-units
    
    # tau-regime shortcut threshold
    _TAU_THIN = 0.01     # below this, F(tau_0) ~ tau_0 * Sum_k w_k exp(-x_k^2)  (linear regime)
    
    # Pre-computed sqrt(pi) sum for thin-regime shortcut (set by _initialize_gauss_quad)
    _GAUSS_QUAD_SQRT_PI_SUM = None  # Sum_k w_k * exp(-x_k^2) ~ sqrt(pi)
    
    # Weak-blend early-exit threshold in Gaussian widths
    # At 4*sigma, overlap contribution is exp(-16) ~ 1.1e-7 -- negligible
    _WEAK_BLEND_THRESHOLD = 4.0
    
    @staticmethod
    def _bb(nu: np.ndarray, T: np.ndarray) -> np.ndarray:
        """Blackbody function for temperatures and frequencies.
        Handles both scalar T and array T for batch operations.
        Uses optimized approximations for accuracy and performance.

        Parameters
        ----------
        nu: np.ndarray
            Frequency in Hz to calculate the blackbody values, shape (n_lines,)
        T: np.ndarray or float
            Temperature(s) in K. Can be scalar or array shape (n_conditions,)

        Returns
        -------
        np.ndarray:
            Blackbody intensity in erg/s/cm**2/sr/Hz
            Shape (n_lines,) if T is scalar, (n_conditions, n_lines) if T is array
        """
        cls = Intensity
        T = np.atleast_1d(T)
        
        # Check if we have multiple temperatures (batch mode)
        if T.size == 1:
            # Single temperature - branchless scalar path
            # np.clip prevents overflow in np.expm1 for x > ~709 (Wien regime)
            # and for x < 1e-5 (Rayleigh-Jeans), expm1 is already accurate.
            T_val = T.item()
            inv_T = 1.0 / T_val
            x = cls._BB_X_FACTOR * inv_T * nu
            
            nu_cu = nu * nu * nu
            
            return cls._BB_COEFF_PLANCK * nu_cu / np.expm1(np.clip(x, 1e-300, 709.0))
        else:
            # Multiple temperatures - vectorized batch path
            # Shape: T is (n_cond,), nu is (n_lines,)
            # Result should be (n_cond, n_lines)
            inv_T = 1.0 / T  # (n_cond,)
            
            # Broadcast: (n_cond, 1) * (n_lines,) -> (n_cond, n_lines)
            x = cls._BB_X_FACTOR * inv_T[:, np.newaxis] * nu[np.newaxis, :]
            
            nu_sq = nu * nu  # (n_lines,)
            nu_cu = nu_sq * nu  # (n_lines,)
            
            out = np.empty(x.shape, dtype=np.float64)
            
            m_rj = x < 1.0e-5
            m_wien = x > 20.0
            m_mid = ~(m_rj | m_wien)
            
            # Compute directly into masked positions to avoid full-size temporaries
            # Use np.broadcast_to for zero-copy broadcast views
            T_2d = np.broadcast_to(T[:, np.newaxis], x.shape)
            nu_sq_2d = np.broadcast_to(nu_sq[np.newaxis, :], x.shape)
            nu_cu_2d = np.broadcast_to(nu_cu[np.newaxis, :], x.shape)
            
            # Rayleigh-Jeans regime
            if np.any(m_rj):
                out[m_rj] = cls._BB_COEFF_RJ * T_2d[m_rj] * nu_sq_2d[m_rj]
            
            # Wien regime
            if np.any(m_wien):
                out[m_wien] = cls._BB_COEFF_PLANCK * nu_cu_2d[m_wien] * np.exp(-x[m_wien])
            
            # Full Planck regime
            if np.any(m_mid):
                out[m_mid] = cls._BB_COEFF_PLANCK * nu_cu_2d[m_mid] / np.expm1(x[m_mid])
            
            return out

    # Pre-computed Gaussian quadrature points and weights for maximum efficiency
    _GAUSS_QUAD_X = None
    _GAUSS_QUAD_W = None
    _GAUSS_QUAD_INITIALIZED = False
    
    @classmethod
    def _initialize_gauss_quad(cls):
        """Initialize Gaussian quadrature points and weights once for all calculations."""
        if not cls._GAUSS_QUAD_INITIALIZED:
            try:
                from numpy.polynomial.legendre import leggauss
                x_quad, w_quad = leggauss(20)
                # Transform from [-1,1] to [-6,6] and pre-compute exp(-x^2)
                #cls._GAUSS_QUAD_X = np.exp(-(6.0 * x_quad) ** 2)  # Pre-compute exp(-x^2)
                cls._GAUSS_QUAD_X = 6.0 * x_quad  # Store actual quadrature points
                cls._GAUSS_QUAD_EXP = np.exp(-(cls._GAUSS_QUAD_X ** 2))
                cls._GAUSS_QUAD_W = 6.0 * w_quad
                cls._GAUSS_QUAD_SQRT_PI_SUM = float(np.dot(cls._GAUSS_QUAD_W, cls._GAUSS_QUAD_EXP))
                cls._GAUSS_QUAD_INITIALIZED = True
            except ImportError:
                # Fallback to scipy method if numpy.polynomial not available
                fixed_quad = _get_scipy()
                # Use a simple tau value to initialize quadrature
                cls._GAUSS_QUAD_X = np.linspace(-6, 6, 20)  # Store actual x values
                cls._GAUSS_QUAD_EXP = np.exp(-(cls._GAUSS_QUAD_X ** 2))
                cls._GAUSS_QUAD_W = np.ones(20) * (12.0 / 20.0)  # Simple uniform weighting
                cls._GAUSS_QUAD_SQRT_PI_SUM = float(np.dot(cls._GAUSS_QUAD_W, cls._GAUSS_QUAD_EXP))
                cls._GAUSS_QUAD_INITIALIZED = True

    def _fint_multi(self, center_tau: np.ndarray, dv_cond: np.ndarray,
                   bb_vals: np.ndarray, freq_ratio: np.ndarray, sqrt_ln2_inv: float,
                   frequencies: np.ndarray) -> np.ndarray:
        """
        Calculate spectral line intensities with proper treatment of overlapping lines.
        
        This method implements the curve-of-growth integral for molecular line intensities,
        handling both isolated and overlapping lines using correct radiative transfer physics.
        For overlapping lines, optical depths are summed before applying the curve-of-growth
        to avoid unphysical intensity enhancement and maintain conservation.
        
        Algorithm Overview:
        1. Identify isolated vs. overlapping lines using velocity separation criteria  
        2. Process isolated lines in batch using vectorized operations
        3. For overlapping lines: sum optical depths -> apply curve-of-growth -> distribute intensity
        
        Physical Basis:
        - Isolated lines: I = (physical_factors) * integral[1 - exp(-tau)] dv
        - Overlapping lines: I_total = integral[1 - exp(-Sum(tau_i))] dv, then distribute proportionally
        
        Performance:
        - Memory complexity: O(n) instead of O(n^2) for overlap detection
        - Time complexity: O(n log n) for sorting + O(n) for processing  
        - Vectorized operations for isolated lines (typically 85%+ of all lines)
        
        Parameters
        ----------
        center_tau : np.ndarray, shape (n_conditions, n_lines)
            Line center optical depths for all physical conditions and lines
        dv_cond : np.ndarray, shape (n_conditions,)
            Doppler broadening velocities (FWHM) in km/s for each condition
        bb_vals : np.ndarray, shape (n_conditions, n_lines)
            Blackbody source function values at line frequencies  
        freq_ratio : np.ndarray, shape (n_conditions, n_lines)
            Frequency ratios (v/c) for unit conversion
        sqrt_ln2_inv : float
            Normalization constant: 1/(2*sqrt(ln(2))) ~ 0.6005 for Gaussian integration
        frequencies : np.ndarray, shape (n_lines,)
            Line rest frequencies in Hz
            
        Returns
        -------
        np.ndarray, shape (n_conditions, n_lines)
            Line-integrated specific intensities following Banzatti et al. 2012 Eq. A1.
            Units: (1/(2*sqrt(ln(2)))) x 10^5 x dv x (nu/c) x B_nu(T) x F(tau_0), where F is the
            curve-of-growth integral. The 10^5 factor converts km/s to cm/s.
            
        Notes
        -----
        - Uses Gaussian quadrature for numerical integration of curve-of-growth
        - Velocity overlap criterion: lines are grouped when their frequency
          separation corresponds to less than the current dv (FWHM)
        """
        cls = self.__class__
        # Initialize quadrature points and weights
        if not cls._GAUSS_QUAD_INITIALIZED:
            cls._initialize_gauss_quad()
        
        # Setup arrays and extract dimensions
        dv_cond = np.atleast_1d(dv_cond)
        n_cond, n_lines = center_tau.shape
        
        # Quadrature setup for curve-of-growth integration
        exp_neg_x_squared = cls._GAUSS_QUAD_EXP  # (n_quad,) - don't add dimension yet
        weights = cls._GAUSS_QUAD_W
        
        # ===============================================================
        # STEP 1: Identify isolated vs. overlapping lines
        # ===============================================================
        
        # Pre-compute physical factors once (used by all code paths below)
        physical_factors = np.empty((n_cond, n_lines), dtype=np.float64)
        np.multiply(freq_ratio, bb_vals, out=physical_factors)
        physical_factors *= (sqrt_ln2_inv * 1e5)
        physical_factors *= dv_cond[:, np.newaxis]

        # Cache sorted order and per-neighbor dv_star thresholds per instance
        # (frequencies are molecule-dependent and stable across calls).
        # dv_star[i] = the minimum dv (km/s) at which sorted lines i and i+1 overlap:
        #   dv_star_i = c * dnu_i / max(nu_i, nu_{i+1})
        if self._sorted_idx is None:
            self._sorted_idx = np.argsort(frequencies)
            self._sorted_freq = frequencies[self._sorted_idx]
            sf = self._sorted_freq
            if sf.size > 1:
                self._cached_dv_star = (c.SPEED_OF_LIGHT_KMS
                    * np.diff(sf) / np.maximum(sf[:-1], sf[1:]))
            else:
                self._cached_dv_star = np.array([], dtype=np.float64)

        sort_indices = self._sorted_idx
        dv_star = self._cached_dv_star  # (n_sorted - 1,)
        
        # Use the effective dv for grouping. In batch mode with varying dv,
        # use max(dv_cond) to find the *superset* of potentially blended pairs.
        # The exact multi-line quadrature in Step 3 handles per-condition
        # offsets correctly, so over-grouping here is conservative, not wrong.
        # For single-condition calls, this is exact.
        dv_group = np.max(dv_cond)
        
        # Fast all-isolated check: if even the widest dv doesn't blend any pair
        if dv_star.size == 0 or np.all(dv_group < dv_star):
            # tau-regime shortcuts: partition lines into thin / normal / thick
            fint_vals = self._cog_with_tau_shortcuts(center_tau, exp_neg_x_squared, weights)
            
            # Final multiplication in-place
            np.multiply(physical_factors, fint_vals, out=physical_factors)
            return physical_factors
        
        # Vectorized overlap grouping using cumsum-based labeling.
        # This avoids np.split (O(n_groups) Python objects) and instead assigns
        # a contiguous integer group label to each sorted line, then uses
        # np.bincount to identify isolated vs blended groups in O(n) time.
        is_boundary = dv_group < dv_star
        
        # Assign group labels via cumulative sum of boundary flags
        group_labels = np.empty(n_lines, dtype=np.int32)
        group_labels[0] = 0
        group_labels[1:] = np.cumsum(is_boundary)
        
        # Count lines per group using bincount
        group_counts = np.bincount(group_labels)
        
        # Map group counts back to each sorted line position
        line_group_counts = group_counts[group_labels]
        
        # Isolated lines: groups with exactly 1 member
        isolated_sorted_mask = line_group_counts == 1
        blended_sorted_mask = ~isolated_sorted_mask
        
        # ===============================================================
        # STEP 2: Process isolated lines
        # ===============================================================
        
        if np.any(isolated_sorted_mask):
            isolated_indices = sort_indices[isolated_sorted_mask]
            tau_batch = center_tau[:, isolated_indices]
            
            # tau-regime shortcuts: partition lines into thin / normal / thick
            fint_vals = self._cog_with_tau_shortcuts(tau_batch, exp_neg_x_squared, weights)
            
            # Multiply curve-of-growth factors in-place into physical_factors
            physical_factors[:, isolated_indices] *= fint_vals
        
        # ===============================================================
        # STEP 3: Process overlapping lines with separation-aware physics
        # ===============================================================
        # Exact multi-line quadrature for blended groups:
        #   Computes tau_tot(x_k) = Sum_j tau_{0,j} exp(-(x_k - delta_j)^2) at each
        #   quadrature point, applies [1 - exp(-tau_tot)] correctly, then
        #   distributes intensity to lines by their fractional opacity.
        #
        # This fixes the co-centered approximation that ignores line separations,
        # which over-saturates peak optical depth and underestimates emission
        # for blends separated by ~0.5-2 line widths.
        
        if not np.any(blended_sorted_mask):
            return physical_factors
        
        x_quad = cls._GAUSS_QUAD_X  # (n_quad,) actual quadrature points in [-6, 6]
        two_sqrt_ln2 = cls._TWO_SQRT_LN2
        
        # Extract blended line data in sorted order
        blended_sort_positions = np.where(blended_sorted_mask)[0]
        blended_orig_indices = sort_indices[blended_sort_positions]
        blended_group_labels = group_labels[blended_sort_positions]
        
        # Find group boundaries within the blended-only arrays
        unique_blend_labels, blend_group_starts = np.unique(blended_group_labels, return_index=True)
        blend_group_ends = np.append(blend_group_starts[1:], len(blended_orig_indices))
        blend_group_sizes = blend_group_ends - blend_group_starts
        
        # ---------------------------------------------------------------
        # STEP 3a: Batch-vectorized processing for common group sizes
        # ---------------------------------------------------------------
        # Groups of the same size can be processed in a single vectorized
        # operation, avoiding per-group Python overhead.
        
        _BATCH_MAX_SIZE = 8  # Batch groups up to this size; larger ones use per-group loop
        
        processed_mask = np.zeros(len(unique_blend_labels), dtype=bool)
        
        for target_size in range(2, _BATCH_MAX_SIZE + 1):
            size_mask = blend_group_sizes == target_size
            if not np.any(size_mask):
                continue
            
            size_group_idx = np.where(size_mask)[0]
            n_groups_this_size = len(size_group_idx)
            
            # Gather indices for all groups of this size into (n_groups, target_size)
            all_indices = np.empty((n_groups_this_size, target_size), dtype=np.intp)
            for k in range(target_size):
                all_indices[:, k] = blended_orig_indices[blend_group_starts[size_group_idx] + k]
            
            # Process each condition
            for ic in range(n_cond):
                tau_batch_g = center_tau[ic][all_indices]     # (n_groups, target_size)
                freq_batch_g = frequencies[all_indices]       # (n_groups, target_size)
                
                # Skip negligible groups
                tau_max_per_group = tau_batch_g.max(axis=1)
                active = tau_max_per_group > 1e-300
                if not np.any(active):
                    continue
                
                tau_bg = tau_batch_g[active]                   # (n_active, target_size)
                freq_bg = freq_batch_g[active]                 # (n_active, target_size)
                active_indices = all_indices[active]           # (n_active, target_size)
                
                # tau-weighted reference frequency
                tau_sums = tau_bg.sum(axis=1)                  # (n_active,)
                safe_sums = np.where(tau_sums > 1e-300, tau_sums, 1.0)
                freq_refs = ((tau_bg * freq_bg).sum(axis=1) / safe_sums)[:, np.newaxis]
                
                # Velocity offsets in quadrature x-units
                dv_kms = (freq_bg - freq_refs) * (c.SPEED_OF_LIGHT_KMS / freq_refs)
                delta_x = (dv_kms / dv_cond[ic]) * two_sqrt_ln2  # (n_active, target_size)
                
                # Shifted Gaussian profiles in-place: (n_quad, n_active, target_size)
                gp = x_quad[:, np.newaxis, np.newaxis] - delta_x[np.newaxis, :, :]
                np.multiply(gp, gp, out=gp)            # gp = shifted_x^2
                np.negative(gp, out=gp)                 # gp = -shifted_x^2
                np.exp(gp, out=gp)                      # gp = exp(-shifted_x^2)
                
                # tau_tot at each quadrature point: (n_quad, n_active)
                tau_at_quad = np.einsum('qgl,gl->qg', gp, tau_bg)
                
                # Curve-of-growth: 1 - exp(-tau_tot)
                # Keep tau_at_quad intact for the division; compute integrand into a separate array
                cog_mask = tau_at_quad > 1e-300
                integrand = np.empty_like(tau_at_quad)
                np.negative(tau_at_quad, out=integrand)         # -tau
                np.expm1(integrand, out=integrand)              # expm1(-tau)
                np.negative(integrand, out=integrand)           # 1-exp(-tau)
                cog_over_tau = np.where(cog_mask, integrand / np.where(cog_mask, tau_at_quad, 1.0), 0.0)
                
                # Distribute intensity proportionally, reusing gp array in-place
                gp *= (cog_over_tau * weights[:, np.newaxis])[:, :, np.newaxis]
                fint = gp.sum(axis=0)                          # (n_active, target_size)
                fint *= tau_bg
                
                # Scatter results back into physical_factors
                physical_factors[ic][active_indices.ravel()] *= fint.ravel()
            
            processed_mask[size_group_idx] = True
        
        # ---------------------------------------------------------------
        # STEP 3b: Per-group loop for remaining large groups
        # ---------------------------------------------------------------
        remaining_idx = np.where(~processed_mask)[0]
        
        if len(remaining_idx) > 0:
            weak_blend_thresh = cls._WEAK_BLEND_THRESHOLD
            tau_thin_thresh = cls._TAU_THIN
            
            for gi in remaining_idx:
                start = blend_group_starts[gi]
                end = blend_group_ends[gi]
                overlapping_indices = blended_orig_indices[start:end]
                
                # Extract optical depths for this blended group
                tau_group = center_tau[:, overlapping_indices]  # (n_cond, n_overlap)
                
                # Early exit: skip quadrature if all tau_0 values are negligible
                if tau_group.max() < 1e-300:
                    continue
                
                freq_group = frequencies[overlapping_indices]
                freq_group_mean = freq_group.mean()
                
                for ic in range(n_cond):
                    tg = tau_group[ic]
                    
                    if tg.max() < 1e-300:
                        continue
                    
                    # tau_0-weighted reference frequency
                    tau_sum_ic = tg.sum()
                    if tau_sum_ic > 1e-300:
                        freq_ref_ic = np.dot(tg, freq_group) / tau_sum_ic
                    else:
                        freq_ref_ic = freq_group_mean
                    
                    dv_kms_ic = (freq_group - freq_ref_ic) * (c.SPEED_OF_LIGHT_KMS / freq_ref_ic)
                    delta_x_ic = (dv_kms_ic / dv_cond[ic]) * two_sqrt_ln2
                    
                    # Weak-blend early exit
                    if np.max(np.abs(delta_x_ic)) > weak_blend_thresh:
                        sorted_dx = np.sort(delta_x_ic)
                        max_neighbor_sep = np.max(np.diff(sorted_dx)) if len(sorted_dx) > 1 else 0.0
                        if max_neighbor_sep > 2.0 * weak_blend_thresh:
                            tau_exp_ic = tg[:, np.newaxis] * exp_neg_x_squared
                            np.negative(tau_exp_ic, out=tau_exp_ic)
                            np.expm1(tau_exp_ic, out=tau_exp_ic)
                            np.negative(tau_exp_ic, out=tau_exp_ic)
                            fint_ic = tau_exp_ic @ weights
                            physical_factors[ic, overlapping_indices] *= fint_ic
                            continue
                    
                    # Full overlap quadrature (in-place shifted_x -> gp)
                    gp = x_quad[:, np.newaxis] - delta_x_ic[np.newaxis, :]
                    np.multiply(gp, gp, out=gp)
                    np.negative(gp, out=gp)
                    np.exp(gp, out=gp)
                    tau_at_quad_ic = gp @ tg
                    
                    if tau_at_quad_ic.max() < tau_thin_thresh:
                        fint_ic = (gp.T @ weights) * tg
                        physical_factors[ic, overlapping_indices] *= fint_ic
                        continue
                    
                    # Curve-of-growth in-place
                    cog_mask_ic = tau_at_quad_ic > 1e-300
                    integrand_ic = np.negative(tau_at_quad_ic)
                    np.expm1(integrand_ic, out=integrand_ic)
                    np.negative(integrand_ic, out=integrand_ic)  # 1 - exp(-tau)
                    cog_over_tau_ic = np.where(cog_mask_ic,
                                              integrand_ic / np.where(cog_mask_ic, tau_at_quad_ic, 1.0), 0.0)
                    # Reuse gp array for weighted sum
                    gp *= (cog_over_tau_ic * weights)[:, np.newaxis]
                    fint_ic = gp.sum(axis=0)
                    fint_ic *= tg
                    physical_factors[ic, overlapping_indices] *= fint_ic

        return physical_factors

    def _cog_with_tau_shortcuts(self, center_tau: np.ndarray,
                                exp_neg_x_squared: np.ndarray,
                                weights: np.ndarray) -> np.ndarray:
        """Curve-of-growth integration with thin-regime shortcut for isolated lines.
        
        Partitions lines into two regimes:
        - Thin (tau_0 < 0.01): F(tau_0) ~ tau_0 * Sum_k w_k exp(-x_k^2)  (linear, no exp needed)
        - Normal: full Gauss-Hermite quadrature (exact for all tau_0)
        
        Parameters
        ----------
        center_tau : np.ndarray, shape (n_cond, n_lines)
            Line center optical depths
        exp_neg_x_squared : np.ndarray, shape (n_quad,)
            Pre-computed exp(-x_k^2) at quadrature points
        weights : np.ndarray, shape (n_quad,)
            Quadrature weights
            
        Returns
        -------
        np.ndarray, shape (n_cond, n_lines)
            Curve-of-growth integral values F(tau_0)
        """
        cls = self.__class__
        n_cond, n_lines = center_tau.shape
        sqrt_pi_sum = cls._GAUSS_QUAD_SQRT_PI_SUM  # ~ sqrt(pi) from quadrature initialization
        tau_thin = cls._TAU_THIN
        
        # Find max tau across all conditions for each line to classify regime
        tau_max = np.max(center_tau, axis=0) if n_cond > 1 else center_tau[0]
        
        # Fast path: all lines optically thin
        all_thin = tau_max.max() < tau_thin if tau_max.size > 0 else False
        
        if all_thin:
            # All lines optically thin: F(tau_0) ~= tau_0 * sqrt_pi_sum
            return center_tau * sqrt_pi_sum
        
        # Fast path: no lines are thin -- full quadrature for all
        all_normal = tau_max.min() >= tau_thin if tau_max.size > 0 else True
        
        if all_normal:
            tau_exp = center_tau[:, :, np.newaxis] * exp_neg_x_squared
            np.negative(tau_exp, out=tau_exp)
            np.expm1(tau_exp, out=tau_exp)
            np.negative(tau_exp, out=tau_exp)
            return tau_exp @ weights
        
        # Mixed: some thin, some normal
        fint_vals = np.empty_like(center_tau)
        
        thin_mask = tau_max < tau_thin
        normal_mask = ~thin_mask
        
        # Thin regime: F ~ tau_0 * sqrt_pi_sum
        if np.any(thin_mask):
            fint_vals[:, thin_mask] = center_tau[:, thin_mask] * sqrt_pi_sum
        
        # Normal regime: full quadrature (handles all tau_0 values exactly)
        if np.any(normal_mask):
            tau_normal = center_tau[:, normal_mask]
            tau_exp = tau_normal[:, :, np.newaxis] * exp_neg_x_squared
            np.negative(tau_exp, out=tau_exp)
            np.expm1(tau_exp, out=tau_exp)
            np.negative(tau_exp, out=tau_exp)
            fint_vals[:, normal_mask] = tau_exp @ weights
        
        return fint_vals

    def _calc_intensity_core(self, t_kin_vals: np.ndarray, n_mol_vals: np.ndarray, 
                            dv_vals: np.ndarray, method: str = "curve_growth") -> tuple:
        """Core intensity calculation using unified vectorized operations.
        
        This method handles all intensity calculations efficiently using a single
        vectorized code path by converting scalar inputs to arrays.
        
        Parameters
        ----------
        t_kin_vals : np.ndarray
            Temperature values (K) - can be scalar or array
        n_mol_vals : np.ndarray  
            Column density values (cm^-2) - can be scalar or array
        dv_vals : np.ndarray
            Line width values (km/s) - can be scalar or array
        method : str
            Calculation method ("curve_growth" or "radex")
            
        Returns
        -------
        tuple
            (intensities, tau_values) as numpy arrays
        """
        m = self._molecule
        lines = m.lines_as_namedtuple
        partition = m.partition
        
        # Ensure inputs are numpy arrays with consistent shape
        t_kin_vals = np.asarray(t_kin_vals)
        n_mol_vals = np.asarray(n_mol_vals) 
        dv_vals = np.asarray(dv_vals)
        
        # Convert scalar inputs to 1-element arrays
        was_scalar = t_kin_vals.ndim == 0 and n_mol_vals.ndim == 0 and dv_vals.ndim == 0
        if was_scalar:
            t_kin_vals = np.atleast_1d(t_kin_vals)
            n_mol_vals = np.atleast_1d(n_mol_vals)
            dv_vals = np.atleast_1d(dv_vals)
        
        # Validate temperature bounds
        t_min, t_max = np.min(partition.t), np.max(partition.t)
        if np.any(t_kin_vals < t_min) or np.any(t_kin_vals > t_max):
            raise ValueError(f'Temperature values outside partition function range [{t_min}, {t_max}]')
        
        # Broadcast to common shape
        t_kin_vals, n_mol_vals, dv_vals = np.broadcast_arrays(t_kin_vals, n_mol_vals, dv_vals)
        output_shape = t_kin_vals.shape
        
        # Flatten for efficient processing
        t_kin_flat = t_kin_vals.ravel()
        n_mol_flat = n_mol_vals.ravel()
        dv_flat = dv_vals.ravel()
        
        # Vectorized partition function
        q_sum_vals: np.ndarray = np.interp(t_kin_flat, partition.t, partition.q)
        
        invT = (1.0 / t_kin_flat).astype(np.float64, copy=False)

        # Use pre-computed class constants for speed
        cls = self.__class__
        
        # 1D condition terms - use pre-computed constant
        cond_term = cls._INV_FGAUSS_1E5 / dv_flat                         # (n_cond,)
        q_sum_inv = 1.0 / q_sum_vals                                       # (n_cond,)
        
        # Cache line-level constants (molecule-dependent only, computed once per Intensity instance)
        if self._cached_line_scalar is None:
            freq_cubed = lines.freq ** 3
            line_term = cls._C3_OVER_8PI / freq_cubed                     # (n_lines,)
            self._cached_line_scalar = line_term * lines.g_up * lines.a_stein
            self._cached_e_delta = lines.e_up - lines.e_low
            self._cached_freq_ratio = (lines.freq / c.SPEED_OF_LIGHT_CGS)[np.newaxis, :]
        
        line_scalar = self._cached_line_scalar
        e_delta = self._cached_e_delta

        # --- Fewer exponentials for Boltzmann difference ---
        # Replace einsum with direct outer product broadcasting (faster for 2D)
        n_cond = len(invT)
        
        # E_low = invT[:, np.newaxis] * lines.e_low[np.newaxis, :]
        # Use np.outer for single condition, broadcasting for multiple
        if n_cond == 1:
            # Single condition: use direct multiplication (avoids 2D overhead)
            E_low = invT[0] * lines.e_low
            E_delta_scaled = invT[0] * e_delta
            exp_low = np.exp(-E_low)
            boltz_diff = exp_low * (-np.expm1(-E_delta_scaled))
            cond_scalar = n_mol_flat[0] * cond_term[0] * q_sum_inv[0]
            center_tau = (cond_scalar * boltz_diff * line_scalar)[np.newaxis, :]
        else:
            # Multiple conditions: use broadcasting
            E_low = np.multiply.outer(invT, lines.e_low)                  # (n_cond, n_lines)
            E_delta_scaled = np.multiply.outer(invT, e_delta)             # (n_cond, n_lines)
            exp_low = np.exp(-E_low)
            boltz_diff = exp_low * (-np.expm1(-E_delta_scaled))
            
            # Compute center_tau: cond_scalar[:, None] * boltz_diff * line_scalar[None, :]
            cond_scalar = n_mol_flat * cond_term * q_sum_inv              # (n_cond,)
            center_tau = (cond_scalar[:, np.newaxis] * boltz_diff) * line_scalar
        
        # Blackbody calculation - vectorized for all temperatures at once
        bb_vals = self._bb(lines.freq, t_kin_flat)
        
        # Use cached frequency ratio (computed once per Intensity instance)
        freq_ratio = self._cached_freq_ratio
        
        # Use pre-computed constant
        sqrt_ln2_inv = cls._SQRT_LN2_INV
        
        if method == "radex":
            intensity = (c.FGAUSS_PREFACTOR * 1e5 * dv_flat[:, np.newaxis] * 
                                freq_ratio * bb_vals * (-np.expm1(-center_tau)))
        elif method == "curve_growth":
            intensity = self._fint_multi(center_tau, dv_flat, 
                                      bb_vals, freq_ratio, sqrt_ln2_inv, lines.freq)
        elif method == "curve_growth_no_overlap":
            # Replicate the old curve_growth method without overlap treatment
            # Initialize quadrature if needed
            cls = self.__class__
            if not cls._GAUSS_QUAD_INITIALIZED:
                cls._initialize_gauss_quad()
            
            # Calculate fint using the old method (pre-computed exp(-x^2) values)
            original_shape = center_tau.shape
            tau_flat = center_tau.ravel()
            integrand = 1.0 - np.exp(-tau_flat[:, np.newaxis] * cls._GAUSS_QUAD_EXP[np.newaxis, :])
            fint_vals = np.dot(integrand, cls._GAUSS_QUAD_W).reshape(original_shape)
            
            # Calculate intensity using the old formula
            intensity = sqrt_ln2_inv * 1e5 * dv_flat[:, np.newaxis] * freq_ratio * bb_vals * fint_vals
        
        if not was_scalar:
            intensity = intensity.reshape(output_shape + (len(lines.freq),))
            center_tau = center_tau.reshape(output_shape + (len(lines.freq),))
        elif was_scalar:
            intensity = intensity.squeeze()
            center_tau = center_tau.squeeze()

        return intensity, center_tau

    def calc_intensity(self, t_kin: Optional[float] = None, n_mol: Optional[float] = None, 
                      dv: Optional[float] = None, method: Literal["curve_growth", "radex", "curve_growth_no_overlap"] = "curve_growth") -> None:
        """Calculate the intensity for a given set of physical parameters. This implements Eq. A1 and A2 in
        Banzatti et al. 2012.

        The calculation method to obtain the intensity from the opacity can be switched between the curve-of-growth
        method used in Banzatti et al. 2012, which considers broadening for high values of tau and the simple expression
        used e.g. in RADEX (van der Tak et al. 2007) which saturates at tau ~ few. For low values (tau < 1), they
        yield the same values.

        Parameters
        ----------
        t_kin: float, optional
            Kinetic temperature in K
        n_mol: float, optional
            Column density in cm**-2
        dv: float, optional
            Intrinsic (turbulent) line width in km/s
        method: Literal["curve_growth", "radex", "curve_growth_no_overlap"], default "curve_growth"
            Calculation method, either "curve_growth" for Eq. A1 or "radex" for less accurate approximation
        """
        start_time = time.perf_counter()
        
        # Check if we can use cached result
        if (self._cache_valid and 
            self._t_kin == t_kin and 
            self._n_mol == n_mol and 
            self._dv == dv and
            self._intensity is not None):
            log_timing(f"Intensity.calc_intensity(cache_hit)", time.perf_counter() - start_time, verbose=False)
            return

        # Validate inputs
        if t_kin is None or n_mol is None or dv is None:
            raise ValueError("t_kin, n_mol, and dv must all be provided")

        # Use optimized core calculation
        intensity, tau = self._calc_intensity_core(
            np.asarray(t_kin), np.asarray(n_mol), np.asarray(dv), method
        )
        
        # Update instance variables and cache
        self._t_kin = t_kin
        self._n_mol = n_mol
        self._dv = dv
        self._tau = tau
        self._intensity = intensity
        self._cache_valid = True
        
        log_timing(f"Intensity.calc_intensity({method})", time.perf_counter() - start_time)

    def calc_intensity_batch(self, t_kin_array: np.ndarray, n_mol_array: np.ndarray, 
                           dv_array: np.ndarray, method: Literal["curve_growth", "radex", "curve_growth_no_overlap"] = "curve_growth") -> np.ndarray:
        """Calculate intensities for multiple parameter combinations using vectorized operations.
        
        This method is optimized for processing many parameter sets simultaneously using the
        unified core calculation engine for maximum efficiency.

        Parameters
        ----------
        t_kin_array: np.ndarray
            Array of kinetic temperatures in K
        n_mol_array: np.ndarray  
            Array of column densities in cm**-2
        dv_array: np.ndarray
            Array of intrinsic line widths in km/s
        method: Literal["curve_growth", "radex", "curve_growth_no_overlap"], default "curve_growth"
            Calculation method

        Returns
        -------
        np.ndarray
            Array with intensities - shape depends on input array broadcasting
        """
        # Use the optimized core calculation
        intensity, _ = self._calc_intensity_core(
            np.asarray(t_kin_array), 
            np.asarray(n_mol_array), 
            np.asarray(dv_array), 
            method
        )
        
        return intensity

    def invalidate_cache(self) -> None:
        """Invalidate the calculation cache, forcing recalculation on next call."""
        self._cache_valid = False
        self._intensity = None
        self._tau = None

    def bulk_parameter_update_vectorized(self, parameter_combinations: list, method: Literal["curve_growth", "radex", "curve_growth_no_overlap"] = "curve_growth") -> np.ndarray:
        """Update multiple parameter combinations and calculate intensities in a vectorized manner.
        
        Parameters
        ----------
        parameter_combinations: list
            List of dictionaries, each containing 't_kin', 'n_mol', and 'dv' keys
        method: Literal["curve_growth", "radex", "curve_growth_no_overlap"], default "curve_growth"
            Calculation method
            
        Returns
        -------
        np.ndarray
            2D array with shape (n_combinations, n_lines) containing intensities
        """
        if not parameter_combinations:
            return np.array([])
        
        # Extract parameter arrays efficiently
        t_kin_vals = np.array([combo['t_kin'] for combo in parameter_combinations])
        n_mol_vals = np.array([combo['n_mol'] for combo in parameter_combinations]) 
        dv_vals = np.array([combo['dv'] for combo in parameter_combinations])
        
        # Use optimized core calculation only
        intensity, _ = self._calc_intensity_core(t_kin_vals, n_mol_vals, dv_vals, method)
        return intensity

    def get_table_in_range(self, lam_min: float, lam_max: float) -> "pandas.DataFrame":
        """Get a table of lines within a specific wavelength window.

        Convenience alias for ``build_table(wavelength_range=(lam_min, lam_max))``.

        Parameters
        ----------
        lam_min : float
            Minimum wavelength in microns.
        lam_max : float
            Maximum wavelength in microns.

        Returns
        -------
        pandas.DataFrame
            DataFrame with line data and computed intensity/tau for
            the requested sub-range.
        """
        return self.build_table(wavelength_range=(lam_min, lam_max))

    def get_lines_in_range_with_intensity(self, lam_min: float, lam_max: float):
        """
        Get MoleculeLine objects in the specified wavelength range with computed intensity and tau values.
        
        Parameters
        ----------
        lam_min : float
            Minimum wavelength in microns
        lam_max : float
            Maximum wavelength in microns
            
        Returns
        -------
        list[tuple[MoleculeLine, float, float]]
            List of ``(MoleculeLine, intensity, tau)`` tuples within the range.
        """
        # If no intensity calculated yet, return empty list
        if self._intensity is None or self._tau is None:
            return []
        
        # Get line data arrays once
        lines_array = self.molecule.lines_as_namedtuple
        
        # Use vectorized mask to find indices in range - O(n) instead of O(n*m)
        lam_arr = lines_array.lam
        mask = (lam_arr >= lam_min) & (lam_arr <= lam_max)
        indices_in_range = np.nonzero(mask)[0]
        
        if len(indices_in_range) == 0:
            return []
        
        # Get lines in range from the molecule line list
        lines_in_range = self.molecule.get_lines_in_range(lam_min, lam_max)
        
        # Build result list - lines_in_range and indices_in_range should correspond
        # since both use the same wavelength filtering
        result = []
        intensity_arr = self._intensity
        tau_arr = self._tau
        
        for line, idx in zip(lines_in_range, indices_in_range):
            result.append((line, intensity_arr[idx], tau_arr[idx]))
        
        return result

    @property
    def tau(self) -> Optional[np.ndarray]:
        """np.ndarray: Opacities per line"""
        return self._tau

    @property
    def intensity(self) -> Optional[np.ndarray]:
        """np.ndarray: Calculated intensity per line in erg/s/cm**2/sr/Hz"""
        return self._intensity

    @property
    def molecule(self) -> 'MoleculeLineList':
        """MoleculeLineList: Molecular line list data used for calculation"""
        return self._molecule

    @property
    def t_kin(self) -> Optional[float]:
        """float: Kinetic temperature in K used for calculation"""
        return self._t_kin

    @property
    def n_mol(self) -> Optional[float]:
        """float: Molecular column density in cm**-2 used for calculation"""
        return self._n_mol

    @property
    def dv(self) -> Optional[float]:
        """float: Line width in km/s used for calculation"""
        return self._dv

    @property
    def wavelength_range(self) -> Optional[tuple]:
        """tuple or None: Active wavelength range ``(lam_min, lam_max)`` in microns."""
        return self._wavelength_range

    @wavelength_range.setter
    def wavelength_range(self, value: Optional[tuple]):
        """Set the wavelength range filter.  Propagates to the underlying
        MoleculeLineList and invalidates all caches."""
        if value != self._wavelength_range:
            self._wavelength_range = value
            self._molecule.wavelength_range = value
            self.invalidate_cache()
            # Also invalidate pre-computed line constants since the set of lines changed
            self._sorted_idx = None
            self._sorted_freq = None
            self._cached_line_scalar = None
            self._cached_e_delta = None
            self._cached_dv_star = None
            self._cached_freq_ratio = None

    def __repr__(self) -> str:
        return f"Intensity(Mol-Name={self.molecule.name}, t_kin={self.t_kin} n_mol={self.n_mol} dv={self.dv}, " \
               f"tau={self.tau}, intensity={self.intensity})"

    def build_table(
        self,
        *,
        full_range: bool = False,
        wavelength_range: Optional[tuple[float, float]] = None,
    ) -> "pandas.DataFrame":
        """Build a DataFrame of line data with computed intensity and tau.

        By default returns only lines within the active wavelength range
        (same behaviour as the legacy ``get_table`` property).  Use the
        keyword arguments to widen or narrow the output.

        Parameters
        ----------
        full_range : bool, optional
            If ``True``, include all lines in the underlying HITRAN
            file regardless of any wavelength filter on the
            ``MoleculeLineList``.  Intensity and tau are computed for
            every line using the current physical parameters.
            Defaults to ``False`` (active wavelength range only).
        wavelength_range : tuple[float, float], optional
            ``(lam_min, lam_max)`` in microns.  When provided, the
            returned table is filtered to this sub-window.
            * If *full_range* is also ``True``, the sub-window is
              applied to the full line list.
            * If *full_range* is ``False``, the sub-window is applied
              to the active (already-filtered) line list.

        Returns
        -------
        pandas.DataFrame
            Columns: ``lev_up``, ``lev_low``, ``lam``, ``tau``,
            ``intens``, ``a_stein``, ``e_up``, ``e_low``, ``g_up``,
            ``g_low``.

        Raises
        ------
        ImportError
            If *pandas* is not installed.

        Examples
        --------
        >>> # Active range only (default)
        >>> df = intensity.build_table()

        >>> # All lines in the HITRAN file with intensity computed for each
        >>> df_full = intensity.build_table(full_range=True)

        >>> # Narrow sub-window of the active range
        >>> df_sub = intensity.build_table(wavelength_range=(4.9, 5.0))
        """
        pd = _get_pandas()
        if pd is None:
            raise ImportError("Pandas required to create table")

        df = self._build_table(pd, full_range=full_range)

        # Apply optional wavelength sub-window
        if wavelength_range is not None:
            lam_min, lam_max = wavelength_range
            mask = (df["lam"] >= lam_min) & (df["lam"] <= lam_max)
            df = df.loc[mask].reset_index(drop=True)

        return df

    # ---- internal table builder -------------------------------------

    # Column order used by all table-building paths
    _TABLE_COLUMNS = ('lev_up', 'lev_low', 'lam', 'tau', 'intens',
                      'a_stein', 'e_up', 'e_low', 'g_up', 'g_low')

    def _build_table(self, pd, *, full_range: bool) -> "pandas.DataFrame":
        """Unified table builder for both active-range and full-range modes.

        Parameters
        ----------
        pd : module
            The pandas module (already imported by caller).
        full_range : bool
            If ``True``, use all lines from the HITRAN file and
            compute intensity and tau for every line, not just those
            in the active wavelength window.  This reproduces the
            behaviour from before wavelength-range filtering was
            introduced.  If ``False``, use only the active lines.

        Returns
        -------
        pandas.DataFrame
        """
        self.molecule._ensure_data_loaded()

        if full_range:
            raw = self.molecule._raw_lines_data
        else:
            raw = self.molecule._get_active_raw_data()

        if raw is None or len(raw) == 0:
            return pd.DataFrame(columns=list(self._TABLE_COLUMNS))

        n_rows = len(raw)

        # ---- intensity / tau mapping --------------------------------
        if not full_range:
            # Active range: computed arrays align 1-to-1 with *raw*
            tau_col = self._tau if self._tau is not None else np.full(n_rows, np.nan)
            intens_col = self._intensity if self._intensity is not None else np.full(n_rows, np.nan)
        else:
            # Full range: compute intensity for ALL lines so that the
            # population diagram (and similar consumers) see the same
            # data they did before wavelength-range filtering existed.
            intens_col, tau_col = self._compute_full_range_intensity(raw)

        return pd.DataFrame({
            'lev_up': raw['lev_up'],
            'lev_low': raw['lev_low'],
            'lam': raw['lam'],
            'tau': tau_col,
            'intens': intens_col,
            'a_stein': raw['a_stein'],
            'e_up': raw['e_up'],
            'e_low': raw['e_low'],
            'g_up': raw['g_up'].astype(int),
            'g_low': raw['g_low'].astype(int),
        })

    def _compute_full_range_intensity(self, raw) -> tuple:
        """Compute intensity and tau for all lines in the HITRAN file.

        Creates a temporary unfiltered MoleculeLineList and Intensity
        object, then runs the standard calculation pipeline with the
        same physical parameters (t_kin, n_mol, dv) stored on this
        instance.  This avoids duplicating any calculation logic.

        Parameters
        ----------
        raw : structured np.ndarray
            The full ``_raw_lines_data`` array from the parent
            MoleculeLineList.

        Returns
        -------
        tuple of (np.ndarray, np.ndarray)
            ``(intensity, tau)`` arrays with length ``len(raw)``.
            Returns arrays filled with NaN if no parameters have
            been set yet.
        """
        n_rows = len(raw)
        if self._t_kin is None or self._n_mol is None or self._dv is None:
            return np.full(n_rows, np.nan), np.full(n_rows, np.nan)

        # Create a temporary unfiltered MoleculeLineList and Intensity
        # so the full calculation pipeline runs over all lines.
        from .MoleculeLineList import MoleculeLineList as _MLL
        tmp_mll = _MLL(
            molecule_id=self.molecule.molecule_id,
            filename=self.molecule._filename,
            wavelength_range=None,
        )
        tmp_intens = Intensity(tmp_mll, wavelength_range=None)
        try:
            tmp_intens.calc_intensity(
                t_kin=self._t_kin,
                n_mol=self._n_mol,
                dv=self._dv,
                method="curve_growth",
            )
            return (tmp_intens._intensity, tmp_intens._tau)
        except Exception:
            return np.full(n_rows, np.nan), np.full(n_rows, np.nan)

    # ---- backward-compatible aliases --------------------------------

    @property
    def get_table(self) -> "pandas.DataFrame":
        """pd.DataFrame: Line data for the active wavelength range.

        .. deprecated::
            Prefer ``build_table()`` which accepts ``full_range`` and
            ``wavelength_range`` keyword arguments.  This property is
            kept for backward compatibility.
        """
        return self.build_table()

    @property
    def table(self) -> "pandas.DataFrame":
        """Alias for ``build_table()`` (active range only)."""
        return self.build_table()

    def _repr_html_(self) -> Optional[str]:
        # noinspection PyProtectedMember
        pd = _get_pandas()
        return self.build_table()._repr_html_() if pd is not None else None