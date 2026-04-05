# -*- coding: utf-8 -*-

"""
This class Spectrum creates a spectrum from intensity instances

* The same algorithm as in the Fortran 90 code are used. This is explained in the appendix of Banzatti et al. 2012.
* After a Spectrum instance is created with the wavelength grid and resolution, different intensity components can be
  added. When they are added, only the area-scaled intensities per line are stored. The convolution with the resolution
  of the instrument is done at the end (a la lazy evaluation at the first call to retrieve the flux). This improves the
  performance very much when different components of the same molecule are added.
* The get the spectral convolution of the raw lines acceptably quick, only a range of about 15 sigma around the line
  center is evaluated. This is the same trick as has been used in the Fortran 90 code.
* Only read access to the fields is granted through properties

- 01/06/2020: SB, initial version

"""

import numpy as np
from typing import Optional, List, Dict, Any, Union

import iSLAT.Constants as c

from ._pandas_import import get_pandas as _get_pandas
from ._mixins import CacheStatsMixin, WavelengthRangeMixin

class Spectrum(CacheStatsMixin, WavelengthRangeMixin):
    """
    Spectrum class for creating and managing spectral data from intensity instances.
    
    Implements lazy evaluation for convolution and caching for performance optimization.
    """
    
    # Pre-computed class-level constants for performance
    _FWHM_TO_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))  # ≈ 2.3548
    _INV_SQRT_2PI = 1.0 / np.sqrt(2.0 * np.pi)  # ≈ 0.3989
    _AU_PC_RATIO_SQ = (c.ASTRONOMICAL_UNIT_CM / c.PARSEC_CM) ** 2  # (AU/pc)²
    _FLUX_JY_FACTOR = 1e19 / c.SPEED_OF_LIGHT_CGS  # Conversion factor for Jy
    
    __slots__ = (
        '_lam_min', '_lam_max', '_dlambda', '_R', '_R_func', '_distance',
        '_lamgrid', '_flux', '_flux_jy', '_I_arrays', '_lam_arrays', '_tau_arrays',
        '_components', '_flux_valid', '_convolution_cache',
        '_kernel_cache', '_cache_stats', '_n_grid_points', '_unique_cache',
        '_wavelength_range'
    )
    
    def __init__(self, lam_min: float = None, lam_max: float = None, 
                 dlambda: float = None, R: float = None, distance: float = None,
                 wavelength_range: Optional[tuple] = None,
                 R_func: Optional[Any] = None):
        """Initialize a spectrum class and prepare it to add intensity components

        Parameters
        ----------
        lam_min: float
            Lower border of the spectrum in micron
        lam_max: float
            Upper border of the spectrum in micron
        dlambda: float
            Resolution of the spectrum to calculate in micron
        R: float
            Spectral resolution of the instrument as a scalar
            R = lambda/delta_lambda.  Used when *R_func* is not provided.
        distance: float
            Distance to the disk in pc
        wavelength_range: tuple of (float, float), optional
            Wavelength range ``(lam_min, lam_max)`` in microns.  When set,
            only lines within this range are included in the convolution.
            Defaults to ``(lam_min, lam_max)`` when not provided.
        R_func: callable, optional
            A function ``R_func(wavelength_um)`` that returns the
            wavelength-dependent resolving power *R* for one or more
            wavelengths (in microns).  When provided, *R_func* is used
            instead of the scalar *R* during convolution, giving each
            spectral line a kernel width that matches the local
            instrument resolution (e.g. JWST MIRI MRS).

            The callable must accept a float or 1-D ``numpy.ndarray``
            and return a value of the same shape.
        """

        # assure valid lambda grid range
        if lam_min >= lam_max:
            raise ValueError('lam_min must be < lam_max')

        # store parameters
        self._lam_min = lam_min
        self._lam_max = lam_max
        self._dlambda = dlambda
        self._R = R
        self._R_func = R_func
        self._distance = distance
        self._wavelength_range = wavelength_range if wavelength_range is not None else (lam_min, lam_max)

        # create wavelength grid with pre-calculated size
        n_points = int(1 + (lam_max - lam_min) / dlambda)
        self._lamgrid = np.linspace(lam_min, lam_max, n_points)
        self._n_grid_points = n_points  # Cache grid length

        # flux array (cached result)
        self._flux = None
        self._flux_jy = None

        # Use lists of numpy arrays for efficient accumulation
        # (concatenating arrays at convolution time is faster than extending lists)
        self._I_arrays = []
        self._lam_arrays = []
        self._tau_arrays = []

        # list with the different intensity components building up the spectrum
        self._components = []
        
        # Cache invalidation flag and caching system
        self._flux_valid = False
        self._convolution_cache = {}
        self._kernel_cache = {}
        self._init_cache_stats()
        self._unique_cache = None  # Cached (key, lam, index_wavelength) from np.unique

    def reset(self):
        """Clear all accumulated intensity components, keeping the grid and kernel caches.

        This is much cheaper than constructing a new Spectrum when only the
        physical parameters (T, N, radius) change between evaluations,
        because the wavelength grid and convolution kernels are preserved.
        """
        self._I_arrays.clear()
        self._lam_arrays.clear()
        self._tau_arrays.clear()
        self._components.clear()
        self._flux = None
        self._flux_jy = None
        self._flux_valid = False
        self._unique_cache = None

    def add_intensity(self, intensity, dA: float):
        """Adds an intensity component to the spectrum

        Parameters
        ----------
        intensity: Intensity
            Intensity structure to add to the spectrum
        dA: float
            Area of the component in au**2
        """

        # Invalidate cached flux when adding new intensity
        self._invalidate_flux_cache()

        # 1. get intensity and wavelength of the lines
        I_all = intensity.intensity
        # Get wavelengths from the new MoleculeLineList structure
        lam_all = intensity.molecule.get_wavelengths()
        
        # Check if we have valid data to process
        if I_all is None or lam_all is None:
            print(f"Warning: No intensity or wavelength data available for molecule {intensity.molecule.name}")
            return
            
        if len(lam_all) == 0 or len(I_all) == 0:
            # No lines to add, just return without error
            return

        # 2. select only lines within the selected wavelength range using vectorized operations
        wr_min, wr_max = self._wavelength_range
        # Use the minimum R across the range for border padding to capture
        # all lines whose wings could contribute to the grid.
        if self._R_func is not None:
            R_border = float(np.nanmin(np.asarray(
                self._R_func(np.array([wr_min, wr_max])), dtype=float)))
        else:
            R_border = self._R
        select_border = 100 * wr_max / R_border
        mask = (lam_all >= wr_min - select_border) & (lam_all <= wr_max + select_border)
        
        if not np.any(mask):
            return  # No lines in range

        # 3. scale for area in au**2 using vectorized operations
        selected_wavelengths = lam_all[mask]
        selected_intensities = I_all[mask] * dA

        # 4. append numpy arrays directly (much faster than list.extend with numpy arrays)
        self._I_arrays.append(selected_intensities)
        self._lam_arrays.append(selected_wavelengths)
        self._tau_arrays.append(intensity.tau[mask])

        # 5. append to components
        self._components.append({'name': intensity.molecule.name, 'fname': getattr(intensity.molecule, 'fname', ''),
                                't_kin': intensity.t_kin, 'n_mol': intensity.n_mol, 'dv': intensity.dv, 'tau': intensity.tau,
                                 'area': dA})
    
    def _invalidate_flux_cache(self):
        """Invalidate cached flux values"""
        self._flux = None
        self._flux_jy = None
        self._flux_valid = False
        self._unique_cache = None
        self._record_cache_invalidation()

    # Number of sigma bins for per-group kernel sizing.
    # Lines are grouped by similar Gaussian width so narrow lines use a tighter
    # kernel instead of the global maximum. More bins = tighter kernels per
    # group (faster), but more loop iterations. 5 is a good balance.
    _N_SIGMA_BINS = 5

    def _convol_flux(self):
        """Internal procedure to carry out the convolution, should never be called directly

        Returns
        -------
        np.ndarray:
            List with fluxes, convolved to the spectral resolution
        """
        
        # Check if we have any intensities to process
        if len(self._I_arrays) == 0:
            # Return a zero flux array if no intensities were added
            return np.zeros_like(self._lamgrid)

        # Concatenate all accumulated arrays at once (much faster than list extend)
        I_array = np.concatenate(self._I_arrays)
        lam_array = np.concatenate(self._lam_arrays)

        # 1. Summarize intensities at the (exactly) same wavelength. This improves
        #    performance because only one convolution kernel needs to be evaluated per
        #    unique line wavelength (independent of how many intensity components share it).
        #
        #    Cache the np.unique result: when the same Spectrum is reused in a fitting
        #    loop (via reset()), the wavelength arrays don't change between calls -- only
        #    the intensities do. Skipping the sort saves ~5-10% of convolution time.
        lam_key = lam_array.data.tobytes()
        if self._unique_cache is not None and self._unique_cache[0] == lam_key:
            lam, index_wavelength = self._unique_cache[1], self._unique_cache[2]
        else:
            lam, index_wavelength = np.unique(lam_array, return_inverse=True)
            self._unique_cache = (lam_key, lam, index_wavelength)

        # np.bincount is ~5-20x faster than np.add.at for scatter-add accumulation
        intens = np.bincount(index_wavelength, weights=I_array,
                             minlength=lam.shape[0]).astype(np.float64)
        
        n_grid = self._n_grid_points

        # 2. Calculate width and normalization of convolution kernel
        # Use wavelength-dependent R(λ) when a callable is provided;
        # otherwise fall back to the fixed scalar R.
        if self._R_func is not None:
            R_per_line = np.atleast_1d(np.asarray(self._R_func(lam), dtype=float))
            fwhm = lam / R_per_line
        else:
            fwhm = lam / self._R
        sigma = fwhm / self._FWHM_TO_SIGMA
        
        # Pre-compute terms that will be reused: norm * intens and 1/(2*sigma^2)
        norm_intens = (self._INV_SQRT_2PI / sigma) * intens  # shape: (n_lines,)
        inv_2sigma_sq = 0.5 / (sigma ** 2)                   # shape: (n_lines,)

        # Grid positions for each unique line wavelength
        lam_grid_position = (n_grid * (lam - self._lam_min) /
                             (self._lam_max - self._lam_min)).astype(np.int32)

        # 3. Per-group kernel sizing
        #    Instead of using the global max sigma for every line (which wastes
        #    computation on narrow lines at short wavelengths), we bin lines by
        #    similar sigma and use a tight kernel range per group.
        flux = np.zeros(n_grid, dtype=np.float64)

        sigma_min = sigma.min()
        sigma_max = sigma.max()

        if sigma_min == sigma_max or lam.shape[0] <= 1:
            # All lines have the same width -- single group, no binning overhead
            sigma_groups = [np.arange(lam.shape[0])]
            sigma_maxes = [sigma_max]
        else:
            # Create logarithmically-spaced bins so each group spans a similar
            # ratio of sigma values (appropriate because sigma is proportional to wavelength).
            bin_edges = np.geomspace(sigma_min * 0.999, sigma_max * 1.001,
                                     num=self._N_SIGMA_BINS + 1)
            group_ids = np.digitize(sigma, bin_edges) - 1  # 0-based bin index
            group_ids = np.clip(group_ids, 0, self._N_SIGMA_BINS - 1)

            unique_groups = np.unique(group_ids)
            sigma_groups = [np.where(group_ids == g)[0] for g in unique_groups]
            sigma_maxes = [sigma[idx].max() for idx in sigma_groups]

        lamgrid = self._lamgrid  # local reference avoids repeated attribute lookup

        for group_idx, group_sigma_max in zip(sigma_groups, sigma_maxes):
            # Tight kernel range for this group
            kernel_range_size = int(15.0 * group_sigma_max / self._dlambda)
            kernel_range = np.arange(-kernel_range_size, kernel_range_size + 1,
                                     dtype=np.int32)

            g_pos = lam_grid_position[group_idx]          # (n_group,)
            g_norm = norm_intens[group_idx]                # (n_group,)
            g_inv2s = inv_2sigma_sq[group_idx]             # (n_group,)
            g_lam = lam[group_idx]                         # (n_group,)

            # 4. 2D broadcasting: compute all (line, offset) pairs
            #    grid_indices shape: (n_group, kernel_size)
            grid_indices = g_pos[:, np.newaxis] + kernel_range[np.newaxis, :]

            # Validity mask for in-bounds indices
            valid_mask = (grid_indices >= 0) & (grid_indices < n_grid)

            # Clip to avoid out-of-bounds; masked values will be zeroed below
            safe_indices = np.clip(grid_indices, 0, n_grid - 1)

            # Wavelength differences -> Gaussian kernel (in-place)
            kernel = lamgrid[safe_indices] - g_lam[:, np.newaxis]  # delta_lam
            np.multiply(kernel, kernel, out=kernel)                 # delta_lam^2
            kernel *= -g_inv2s[:, np.newaxis]                       # -delta_lam^2 / (2*sigma^2)
            np.exp(kernel, out=kernel)                              # exp(...)
            kernel *= g_norm[:, np.newaxis]                         # norm * exp(...)

            # Zero out invalid (out-of-grid) positions
            kernel *= valid_mask

            # 6. Scatter-add via np.bincount (much faster than np.add.at)
            flat_idx = safe_indices.ravel()
            flat_wt = kernel.ravel()
            flux += np.bincount(flat_idx, weights=flat_wt, minlength=n_grid)

        # 7. Scale for distance and correct units for the area
        inv_dist_sq = 1.0 / (self._distance ** 2)
        scaled_flux = flux * (self._AU_PC_RATIO_SQ * inv_dist_sq)

        return scaled_flux

    @property
    def flux(self) -> np.ndarray:
        """np.ndarray: Flux density in erg/s/cm^2/micron"""
        if self._flux is None or not self._flux_valid:
            self._record_cache_miss()
            self._flux = self._convol_flux()
            self._flux_valid = True
        else:
            self._record_cache_hit()
        return self._flux

    @property
    def flux_jy(self) -> np.ndarray:
        """np.ndarray: Flux density in Jy/micron"""
        if self._flux_jy is None:
            flux_data = self.flux  # This triggers flux calculation if needed
            # Use pre-computed conversion factor
            self._flux_jy = flux_data * self._FLUX_JY_FACTOR * (self._lamgrid ** 2)
        return self._flux_jy

    # wavelength_range property provided by WavelengthRangeMixin

    def _on_wavelength_range_changed(self, old, new):
        """Hook called by WavelengthRangeMixin when wavelength_range changes."""
        self._invalidate_flux_cache()

    @property
    def lamgrid(self) -> np.ndarray:
        """np.ndarray: Wavelength grid in micron"""
        return self._lamgrid

    @property
    def components(self) -> List[Dict[str, Any]]:
        """list of dict: Intensity components added to the spectrum"""
        return self._components

    @property
    def R_func(self):
        """callable or None: Wavelength-dependent resolving-power function."""
        return self._R_func

    @R_func.setter
    def R_func(self, value):
        """Set a new R(λ) function.  Invalidates cached flux."""
        if value is not self._R_func:
            self._R_func = value
            self._invalidate_flux_cache()

    def resample_to(self, target_wavelengths: np.ndarray,
                    unit: str = 'jy',
                    rv_shift: float = 0.0,
                    fill: float = 0.0) -> np.ndarray:
        """Flux-conserving resampling of the spectrum onto a target grid.

        Uses the same ``spectres``-style algorithm employed by
        :func:`iSLAT.Modules.DataTypes.Molecule._spectres` to perform
        flux-conserving resampling that correctly handles uneven pixel
        sampling (e.g. JWST MIRI MRS).

        Parameters
        ----------
        target_wavelengths : np.ndarray
            1-D array of desired output wavelength positions (µm).
        unit : str, optional
            ``'jy'`` (default) for Jy, ``'cgs'`` for erg/s/cm²/µm.
        rv_shift : float, optional
            Radial-velocity shift in km/s to apply to the model grid
            **before** resampling.  Positive values red-shift the model
            (i.e. ``λ_obs = λ_rest x (1 + rv / c)``).
        fill : float, optional
            Value used for target pixels that fall outside the model grid.

        Returns
        -------
        np.ndarray
            Resampled flux array with the same length as
            *target_wavelengths*.
        """
        import importlib as _il
        _su = _il.import_module('iSLAT.Modules.DataProcessing.spectral_utils')
        _spectres = _su.spectres

        source_wave = self._lamgrid
        if rv_shift != 0.0:
            doppler = 1.0 + rv_shift / c.SPEED_OF_LIGHT_KMS
            source_wave = source_wave * doppler

        if unit == 'jy':
            source_flux = self.flux_jy
        else:
            source_flux = self.flux

        return _spectres(target_wavelengths, source_wave, source_flux,
                         fill=fill)

    @property
    def get_table(self):
        """pd.DataFrame: Pandas dataframe"""
        pd = _get_pandas()
        # Pre-compute all data to avoid repeated property access
        flux_data = self.flux
        flux_jy_data = self.flux_jy
        
        return pd.DataFrame({
            'lam': self._lamgrid,
            'flux': flux_data,
            'flux_jy': flux_jy_data
        })

    # get_cache_stats() provided by CacheStatsMixin

    def _repr_html_(self):
        # noinspection PyProtectedMember
        return self.get_table._repr_html_()