"""
SpectralPanel -- Abstract base class for a single panel of spectral data.

Represents a rectangular axes view of observed spectral data within a
bounded wavelength range ``[xmin, xmax]``.  Concrete subclasses (e.g.
:class:`LineInspectionPlot`) implement the actual rendering logic, but
every spectral panel shares the same fundamental contract:

* Store observed *wave_data* / *flux_data* (and optional *error_data*).
* Define a wavelength window via *xmin* / *xmax*.
* Be embeddable in a larger figure via an external *ax* parameter.
* Provide a :meth:`generate_plot` that draws content onto its axes.

The :class:`StackedSpectralPanel` composer creates and manages many
``SpectralPanel`` instances to build multi-panel stacked figures.
"""

from abc import abstractmethod
from typing import Optional, Tuple, Dict, Any, TYPE_CHECKING

import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure as MplFigure

from .BasePlot import BasePlot

if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict


class SpectralPanel(BasePlot):
    """
    Abstract base for a single spectral-data panel.

    Subclasses must implement :meth:`generate_plot` (inherited from
    :class:`BasePlot`) and optionally override :meth:`compute_ylim` to
    customise how vertical limits are determined.

    Parameters
    ----------
    wave_data : np.ndarray
        Full observed wavelength array.
    flux_data : np.ndarray
        Full observed flux array matching *wave_data*.
    xmin : float
        Left wavelength bound of the panel.
    xmax : float
        Right wavelength bound of the panel.
    error_data : np.ndarray, optional
        Flux uncertainties matching *wave_data*.
    molecules : MoleculeDict, optional
        Molecule collection -- visible molecules are plotted.
    figsize : tuple, optional
        Figure size.  Defaults to ``(10, 4)``.
    ax : Axes, optional
        Pre-existing axes for embedding in a larger figure.  When
        provided the panel does **not** create its own figure.
    **kwargs
        Forwarded to :class:`BasePlot`.
    """

    def __init__(
        self,
        wave_data: np.ndarray,
        flux_data: np.ndarray,
        xmin: float,
        xmax: float,
        error_data: Optional[np.ndarray] = None,
        molecules: Optional["MoleculeDict"] = None,
        figsize: Optional[Tuple[float, float]] = None,
        ax: Optional[Axes] = None,
        **kwargs,
    ):
        super().__init__(figsize=figsize or (10, 4), **kwargs)
        self.wave_data = np.asarray(wave_data)
        self.flux_data = np.asarray(flux_data)
        self.xmin = xmin
        self.xmax = xmax
        self.error_data = (
            np.asarray(error_data) if error_data is not None else None
        )
        self.molecules = molecules
        self._external_ax: Optional[Axes] = ax
        self._ax: Optional[Axes] = None

    # ------------------------------------------------------------------
    # Panel range property
    # ------------------------------------------------------------------
    @property
    def xlim(self) -> Tuple[float, float]:
        """Wavelength limits of this panel."""
        return (self.xmin, self.xmax)

    # ------------------------------------------------------------------
    # Axes access
    # ------------------------------------------------------------------
    @property
    def ax(self) -> Optional[Axes]:
        """The matplotlib Axes used by this panel (may be *None* before
        :meth:`generate_plot` is called)."""
        return self._ax

    # ------------------------------------------------------------------
    # Masking helper
    # ------------------------------------------------------------------
    def get_panel_mask(self) -> np.ndarray:
        """Boolean mask selecting data points inside ``[xmin, xmax]``."""
        return (self.wave_data >= self.xmin) & (self.wave_data <= self.xmax)

    def get_panel_data(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
        """Return ``(wave, flux, error)`` for points inside the panel range.

        *error* is ``None`` when :attr:`error_data` was not provided.
        """
        mask = self.get_panel_mask()
        err = self.error_data[mask] if self.error_data is not None else None
        return self.wave_data[mask], self.flux_data[mask], err

    # ------------------------------------------------------------------
    # Default y-limit computation (overridable)
    # ------------------------------------------------------------------
    def compute_ylim(
        self,
        ymax_factor: float = 0.2,
    ) -> Tuple[float, float]:
        """Compute ``(ymin, ymax)`` for this panel.

        The default implementation examines the observed flux within the
        panel range and adds *ymax_factor* headroom.  Subclasses may
        override this to incorporate model data, residuals, etc.

        Parameters
        ----------
        ymax_factor : float
            Fractional padding above the peak flux (0.2 = 20 %).

        Returns
        -------
        tuple[float, float]
        """
        mask = self.get_panel_mask()
        if np.any(mask):
            peak = float(np.nanmax(self.flux_data[mask]))
            return (-0.005, peak + peak * ymax_factor)
        return (-0.005, 0.1)

    # ------------------------------------------------------------------
    # Convenience: change range and regenerate
    # ------------------------------------------------------------------
    def set_range(self, xmin: float, xmax: float) -> None:
        """Update the panel wavelength range and regenerate the plot."""
        self.xmin = xmin
        self.xmax = xmax
        self.generate_plot()

    # ------------------------------------------------------------------
    # Resolve / prepare axes
    # ------------------------------------------------------------------
    def _resolve_axes(self) -> Axes:
        """Return an Axes to draw on.

        If an external axes was supplied at construction time it is
        returned directly (and cleared).  Otherwise a new figure + axes
        pair is created.
        """
        if self._external_ax is not None:
            self._ax = self._external_ax
        else:
            self._ensure_figure()
            self.fig.clf()
            self._ax = self.fig.add_subplot(111)
        self._ax.clear()
        return self._ax
