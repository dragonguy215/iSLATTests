"""
CompositePlot — Abstract base for multi-panel composite plot objects.

A :class:`CompositePlot` owns an ordered collection of :class:`BasePlot`
child panels that share a single figure.  Concrete subclasses supply:

* :meth:`_build_layout` — creates the axes / gridspec layout on the
  figure and returns it.
* :meth:`_create_panels` — instantiates the child :class:`BasePlot`
  objects, attaches axes from the layout, and registers each one via
  :meth:`_register_panel`.

The default :meth:`generate_plot` orchestrates the full build cycle::

    _ensure_figure() → fig.clf() → child_panels.clear()
    → _build_layout() → _create_panels()
    → for each panel: panel.generate_plot()
    → _post_render() → apply_theme_to_figure()

Subclasses that need different sequencing (e.g. borrowed-axes mode in
:class:`MainPlotGrid`, or per-row y-limit normalisation in
:class:`StackedSpectralPanel`) should override :meth:`generate_plot`
directly — the shared child-panel registry and helper methods remain
available regardless.
"""

from abc import abstractmethod
from typing import Any, Dict, Iterator, List, Optional, TYPE_CHECKING

from matplotlib.axes import Axes

from .BasePlot import BasePlot

if TYPE_CHECKING:
    pass  # avoid circular imports; annotations below use string literals

class CompositePlot(BasePlot):
    """Abstract base for composite multi-panel plot objects.

    Adds a named child-panel registry on top of :class:`BasePlot` and
    provides a default :meth:`generate_plot` that sequentially builds the
    layout, creates panels, renders each one, runs a post-render hook, and
    applies the figure theme.

    Parameters
    ----------
    **kwargs
        Forwarded verbatim to :class:`BasePlot`.
    """

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        #: Ordered mapping of ``name → panel`` for every child panel
        #: registered by :meth:`_register_panel`.
        self.child_panels: Dict[str, "BasePlot"] = {}

    # ------------------------------------------------------------------
    # Abstract interface (subclasses must implement)
    # ------------------------------------------------------------------

    @abstractmethod
    def _build_layout(self) -> Any:
        """Create the axes / gridspec layout for this composite figure.

        Called after :meth:`~BasePlot._ensure_figure` and ``fig.clf()``
        in the default :meth:`generate_plot`.  The return value is passed
        directly to :meth:`_create_panels`.

        Returns
        -------
        Any
            Whatever the subclass needs to attach axes to its child
            panels — a tuple of :class:`~matplotlib.axes.Axes`, a
            :class:`~matplotlib.gridspec.GridSpec`, a 2-D NumPy array of
            :class:`~matplotlib.axes.Axes`, etc.
        """
        ...

    @abstractmethod
    def _create_panels(self, layout: Any) -> None:
        """Instantiate child panels and register them.

        Called with the value returned by :meth:`_build_layout`.
        Implementations should create :class:`BasePlot` instances (or
        subclasses), attach the axes from *layout*, and call
        :meth:`_register_panel` for each one.

        .. important::
           Do **not** call ``panel.generate_plot()`` here.  The default
           :meth:`generate_plot` loop handles rendering after this method
           returns.

        Parameters
        ----------
        layout : Any
            The object returned by :meth:`_build_layout`.
        """
        ...

    # ------------------------------------------------------------------
    # Default generate_plot orchestration
    # ------------------------------------------------------------------

    def generate_plot(self, **kwargs) -> None:
        """Build the composite figure.

        The default implementation follows the sequence::

            _ensure_figure → fig.clf → child_panels.clear
            → _build_layout → _create_panels
            → for each panel: panel.generate_plot
            → _post_render → apply_theme_to_figure

        Subclasses with non-standard layouts (e.g. borrowed-axes mode)
        should override this method.
        """
        self._ensure_figure()
        self.fig.clf()
        self.child_panels.clear()
        layout = self._build_layout()
        self._create_panels(layout)
        for panel in self.child_panels.values():
            panel.generate_plot(**kwargs)
        self._post_render(**kwargs)
        self.apply_theme_to_figure()

    # ------------------------------------------------------------------
    # Post-render hook
    # ------------------------------------------------------------------

    def _post_render(self, **kwargs) -> None:
        """Called after all child panels have been rendered.

        The default implementation is a no-op.  Subclasses override this
        to apply cross-panel decorations such as shared axis labels,
        figure-level annotations, or per-panel overlays that depend on
        finalised y-limits.
        """

    # ------------------------------------------------------------------
    # Panel registry helpers
    # ------------------------------------------------------------------

    def _register_panel(self, name: str, panel: "BasePlot") -> None:
        """Add *panel* to the child-panel registry under *name*.

        If *name* already exists the entry is overwritten (same-object
        re-registration is harmless).

        Parameters
        ----------
        name : str
            Unique identifier for the panel (e.g. ``"spectrum"``,
            ``"cell_0_0"``).
        panel : BasePlot
            The panel instance to register.
        """
        self.child_panels[name] = panel

    def get_panel(self, name: str) -> Optional["BasePlot"]:
        """Return the child panel registered as *name*, or *None*.

        Parameters
        ----------
        name : str
            The key passed to :meth:`_register_panel`.
        """
        return self.child_panels.get(name)

    def iter_panels(self) -> Iterator["BasePlot"]:
        """Iterate over all registered child panels in insertion order."""
        return iter(self.child_panels.values())

    def get_axes(self) -> List[Axes]:
        """Return a flat list of all axes held by registered child panels.

        Checks for a ``.ax`` property (as defined by
        :class:`~iSLAT.Modules.Plotting.SpectralPanel.SpectralPanel`)
        and falls back to ``._external_ax`` for other :class:`BasePlot`
        subtypes.  Duplicate axes objects are de-duplicated while
        preserving insertion order.
        """
        seen: set = set()
        axes: List[Axes] = []
        for panel in self.child_panels.values():
            ax = getattr(panel, "ax", None) or getattr(panel, "_external_ax", None)
            if ax is not None and id(ax) not in seen:
                seen.add(id(ax))
                axes.append(ax)
        return axes

    # ------------------------------------------------------------------
    # Convenience refresh
    # ------------------------------------------------------------------

    def refresh(self) -> None:
        """Re-render all registered child panels without rebuilding the layout.

        Each panel's :meth:`~BasePlot.generate_plot` is called in
        registration order.  After all panels have been re-rendered the
        theme is re-applied to the figure.
        """
        for panel in self.child_panels.values():
            panel.generate_plot()
        self.apply_theme_to_figure()