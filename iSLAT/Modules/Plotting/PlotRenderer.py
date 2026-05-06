from typing import Optional, List, Dict, Any, Tuple, Union, TYPE_CHECKING
#from matplotlib import lines
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
#from matplotlib.collections import PolyCollection

from matplotlib.axes import Axes
from .BasePlot import BasePlot
from .LineInspectionPlot import LineInspectionPlot
from .PopulationDiagramPlot import PopulationDiagramPlot
from .SpectrumPanel import SpectrumPanel

# Import debug configuration
try:
    from iSLAT.Modules.Debug import debug_config, DebugLevel
except ImportError:
    # Fallback if debug module is not available
    class DebugLevel:
        NONE = 0
        ERROR = 1
        WARNING = 2
        INFO = 3
        VERBOSE = 4
        TRACE = 5
    
    class FallbackDebugConfig:
        def __init__(self):
            self.level = DebugLevel.WARNING
        
        def should_log(self, component, level):
            return level <= self.level
        
        def log(self, component, level, message, **kwargs):
            if self.should_log(component, level):
                print(f"[{component.upper()}] {message}")
        
        def verbose(self, component, message, **kwargs):
            self.log(component, DebugLevel.VERBOSE, message, **kwargs)
        
        def info(self, component, message, **kwargs):
            self.log(component, DebugLevel.INFO, message, **kwargs)
        
        def warning(self, component, message, **kwargs):
            self.log(component, DebugLevel.WARNING, message, **kwargs)
        
        def error(self, component, message, **kwargs):
            self.log(component, DebugLevel.ERROR, message, **kwargs)
    
    debug_config = FallbackDebugConfig()

# Import actual data types for proper type hinting
if TYPE_CHECKING:
    from iSLAT.Modules.DataTypes.Molecule import Molecule
    from iSLAT.Modules.DataTypes.MoleculeDict import MoleculeDict
    #from iSLAT.Modules.DataTypes.MoleculeLineList import MoleculeLineList
    from iSLAT.Modules.DataTypes.MoleculeLine import MoleculeLine
    #from iSLAT.Modules.DataTypes.Intensity import Intensity
    #from iSLAT.Modules.DataTypes.Spectrum import Spectrum

class PlotRenderer(BasePlot):
    """
    Handles all plot rendering and visual updates for the iSLAT spectroscopy tool.
    
    Inherits from :class:`BasePlot` so that all shared helpers
    (``_plot_observed_spectrum``, ``_plot_summed_spectrum``,
    ``_get_theme_value``, ``get_molecule_spectrum_data``, etc.) are
    available as regular instance / inherited methods instead of
    being called through error-prone duck-typing.
    
    This class provides comprehensive rendering of:
    - Observed spectrum data with error bars
    - Individual molecule model spectra leveraging molecule caching
    - Summed model spectra
    - Population diagrams using molecule cached data
    - Line inspection plots
    - Saved line markers
    
    Features:
    - Efficient molecule visibility filtering
    - Direct use of molecule cached calculations for optimal performance
    - Memory-conscious plotting with line limit management
    - Batch operations for better performance with many molecules
    - Zero cache conflicts by relying entirely on molecule caching
    """
    
    def __init__(self, plot_manager: Any) -> None:
        # Initialise BasePlot with the pre-existing GUI figure and theme.
        # ``generate_plot`` is not meaningful for PlotRenderer (it renders
        # on demand), so we provide a no-op implementation.
        super().__init__(fig=plot_manager.fig, theme=plot_manager.theme)

        self.plot_manager = plot_manager
        self.islat = plot_manager.islat
        
        self.ax1: Axes = plot_manager.ax1
        self.ax2: Axes = plot_manager.ax2
        self.ax3: Axes = plot_manager.ax3
        self.canvas = plot_manager.canvas
        
        self.model_lines: List[Line2D] = []
        self.active_lines: List[Line2D] = []
        
        self.render_out = False
        
        # Reusable standalone plot delegates (render into GUI axes)
        self._line_inspection_plot: Optional[LineInspectionPlot] = None
        self._population_diagram_plot: Optional[PopulationDiagramPlot] = None
        
        # Population diagram caching
        self._pop_diagram_molecule = None
        self._pop_diagram_cache_key = None
        
        # Simplified stats - only for performance monitoring, no data caching
        self._plot_stats = {
            'renders_count': 0,
            'molecules_rendered': 0
        }
    
    # ------------------------------------------------------------------
    # BasePlot abstract method — PlotRenderer renders on demand, not via
    # a single ``generate_plot`` call, so this is intentionally a no-op.
    # ------------------------------------------------------------------
    def generate_plot(self, **kwargs) -> None:  # pragma: no cover
        """No-op.  PlotRenderer renders incrementally via dedicated methods."""
        pass

    # Helper methods for common operations — delegate to BasePlot where possible
    @staticmethod
    def _get_molecule_display_name(molecule: 'Molecule') -> str:
        """Get display name for a molecule (delegates to :class:`BasePlot`)."""
        return BasePlot.get_molecule_display_name(molecule)
    
    # _get_theme_value is inherited from BasePlot

    def _get_molecule_color(self, molecule: 'Molecule') -> str:
        """Get color for molecule from theme or molecule properties.

        Extends the static :meth:`BasePlot.get_molecule_color` with
        theme-based colour assignment so that the GUI can assign
        per-molecule colours from the theme palette.
        """
        mol_name = self._get_molecule_display_name(molecule)
        
        # Check molecule's own color first
        if hasattr(molecule, 'color') and molecule.color:
            return molecule.color
            
        # Then check theme colors
        molecule_colors = self._get_theme_value('molecule_colors', {})
        if mol_name in molecule_colors:
            molecule.color = molecule_colors[mol_name]
            return molecule_colors[mol_name]
            
        # Default fallback
        molecule.color = self._get_theme_value('default_molecule_color', 'blue')
        return self._get_theme_value('default_molecule_color', 'blue')

    def clear_model_lines(self, ax: "Axes" = None, lines: List["Line2D"] = None, do_clear_self: bool = True) -> None:
        """Clear only the model spectrum lines from the main plot"""
        if ax is None:
            ax = self.ax1
        if lines is None:
            lines = self.model_lines
        
        for line in lines:
            if line in ax.lines:
                line.remove()
        
        if do_clear_self:
            self.model_lines.clear()

    def remove_molecule_lines(self, molecule_name: str, ax: "Axes" = None, lines: List["Line2D"] = None, update_legend: bool = True) -> None:
        """
        Remove lines associated with a specific molecule
        
        Parameters
        ----------
        molecule_name : str
            Name of molecule whose lines should be removed
        ax : Axes, optional
            Axes to remove from (default: ax1)
        lines : List[Line2D], optional
            List of lines to search in (default: self.model_lines)
        update_legend : bool, optional
            Whether to update legend after removal (default: True)
        """
        print(f"removing lines from {molecule_name}")

        if ax is None:
            ax = self.ax1
        if lines is None:
            lines = self.model_lines

        lines_to_remove = []
        for line in lines:
            # Check if line belongs to this molecule (by stored metadata first)
            if hasattr(line, '_molecule_name') and line._molecule_name == molecule_name:
                print(f"found lines matching with {line._molecule_name} (metadata)")
                lines_to_remove.append(line)

        # Remove lines from plot and list
        for line in lines_to_remove:
            if line in ax.lines:
                line.remove()
            if line in lines:
                lines.remove(line)
        
        # Only update legend if requested
        if update_legend and lines_to_remove:
            BasePlot._update_legend(ax)
            # Respect the legend toggle state from the controller
            if hasattr(self, 'plot_manager') and hasattr(self.plot_manager, 'legend_toggle'):
                leg = ax.get_legend()
                if leg is not None:
                    leg.set_visible(self.plot_manager.legend_toggle)
    
    def set_molecule_visibility(self, molecule_name: str, visible: bool, ax: "Axes" = None, lines: List["Line2D"] = None) -> bool:
        """
        Toggle molecule line visibility without removing/recreating (fastest method).
        
        This is much faster than remove/recreate as it just changes visibility state.
        
        Parameters
        ----------
        molecule_name : str
            Name of molecule to show/hide
        visible : bool
            True to show, False to hide
        ax : Axes, optional
            Axes containing the lines (default: ax1)
        lines : List[Line2D], optional
            List of lines to search in (default: self.model_lines)
            
        Returns
        -------
        bool
            True if lines were found and updated, False otherwise
        """
        if ax is None:
            ax = self.ax1
        if lines is None:
            lines = self.model_lines
        
        found_lines = False
        for line in lines:
            if hasattr(line, '_molecule_name') and line._molecule_name == molecule_name:
                line.set_visible(visible)
                found_lines = True
        
        return found_lines
    
    def _render_fit_results_in_line_inspection(self, fit_result: Any, xmin: float, xmax: float, max_y: float) -> None:
        """Helper method to render fit results in the line inspection plot."""
        
        # Clear old fit results if setting is enabled (but preserve data and molecule plots)
        if self._should_clear_old_fits():
            self._clear_old_fit_results_in_range(xmin, xmax)
        
        try:
            gauss_fit, fitted_wave, fitted_flux = fit_result
            if gauss_fit is not None and fitted_wave is not None and fitted_flux is not None:
                # Filter fit data to range
                fit_mask = (fitted_wave >= xmin) & (fitted_wave <= xmax)
                if np.any(fit_mask):
                    fit_line = self.ax2.plot(fitted_wave[fit_mask], fitted_flux[fit_mask], 
                                color='red', linewidth=1, label='Total Fit', linestyle='--')[0]
                    # Mark as fit result for future removal
                    fit_line._islat_fit_result = True
                    
                    # Check if this is a multi-component fit by looking at the fit result structure
                    if hasattr(gauss_fit, 'params') and gauss_fit.params:
                        # Count components by looking for numbered prefixes (g1_, g2_, etc.)
                        component_prefixes = set()
                        for param_name in gauss_fit.params:
                            if '_' in param_name:
                                prefix = param_name.split('_')[0] + '_'
                                if prefix.startswith('g') and prefix[1:-1].isdigit():
                                    component_prefixes.add(prefix)
                        
                        # If multi-component, plot individual components
                        if len(component_prefixes) > 1:
                            try:
                                components = gauss_fit.eval_components(x=fitted_wave[fit_mask])
                                for i, prefix in enumerate(sorted(component_prefixes)):
                                    if prefix in components:
                                        component_flux = components[prefix]
                                        comp_line = self.ax2.plot(fitted_wave[fit_mask], component_flux, 
                                                    linestyle='--', linewidth=1, 
                                                    label=f"Component {i+1}")[0]
                                        # Mark as fit result for future removal
                                        comp_line._islat_fit_result = True
                            except Exception as e:
                                debug_config.warning("plot_renderer", f"Could not plot fit components: {e}")
                        else:
                            # Single component fit, fill uncertainty area
                            dely = gauss_fit.eval_uncertainty(sigma = self.islat.user_settings.get('fit_line_uncertainty', 1.0))
                            fill_collection = self.ax2.fill_between(fitted_wave, fitted_flux - dely, fitted_flux + dely,
                                                color='gray', alpha=0.3, label=r'3-$\sigma$ uncertainty band')
                            # Mark as fit result for future removal
                            fill_collection._islat_fit_result = True

        except Exception as e:
            debug_config.warning("plot_renderer", f"Could not render fit results: {e}")
        
        handles, labels = self.ax2.get_legend_handles_labels()
        if handles:
            self.ax2.legend()
            # Respect the legend toggle state from the controller
            if hasattr(self, 'plot_manager') and hasattr(self.plot_manager, 'legend_toggle'):
                leg = self.ax2.get_legend()
                if leg is not None:
                    leg.set_visible(self.plot_manager.legend_toggle)
        # Don't call canvas.draw_idle() here - let caller batch it
    
    def _should_clear_old_fits(self) -> bool:
        """Check if old fit results should be cleared when making new selections."""
        try:
            if hasattr(self.islat, 'user_settings'):
                return self.islat.user_settings.get('clear_old_fits', True)
            return True  # Default to True if setting not found
        except Exception as e:
            debug_config.warning("plot_renderer", f"Could not check clear_old_fits setting: {e}")
            return True
    
    def _clear_old_fit_results_in_range(self, xmin: float, xmax: float) -> None:
        """Clear old fit results that overlap with the new selection range."""
        # Remove existing fit lines from the line inspection plot
        lines_to_remove = []
        for line in self.ax2.lines:
            if hasattr(line, '_islat_fit_result') or (hasattr(line, 'get_label') and 
                                                     line.get_label() and 
                                                     ('Fit' in line.get_label() or 'Component' in line.get_label())):
                # Check if the fit line overlaps with the current selection
                if hasattr(line, 'get_xdata'):
                    line_xdata = line.get_xdata()
                    if len(line_xdata) > 0:
                        line_xmin = np.min(line_xdata)
                        line_xmax = np.max(line_xdata)
                        # Check for overlap
                        if (line_xmin <= xmax and line_xmax >= xmin):
                            lines_to_remove.append(line)
        
        # Remove existing fit collections (fill_between objects) from the line inspection plot
        collections_to_remove = []
        for collection in self.ax2.collections:
            if hasattr(collection, '_islat_fit_result') or (hasattr(collection, 'get_label') and 
                                                           collection.get_label() and 
                                                           ('uncertainty' in collection.get_label().lower() or 'sigma' in collection.get_label().lower())):
                # For collections, we need to check their path bounds
                try:
                    paths = collection.get_paths()
                    if paths:
                        # Get bounds from the first path
                        bounds = paths[0].get_extents()
                        coll_xmin = bounds.xmin
                        coll_xmax = bounds.xmax
                        # Check for overlap
                        if (coll_xmin <= xmax and coll_xmax >= xmin):
                            collections_to_remove.append(collection)
                except:
                    # If we can't determine bounds, remove it to be safe
                    collections_to_remove.append(collection)
        
        # Remove the overlapping fit lines
        for line in lines_to_remove:
            line.remove()
            debug_config.trace("plot_renderer", f"Removed old fit result line: {line.get_label()}")
            
        # Remove the overlapping fit collections
        for collection in collections_to_remove:
            collection.remove()
            debug_config.trace("plot_renderer", f"Removed old fit result collection: {collection.get_label()}")
    
    def get_molecule_spectrum_data(self, molecule: 'Molecule', wave_data: np.ndarray,
                                    interpolate_to_input: bool = False,
                                    target_wavelengths: Optional[np.ndarray] = None,
                                    ) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
        """
        Get spectrum data from molecule's caching system.
        
        Thin override that adds debug logging around the inherited
        :meth:`BasePlot.get_molecule_spectrum_data` static method.
        """
        if molecule is None or wave_data is None:
            return None, None
        
        result_wavelengths, result_flux = BasePlot.get_molecule_spectrum_data(
            molecule, wave_data,
            interpolate_to_input=interpolate_to_input,
            target_wavelengths=target_wavelengths,
        )
        
        if result_wavelengths is not None and result_flux is not None and len(result_flux) > 0:
            debug_config.verbose("plot_renderer",
                               f"Retrieved flux data for {self._get_molecule_display_name(molecule)}",
                               data_points=len(result_flux))
            return result_wavelengths, result_flux
        
        debug_config.warning("plot_renderer", f"No flux data available for {self._get_molecule_display_name(molecule)}")
        return None, None
    
    def render_individual_molecule_spectrum(self, molecule: 'Molecule', wave_data: np.ndarray, 
                                         plot_name: Optional[str] = None, subplot: Optional[Axes] = None,
                                         update_legend: bool = True) -> bool:
        """
        Render a single molecule spectrum using the molecule's cached data.
        
        Works with any axes - searches for existing lines on the target axes and updates
        them in place, or creates new ones. Does not depend on self.model_lines.
        
        Parameters
        ----------
        molecule : Molecule
            Molecule object with internal caching
        wave_data : np.ndarray
            Wavelength array
        plot_name : Optional[str]
            Custom name for plotting
        subplot : Optional[Axes]
            Subplot to render on (default: ax1)
        update_legend : bool
            Whether to update the legend after rendering (default: True)
            Set to False when rendering multiple molecules to update legend once at the end
            
        Returns
        -------
        bool
            True if successfully plotted, False otherwise
        """
        try:
            plot = subplot if subplot else self.ax1
            # Increment render stats
            self._plot_stats['renders_count'] += 1
            
            molecule_name = plot_name or self._get_molecule_display_name(molecule)
            mol_identifier = getattr(molecule, 'name', molecule_name)
            
            # Search for existing line with this molecule name directly on the target axes
            existing_line = None
            for line in plot.lines:
                if (hasattr(line, '_molecule_name') and 
                    line._molecule_name == mol_identifier):
                    existing_line = line
                    break
            
            # When match_spectral_sampling is enabled, resample each molecule's
            # model onto the data pixel grid (corrected for stellar RV).
            # get_matched_sampling_wavelengths expects observer-frame
            # wavelengths; use wave_data_original when available.
            use_interp = False
            target_wave = None
            if hasattr(self, 'islat') and hasattr(self.islat, 'molecules_dict'):
                mol_dict = self.islat.molecules_dict
                if hasattr(mol_dict, 'get_matched_sampling_wavelengths'):
                    wave_obs = (self.islat.wave_data_original
                                if hasattr(self.islat, 'wave_data_original')
                                else wave_data)
                    if wave_obs is not None:
                        use_interp, target_wave = mol_dict.get_matched_sampling_wavelengths(wave_obs)
                        if not use_interp:
                            target_wave = None

            # Get spectrum data directly from molecule's caching system
            plot_lam, plot_flux = self.get_molecule_spectrum_data(
                molecule, wave_data,
                interpolate_to_input=use_interp,
                target_wavelengths=target_wave,
            )
            
            if plot_lam is None or plot_flux is None:
                print(f"No spectrum data available for {molecule_name}")
                return False
            
            # Check if we actually have meaningful flux data
            if len(plot_flux) == 0 or np.all(plot_flux == 0):
                print(f"Spectrum data is empty or all zeros for {molecule_name}")
                return False
            
            # Get molecule properties
            color = self._get_molecule_color(molecule)
            label = getattr(molecule, 'displaylabel', molecule_name)
            
            lw = self._get_theme_value("model_linewidth", 1)
            alpha = self._get_theme_value("model_alpha", 1)
            
            # Update existing line if found, otherwise create new
            if existing_line is not None:
                # Update existing line data - much faster than remove/recreate
                existing_line.set_data(plot_lam, plot_flux)
                existing_line.set_color(color)
                existing_line.set_alpha(alpha)
                existing_line.set_linewidth(lw)
                existing_line.set_label(label)
                line = existing_line
            else:
                # Create new line only if it doesn't exist on this axes
                line, = plot.plot(
                    plot_lam,
                    plot_flux,
                    linestyle='--',
                    color=color,
                    alpha=alpha,
                    linewidth=lw,
                    label=label,
                    zorder=self._get_theme_value("zorder_model", 3)
                )
                
                # Store molecule name in line metadata for selective removal
                line._molecule_name = mol_identifier
                
                # Only track in self.model_lines if this is the main plot (ax1)
                if plot == self.ax1:
                    self.model_lines.append(line)
    
            # Only update legend if requested (batch operations can skip this)
            if update_legend:
                BasePlot._update_legend(plot)
                # Respect the legend toggle state from the controller
                if hasattr(self, 'plot_manager') and hasattr(self.plot_manager, 'legend_toggle'):
                    leg = plot.get_legend()
                    if leg is not None:
                        leg.set_visible(self.plot_manager.legend_toggle)
            
            self._plot_stats['molecules_rendered'] += 1
            
            return True
            
        except Exception as e:
            print(f"Error plotting molecule {self._get_molecule_display_name(molecule)}: {e}")
            import traceback
            traceback.print_exc()
            return False

    def plot_fitted_saved_lines(self, fit_data, ax: Optional[plt.Axes] = None) -> None:
        """
        Plot the fitted saved lines on the provided axes using flux integral calculation.

        Parameters
        ----------
        fit_data : Any
            The data to plot.
        ax : Optional[plt.Axes]
            The axes to plot on. If None, uses the current axes.
        """
        if ax is None:
            ax = plt.gca()

        #print(f"fit data: {fit_data}")

        # Unpack the fit_data tuple
        gauss_fits, fitted_waves, fitted_fluxes = fit_data
        
        # Iterate through each fit
        for i, (gauss_fit, fitted_wave, fitted_flux) in enumerate(zip(gauss_fits, fitted_waves, fitted_fluxes)):
            # Skip failed fits where results are None
            if gauss_fit is None or fitted_wave is None or fitted_flux is None:
                continue
            
            lam_min = np.min(fitted_wave)
            lam_max = np.max(fitted_wave)
            
            # Plot the fit result using SpectrumPanel helper
            uncertainty_sigma = self.islat.user_settings.get('fit_line_uncertainty', 3.0)
            SpectrumPanel.plot_gaussian_fit(
                ax, gauss_fit, fitted_wave, fitted_flux,
                color='lime', uncertainty_sigma=uncertainty_sigma,
            )
            
            # plot the xmin and xmax for each line
            ax.vlines([lam_min, lam_max], -2, 10, colors='lime', alpha=0.5)

            #ax.legend()
            # Note: caller is responsible for calling canvas.draw_idle()