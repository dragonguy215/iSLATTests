import matplotlib
# matplotlib.use('Agg')
matplotlib.use ("TKAgg")
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button, SpanSelector, TextBox, CheckButtons
from matplotlib.artist import Artist
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

class iSLATPlot:
    """
    iSLATPlot class to handle the iSLAT plotting functionalities.
    This class is responsible for creating the main plot and managing the
    graphical user interface (GUI) components related to plotting.
    """

    def __init__(self, window, input_spectrum_data, molecules_data, xp1, xp2, wave_data, flux_data, config):
        self.window = window
        self.input_spectrum_data = input_spectrum_data
        self.molecules_data = molecules_data
        self.xp1 = xp1
        self.xp2 = xp2
        self.wave_data = wave_data
        self.flux_data = flux_data
        self.config = config

        # Initialize the plot
        self.canvas = self.create_plot()

    def create_plot(self):
        """
        Create the main plot for iSLAT.
        This method sets up the figure, axes, and various components of the plot.
        Returns a plot canvas that can be embedded into other windows.
        """
        # Extract user settings
        user_settings = self.config['user_settings']
        foreground = user_settings["theme"]["foreground"]
        background = user_settings["theme"]["background"]
        iSLAT_version = self.config['iSLAT_version']

        # Prepare data for plotting
        input_spectrum_data = self.input_spectrum_data
        xp1 = self.xp1
        xp2 = self.xp2
        wave_data = self.wave_data
        flux_data = self.flux_data

        # Create the figure and axes
        fig = Figure(figsize=(15, 8.5))
        gs = GridSpec(nrows=2, ncols=2, width_ratios=[1, 1], height_ratios=[1, 1.5], figure=fig)
        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[1, 1])

        # Configure axes
        ax1.set_ylabel('Flux density (Jy)')
        ax2.set_xlabel('Wavelength (μm)')
        ax2.set_ylabel('Flux density (Jy)')
        ax1.set_xlim(xmin=xp1, xmax=xp2)
        ax1.set_ylim(ymin=min(flux_data), ymax=max(flux_data) + (max(flux_data) / 8))
        ax2.set_title('Line inspection plot', fontsize='medium')

        # Plot data
        data_line, = ax1.plot(wave_data, flux_data, color=foreground, linewidth=1, label='Data')

        # Plot molecule data
        for mol_name, mol_filepath, mol_label in self.molecules_data:
            molecule_name_lower = mol_name.lower()
            if molecule_name_lower == 'h2o':
                ax1.plot([], [], alpha=0.8, linewidth=1, label=mol_label)
            else:
                ax1.plot([], [], alpha=0.8, linewidth=1, label=mol_label)

        # Add legend
        ax1.legend()

        # Adjust layout
        fig.subplots_adjust(left=0.06, right=0.97, top=0.97, bottom=0.09)

        # Create a canvas for embedding
        canvas = FigureCanvasTkAgg(fig, master=self.window)
        canvas.draw()

        # Return the canvas
        return canvas

    def populate_population_diagram(self):
        """
        Populate the population diagram graph.
        This method handles the creation and configuration of the population diagram.
        """
        # Placeholder for population diagram logic
        pass

    def create_visibility_buttons(self):
        """
        Create visibility buttons for molecules.
        This method handles the creation of buttons for toggling molecule visibility.
        """
        num_rows = 9
        row_height = 0.035
        row_width = 0.19
        total_height = row_height * num_rows
        start_y = 0.52 + (0.45 - total_height) / 2  # Center vertically

        column_labels = ['Molecule', 'Temp.', 'Radius', 'Col. Dens', 'On', 'Del.', 'Color']
        vis_buttons_dict = {}

        # Placeholder for visibility button creation logic
        pass