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
        self.create_plot()
    
    def create_plot(self):
        """
        Create the main plot for iSLAT.
        This method sets up the figure, axes, and various components of the plot.
        """
        # Extract user settings
        user_settings = self.config['user_settings']
        foreground = user_settings["theme"]["foreground"]
        background = user_settings["theme"]["background"]
        iSLAT_version = self.config['iSLAT_version']

        # Set the matplotlib style
        plt.style.use(user_settings["theme"]["matplotlib_style"])

        # Prepare data for plotting
        input_spectrum_data = self.input_spectrum_data
        xp1 = self.xp1
        xp2 = self.xp2
        wave_data = self.wave_data
        flux_data = self.flux_data

'''# Creating the graph
fig = plt.figure(figsize=(15, 8.5))
# fig = plt.figure()
gs = GridSpec (nrows=2, ncols=2, width_ratios=[1, 1], height_ratios=[1, 1.5])
ax1 = fig.add_subplot (gs[0, :])
ax2 = fig.add_subplot (gs[1, 0])
ax3 = fig.add_subplot (gs[1, 1])

# Create a text box for streaming useful information in the third section
# text_box = fig.add_subplot(gs[1, 0])  # This is the third section
# text_box_data = TextBox(text_box, label='', color = background, hovercolor= background)
# ax3.set_ylabel(r'ln(4πF/(hν$A_{u}$$g_{u}$))')
# ax3.set_xlabel(r'$E_{u}$')
ax2.set_xlabel ('Wavelength (μm)')
ax1.set_ylabel ('Flux density (Jy)')
ax2.set_ylabel ('Flux density (Jy)')
ax1.set_xlim (xmin=xp1, xmax=xp2)
plt.rcParams['font.size'] = 10
data_line, = ax1.plot (wave_data, flux_data, color=foreground, linewidth=1)

for mol_name, mol_filepath, mol_label in molecules_data:
    molecule_name_lower = mol_name.lower ()

    if molecule_name_lower == 'h2o':
        exec (
            f"{molecule_name_lower}_line, = ax1.plot({molecule_name_lower}_spectrum.lamgrid, fluxes_{molecule_name_lower}, alpha=0.8, linewidth=1)",
            globals ())
    else:
        exec (f"{molecule_name_lower}_line, = ax1.plot([], [], alpha=0.8, linewidth=1)", globals ())
    exec (f"{molecule_name_lower}_line.set_label('{mol_label}')", globals ())
data_line.set_label ('Data')
sum_line, = ax1.plot ([], [], color='purple', linewidth=1)
sum_line.set_label ('Sum')
ax1.legend ()

ax2.set_frame_on (False)
ax3.set_frame_on (False)

# make empty lines for the second plot
ax2.set_title ('Line inspection plot', fontsize='medium')
data_line_select, = ax2.plot ([], [], color=foreground, linewidth=1)

# Scaling the y-axis based on tallest peak of data
range_flux_cnts = input_spectrum_data[(input_spectrum_data['wave'] > xp1) & (input_spectrum_data['wave'] < xp2)]
range_flux_cnts.index = range (len (range_flux_cnts.index))
fig_height = np.nanmax (range_flux_cnts.flux)
fig_bottom_height = np.min (range_flux_cnts.flux)
ax1.set_ylim (ymin=fig_bottom_height, ymax=fig_height + (fig_height / 8))

# adjust the plots to make room for the widgets
fig.subplots_adjust (left=0.06, right=0.97, top=0.97, bottom=0.09)

# Populating the population diagram graph
pop_diagram ()

num_rows = 9

# Calculate the height and width of each row
row_height = 0.035
row_width = 0.19

# Calculate the total height of all rows
total_height = row_height * num_rows

# Calculate the starting y-position for the first row within the control_border
start_y = 0.52 + (0.45 - total_height) / 2  # Center vertically

# Define the column labels
column_labels = ['Molecule', 'Temp.', 'Radius', 'Col. Dens', 'On', 'Del.', 'Color']

# Create a dictionary to store the visibility buttons
vis_buttons_dict = {}

# Create a tkinter window
window = tk.Tk ()
window.title ("iSLAT " + iSLAT_version)'''