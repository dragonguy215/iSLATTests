from iSLAT.ir_model.spectrum import Spectrum
from iSLAT.ir_model.moldata import MolData
from iSLAT.ir_model.intensity import Intensity
from iSLAT.ir_model.constants import constants as c
import numpy as np

class Molecule:
    def __init__(self, name, intrinsic_line_width, model_pixel_res, model_line_width, distance, wavelength_range, hitran_data, is_visible = True, color= None, displaylabel= None, filepath= None, initial_molecule_parameters = None, temp = None, n_mol= None, radius= None):
        """Initialize a molecule with its parameters.

        Parameters
        ----------
        name: str
            Name of the molecule
        temp: float
            Kinetic temperature in Kelvin
        n_mol: float
            Number density of the molecule in cm^-3
        radius: float
            Radius
        area: float
            Area scaling factor for the intensity component
        """
        self.name = name
        #self.displaylabel = displaylabel
        self.displaylabel = displaylabel if displaylabel is not None else name  # Display label for the molecule
        self.filepath = filepath
        self.color = color

        self.mol_data = MolData(name, filepath)  # Load molecule data from file

        self.hitran_data = hitran_data

        #print("Initial parameters for molecule:", initial_molecule_parameters)

        self.t_kin = initial_molecule_parameters.get('t_kin', temp)  # Kinetic temperature in Kelvin
        self.scale_exponent = initial_molecule_parameters.get('scale_exponent', 1.0)  # Scaling exponent for the intensity
        self.scale_number = initial_molecule_parameters.get('scale_number', 1.0)  # Scaling number for the intensity
        self.radius_init = initial_molecule_parameters.get('radius_init', radius)  # Initial radius

        self.temp = temp if temp is not None else self.t_kin  # Kinetic temperature in Kelvin
        self.radius = radius if radius is not None else self.radius_init  # Radius of the molecule
        self.n_mol_init = float(self.scale_number * (10**self.scale_exponent))
        self.n_mol = n_mol if n_mol is not None else self.n_mol_init
        self.is_visible = is_visible

        self.intrinsic_line_width = intrinsic_line_width
        self.model_pixel_res = model_pixel_res  # Pixel resolution for the model spectrum
        self.model_line_width = model_line_width  # Line width for the model spectrum
        self.distance = distance
        self.wavelength_range = wavelength_range #if wavelength_range is not ((None, None) or None) else (0.3, 1000)
        
        self.intensity = Intensity(self.mol_data)
        self.calculate_intensity()

        self.spectrum = Spectrum(
            lam_min = self.wavelength_range[0],
            lam_max = self.wavelength_range[1],
            dlambda = self.model_pixel_res,
            R = self.model_line_width, 
            distance = self.distance
        )

        self.spectrum.add_intensity(
            intensity=self.intensity,
            dA=self.radius * 2 ** np.pi
        )

    def calculate_intensity(self):
        print(f"Calculating intensity for {self.name} with T={self.t_kin}, n_mol={self.n_mol}, dv={self.intrinsic_line_width}")
        self.intensity.calc_intensity(
            t_kin=self.t_kin,
            n_mol=self.n_mol,
            dv=self.intrinsic_line_width
        )
    
    def get_flux(self, wavelength_array):
        lam_grid = self.spectrum._lamgrid
        flux_grid = self.spectrum.flux
        return np.interp(wavelength_array, lam_grid, flux_grid)