from iSLAT.COMPONENTS.Molecule import Molecule

class MoleculeDict(dict):
    """A dictionary to store Molecule objects with their names as keys, and to preform operations on the collection of molecules."""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fluxes = {}

    def add_molecule(self, mol_entry, intrinsic_line_width, wavelength_range, model_pixel_res, model_line_width, dist):
        """Add a new molecule to the dictionary using molecule entry data."""
        mol_name = mol_entry["name"]

        # Create a Molecule instance
        molecule = Molecule(
            name=mol_name,
            filepath=mol_entry["file"],
            displaylabel=mol_entry["label"],
            initial_molecule_parameters=self.initial_molecule_parameters.get(mol_name, {}),
            wavelength_range = wavelength_range,
            intrinsic_line_width=self.save_file_data[mol_name]["Broad"],
            model_pixel_res=model_pixel_res,
            model_line_width=model_line_width,
            
            distance = self.save_file_data.get(mol_name, {}).get("Dist", dist),
            radius = self.save_file_data.get(mol_name, {}).get("Rad", None),
            temp = self.save_file_data.get(mol_name, {}).get("Temp", None),
            n_mol = self.save_file_data.get(mol_name, {}).get("N_Mol", None),
            is_active=self.save_file_data.get(mol_name, {}).get("Vis", True)
        )

        # Store the molecule in the dictionary
        self[mol_name] = molecule

        print(f"Molecule Initialized: {mol_name}")
        return molecule

    def load_molecules_data(self, molecules_data, initial_molecule_parameters, save_file_data, wavelength_range, intrinsic_line_width, model_pixel_res, model_line_width, dist):
        """Load multiple molecules data into the dictionary."""
        self.initial_molecule_parameters = initial_molecule_parameters
        self.save_file_data = save_file_data
        for mol_entry in molecules_data:
            self.add_molecule(
                mol_entry,
                intrinsic_line_width=intrinsic_line_width,
                wavelength_range=wavelength_range,
                model_pixel_res=model_pixel_res,
                model_line_width=model_line_width,
                dist=dist
            )