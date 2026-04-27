""" Module to create sample lists for the OT2
Author: Kaarthik, Takumi

"""
import csv
import datetime
import inspect
import os
import datetime
import copy
import shutil
import sys
import math
from io import StringIO
import platform
import pathlib
from pathlib import Path
import pandas as pd
import pickle
import numpy as np
import re
import ast
import itertools
import rotbot_control # Otto Control module

# Necessary for Path library to work on Windows and Linux
plt = platform.system()
if plt == 'Linux': pathlib.WindowsPath = pathlib.PosixPath

# FILE ARCHITECTURE
master_dir = Path(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sample_dir = master_dir / 'samples'
out_dir = master_dir / 'archive'

# AVAIlABLE LAB WARE
plate_types = ['biorad_96_wellplate_200ul_pcr', 'thermoscientificnunc_96_wellplate_2000ul']
pipette_types = ['p20_single_gen2', 'p300_single_gen2', 'p1000_single_gen2']
tiprack_types = ['opentrons_96_tiprack_20ul', 'opentrons_96_tiprack_300ul', 'opentrons_96_tiprack_1000ul']

# Well locations for the 96-well plate- [A1, A2, ..., A12, B1, ..., B12]
_wells = itertools.cycle([f"{letter}{number}" for letter in 'ABCDEFGH' for number in range(1, 13)]) # 96-well plate

class Sample:
    """
    A class to represent a sample that needs to be prepared by the OT2.
    """

    def __init__(self, composition, order_mixing=[], stock_dict={}, volume=None, solvent='water',
                 viscous_species=['peg'], special_species=['bsa', 'fus', 'bik1'],
                 skip_species=[], name="", id=1, metadata="",
                 savedir=sample_dir):
        # ATTRIBUTES
        ## BASIC ATTRIBUTES
        self.name = name # Name of the sample
        self.id = id # ID of the sample
        self.solvent = solvent # Solvent used in the sample, default is water
        self.composition = composition  # dictionary containing the composition of the sample
        self.order_mixing = order_mixing # list containing the order of mixing the chemicals
        ### Sanity checks
        self._check_composition()  # Check if the solvent is in the composition
        self._check_order_mixing()  # Check if the solvent is in the composition
        # Boolean list which decides if the robot needs to mix thesample after dispensing each chemical
        self.mixAfterDispence = {chemical: False for chemical in self.composition.keys()}
        self.mixAfterDispence[self.solvent] = True # Always mix when the solvent is added

        ## Dictionary containing the relative viscosity of the chemicals
        self.relative_viscosity = self._set_relative_viscosity(viscous_species=viscous_species, 
                                                               special_species=special_species)
        self.chemicals = list(composition.keys())  # List of chemicals in the sample
        self.volume = volume  # Volume of the sample in uL
        # Stock solutions (mM or g/L), dict(chemical: concentration)
        self.stock_dict = stock_dict  # e.g. {'bsa': 46.7, 'peg': 266.0, 'kp': 100.0, 'kcl': 200.0, 'nh3': 0.0, 'water': 0}
        self.recipe = {}  # Recipe for the sample- {chemical1: volume in uL, chemical2:, volume in uL, ...}
        self.skip_species = skip_species # List of species to skip during the mixing step

        ## METADATA
        self.prepared = 0 # 0: not prepared, 1: prepared
        self.metadata = metadata # Metadata of the sample
        self.timestamp = '1900-01-01 12:00:00' # Default directory to save the samples
        self.filepath = ''

        # If the stock_dict and volume are given, update the recipe
        if self.stock_dict != {}:
            self._check_stock_dict()
            if self.volume is not None:
                self._update_recipe()

        # LABWARE
        self.plate_deck_loc = 10 # Deck number of the well plate
        self.well = "" # Well location of the sample. e.g. A1, B2, ...
        self.subwell = None  # dict (chemical: well)- used when the sample prep requires multiple wells
        self.well_obj = None # opentrons.protocol_api.labware.Well object will be stored once the sample is prepared
        self.current_volume = 0
        self.stock_deck_loc = 6 # Deck number of the stock tube rack


    # METHODS
    def _check_composition(self):
        """
        This function checks if the composition contains the solvent.
        Returns
        -------

        """
        if self.solvent not in self.composition.keys():
            self.composition[self.solvent] = 0
            print(f'{self.solvent} is not in the composition. Added {self.solvent} to the composition.')

    def _check_order_mixing(self, ):
        """
        This function checks if the order of mixing contains all the chemicals in the composition.
        ... If the order of mixing is not given, it sets the order of mixing as the sorted list of chemicals.

        Returns
        -------

        """
        # If order of mixing is not given, set the order of mixing as the sorted list of chemicals
        if self.order_mixing == []:
            self.order_mixing = sorted(list(self.composition.keys()))

        # If solvent is not in the order of mixing, add it to the beginning
        if self.solvent not in flatten_list(self.order_mixing):
            self.order_mixing = [self.solvent] + self.order_mixing # Add solvent first by default
            print(f'{self.solvent} is not in the order of mixing. Added {self.solvent}'
                  f' to first place of the order of mixing.')
        # Allow only the chemicals in the composition
        if len(flatten_list(self.order_mixing)) > len(self.composition):
            new_order = self.order_mixing
            for chem in self.composition.keys():
                new_order = str(new_order).replace("'" + chem + "'" + ',', '')
                new_order = str(new_order).replace("'" + chem + "'", '')
            ast.literal_eval(new_order)
            self.order_mixing = new_order
        if len(flatten_list(self.order_mixing)) != len(self.composition):
            print('Sample: order_mixing:', self.order_mixing, 'composition:', self.composition.keys())
            raise ValueError('The order of mixing must contain all the chemicals in the composition.')
    def _check_stock_dict(self):
        if self.solvent not in self.stock_dict.keys():
            self.stock_dict[self.solvent] = 0
            print(f'{self.solvent} is not in the stock_dict. Added {self.solvent} to the stock_dict.')


    def _update_composition(self, composition, **kwargs):
        """
        This function updates the composition of the sample.
        ... It also updates the relative viscosity based on the new composition.

        Parameters
        ----------
        composition: dict, dictionary containing the composition of the sample
        kwargs: dict, keyword arguments to update the relative viscosity

        Returns
        -------

        """
        self.composition = composition
        self._set_relative_viscosity(**kwargs)  # Update relative viscosity based on new composition


    def _set_relative_viscosity(self, viscous_species=['peg'], special_species=['bsa']):
        relative_viscosity = get_relative_viscosity(self.composition.keys(), viscous=viscous_species, special=special_species)
        self.relative_viscosity = relative_viscosity
        return relative_viscosity


    def _update_recipe(self):
        self.recipe = get_recipe(self.volume, self, self.stock_dict, self.solvent)

    def _update_metadata(self, metadata):
        self.metadata = metadata

    def _update_name(self, name):
        self.name = name

    def toDataFrame(self):
        """
        Converts attributes of the class instance into a pandas DataFrame.

        Returns:
            DataFrame: A pandas DataFrame with attributes as columns and their values as rows.
        """
        # Get dictionary of instance attributes
        attributes = self.__dict__

        # Convert to DataFrame
        df = pd.DataFrame([attributes.values()], columns=attributes.keys())
        return df


    def dilute(self, dilution_factor, verbose=False):
        """
        This function dilutes the sample by a given dilution factor.
        ... The volume of the sample is reduced by the dilution factor.

        Parameters
        ----------
        dilution_factor: float, dilution factor by which the sample needs to be diluted

        Returns
        -------
        None
        """
        self.volume = self.volume * dilution_factor
        new_composition = {species: conc / dilution_factor for species, conc in self.composition.items() if species != self.solvent}
        self.composition = new_composition
        self._update_recipe()
        if verbose:
            print(f'Sample diluted by a factor of {dilution_factor}. New volume: {self.volume} uL')

    def copy(self):
        """
        This function creates a copy of the Sample object.

        Returns
        -------
        sample: Sample object, copy of the Sample object
        """
        return copy.deepcopy(self)

    def save(self, savedir, filename=None,verbose=True):
        """
        This function saves the Sample object as a pickle file.
        ... The filepath is 'savedir/sample_<id>_<name>.pkl'
        ... This also updates the filepath attribute of the Sample object.

        Parameters
        ----------
        savedir: str or Path object, directory to save the Sample object
        verbose: bool, whether to print the filepath or not

        Returns
        -------
        filepath: Path object, path of the saved file
        """
        if isinstance(savedir, str):
            savedir = Path(savedir)
        if not os.path.exists(savedir):
            os.makedirs(savedir)

        if filename is None:
            # Default filename
            filename = f'sample_{self.id:03d}_{self.name}.pkl'

        filepath = savedir / filename
        with open(filepath, 'wb') as f:
            pickle.dump(self, f)
        if verbose:
            print(f'Sample saved as {filepath}')
        self.filepath = filepath
        return filepath

# STEP 1: DEFINE STOCK CONCENTRATIONS
# STOCK SOLUTIONS
def write_stocks(stock_dict, savedir='./preparation', filename='stocks.csv'):
    """
    Write chemical information of stock solutions to a CSV file.

    Parameters:
    stock_dict: dict, dictionary containing the stock concentrations of the chemicals
        ... e.g. {'fus_buffer': 50., 'NaHCO3': 25., 'water': 0.} # given in mM or g/L
    savedir: str or Path object, directory to save the CSV file
    filename: str, name of the CSV file, default is 'stocks.csv'

    Returns:
    filepath: Path object, path of the saved file
    """
    # Ensure savedir is a Path object
    if not isinstance(savedir, Path):
        savedir = Path(savedir)

    # Create a directory if it does not exist
    savedir.mkdir(parents=True, exist_ok=True)

    # Define the filepath
    filepath = savedir / filename

    # Open the CSV file in write mode
    with open(filepath, 'w', newline='') as csvfile:
        # Define the fieldnames (column names) for the CSV
        fieldnames = ['Chemical', 'Concentration', 'Comments']

        # Create a CSV writer object
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write the header row
        writer.writeheader()

        # Write each chemical and its concentration as a row in the CSV
        for i, (chemical, concentration) in enumerate(stock_dict.items()):
            writer.writerow({'Chemical': chemical, 'Concentration': concentration, 'Comments': ''})

    print(f"write_stocks:\n...Stock solution file was created at {filepath}.")
    return filepath

def read_stocks(filepath='./preparation/stocks.csv'):
    """
    Read the stock solutions from a CSV file.
    Parameters
    ----------
    filepath: str or Path object, path to the CSV file

    Returns
    -------
    stock_dict: dict, dictionary containing the stock concentrations of the chemicals
        ... e.g. {'fus_buffer': 50., 'NaHCO3': 25., 'water': 0.} # given in mM or g/L
    """
    # Ensure filepath is a Path object
    if not isinstance(filepath, Path):
        filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found at {filepath}. "
                                f"Define the stock solutions like "
                                "{'fus_buffer': 50., 'NaHCO3': 25., 'water': 0}.")

    # Open the CSV file in read mode
    with open(filepath, 'r', newline='') as csvfile:
        stock_dict = {}
        # Create a CSV reader object
        reader = csv.DictReader(csvfile)
        stock_dict = {row['Chemical']: float(row['Concentration']) for row in reader}
    return stock_dict

# STEP 2: CREATE SAMPLE OBJECT(S)
## Single sample
def create_sample(composition={'fus_buffer': 50., 'proline': 25, 'water': 0},
                  order_mixing=['proline', 'fus_buffer', 'water'],
                  solvent='water',
                  stock_dict={}, volume=None, # If given, the recipe for the samples will be updated!
                  id=1,
                  name='', metadata='',
                  save=False, dst_dir='./samples', # If True, save the samples to a pickle file
                  **kwargs
                  ):
    """
    This function creates a sample with the given composition and order of mixing.
    ... If 'save' is True, the Sample object is saved as a pickle file.

    Parameters
    ----------
    composition: dict, dictionary containing the composition of the sample
        ... e.g. {'fus_buffer': 50., 'NaHCO3': 25., 'water': 0.} # given in mM or g/L
    order_mixing: list, list containing the order of mixing the chemicals
        ... e.g. ['NaHCO3', 'fus_buffer', 'water']
    stock_dict: dict, dictionary containing the stock concentrations of the chemicals (optional but recommended)
        ... e.g. {'fus_buffer': 50., 'NaHCO3': 25., 'water': 0.} # given in mM or g/L
        ... e.g. {'bsa': 46.7, 'peg': 266.0, 'kp': 100.0, 'kcl': 200.0, 'nh3': 0.0, 'water': 0}
        ... If given with volume, the recipe for the samples will be updated!
    volume: float, volume of the sample in uL (optional but recommended)
        ... e.g. 118 # uL
    id: int, ID of the sample
        ... If given, the id of the sample will be updated! default: 1
    name: str, name of the sample, default is ''
        ... Adds name to the Sample object. sample.name = name
        ... Used in the filename if save is True
    metadata: str, metadata of the sample
        ... Adds metadata to the Sample object. sample.metadata = metadata
    save: bool, whether to save the samples to a pickle file or not
    dst_dir: str, directory to save the samples, default is './samples'
        ... If save is True, the samples will be saved in the directory 'dst_dir'
    Returns
    -------

    """
    comp_ = copy.deepcopy(composition)
    comp_[solvent] = 0  # Always set the solvent conc 0
    sample = Sample(comp_, order_mixing, solvent=solvent,
                    stock_dict=stock_dict, volume=volume,
                    id=id, name=name, metadata=metadata, **kwargs)
    if save:
        sample.save(dst_dir)
    return sample
        
# Multiple samples (variation of a chemical)
def create_samples(composition={'fus_buffer': 50., 'proline': 0, 'water': 0},
                   vary='proline', concentrations=np.linspace(0, 500, 6), solvent='water',
                   order_mixing=['proline', 'fus_buffer', 'water'],
                   stock_dict={}, volume=None,
                   save=False, # If True, save the samples to a pickle file
                   name='', metadata='',
                   repeat=1,
                   **kwargs,):
    """
    This function creates a list of samples with variation of a given chemical.
    ... The user has to input various details of the samples such as the chemicals, their composition, etc.
    ... The sample list is saved as a pickle file.

    Parameters
    ----------
    composition: dict, dictionary containing the composition of the sample
        ... e.g. {'fus_buffer': 50., 'NaHCO3': 25., 'water': 0.} # given in mM or g/L
    vary: str, chemical that needs to be varied
        ... e.g. 'proline'
    concentrations: list, list of concentrations of the chemical that needs to be varied
        ... e.g. np.linspace(0, 500, 6) for 0, 100, 200, 300, 400, 500
        ... e.g. [0, 25, 50, 100, 200, 300, 400, 500, 600, 700]
    solvent: str, solvent used in the samples, default is 'water'
    order_mixing: list, list containing the order of mixing the chemicals
        ... e.g. ['NaHCO3', 'fus_buffer', 'water']
    stock_dict: dict, dictionary containing the stock concentrations of the chemicals (optional but recommended)
        ... e.g. {'fus_buffer': 50., 'NaHCO3': 25., 'water': 0.} # given in mM or g/L
        ... e.g. {'bsa': 46.7, 'peg': 266.0, 'kp': 100.0, 'kcl': 200.0, 'nh3': 0.0, 'water': 0}
        ... If given with volume, the recipe for the samples will be updated!
    volume: float, volume of the sample in uL (optional but recommended)
        ... e.g. 118 # uL
        ... If given with stock_dict, the recipe for the samples will be updated!
    name: str, name of the sample list
        ... Adds name to the Sample object(s). sample.name = name
        ... If save is True, the samples will be saved as 'sample_list-<name>-<timestamp>.pkl'
    metadata: str, metadata of the samples
        ... Adds metadata to the Sample object(s). sample.metadata = metadata
    save: bool, whether to save the samples to a pickle file or not
    sample_dir: str, directory to save the samples
        ... If save is True, the samples will be saved in the directory 'sample_dir'

    Returns
    -------
    samples: list, list of Sample objects

    Examples
    --------
    samples = otto.create_samples(composition = {'fus_buffer': 50., NaHCO3: 0, 'water': 0},
                             vary=NaHCO3, concentrations=[0, 25, 50, 100, 200, 300, 400, 500, 600, 700],
                             solvent='water',
                             order_mixing=[NaHCO3, 'fus_buffer', 'water'],
                             stock_dict=stock_dict, volume=118, # Provide with stock_dict and sample volume to create a recipe
                             save=True,
                            )
    """
    if set(composition.keys()) != set(flatten_list(order_mixing)):
        print(list(set(flatten_list(order_mixing))))
        raise ValueError(f'order_mixing must contain {composition.keys()}')

    samples = []
    count = 1
    for i, conc in enumerate(concentrations):
        for n in range(repeat):
            if repeat > 1:
                footer = f'_n{n + 1:03d}'
            else:
                footer = ''
            comp_ = copy.deepcopy(composition)
            comp_[vary] = conc
            comp_[solvent] = 0  # Always set the solvent conc 0
            sample = Sample(comp_, order_mixing, stock_dict=stock_dict, volume=volume,
                            solvent=solvent,
                            id=count, name=name + f'{i+1:03d}{footer}', metadata=metadata, **kwargs,
                            )
            samples.append(sample)
            count += 1

    if save:
        save_samples(samples, name=name + '_variation', dst_dir=sample_dir)
    return samples


def create_samples_from_csv(csv_file, stock_csv_file=r'./samples/csv/stocks.csv', savedir=r'./samples'):
    """
    This function creates samples from a CSV file.
    ... The CSV file must contain the composition of the samples, the volume, and the order of mixing.

    Parameters
    ----------
    csv_file: str or Path object, path to the CSV file containing the samples
    stock_csv_file: str or Path object, path to the CSV file containing the stock solutions
    savedir: str or Path object, directory to save the samples

    Returns
    -------
    samples: list, list of Sample objects
    """
    def format_order_mixing(s):
        ss = ""
        s = s.replace(' ', "")
        for elem in s.split(','):
            if elem[0] == '[':
                ss += '[' + f"'{elem.strip()[1:]}', "

            elif elem[-1] == ']':
                ss += f"'{elem.strip()[:-1]}'" + '], '
            else:
                ss += f"'{elem.strip()}', "

        ss = ss[:-2]
        s_as_list = list(ast.literal_eval(ss))
        return s_as_list

    df = pd.read_csv(csv_file)
    n = len(df)
    samples = []
    for i in range(n):
        chemicals = [name for name in df.columns if name.lower() not in ['name', 'volume', 'order of mixing', 'comments']]
        composition = df[chemicals].iloc[i].to_dict()
        order = format_order_mixing(df.iloc[i]['Order of mixing'])
        volume = df.iloc[i]['volume']
        name = df.iloc[i]['Name']

        stock_dict = read_stocks(filepath=stock_csv_file)
        s = create_sample(composition=composition,
                            order_mixing=order,
                            solvent='water',
                            stock_dict=stock_dict,
                            volume=volume,
                            id=1,
                            name=name,
                            metadata='',
                            save=False,
                            dst_dir=savedir,
                            repeat=1)
        samples.append(s)
    save_samples(samples, filepath=Path(savedir) / Path(csv_file).stem)

# Saving Sample objects
def save_sample(sample, name='', metadata='', dst_dir=sample_dir, verbose=True):
    """
    This function writes the sample to a pickle file.
    ... The sample is saved as a Sample object.

    Parameters
    ----------
    sample: Sample object
    name: str, name of the sample
    metadata: str, metadata of the sample
    dst_dir: str or Path object, directory to save the sample
    verbose: bool, whether to print the filepath or not

    Returns
    -------
    filepath: Path object, path of the saved file
    """
    if name !='':
        print(f'Updating name of the sample to {name}')
        sample.name = name
    if metadata != '':
        print(f'Updating metadata of the sample')
        sample.metadata = metadata
    filepath = sample.save(dst_dir, verbose=verbose)
    return filepath

def save_samples(samples, filepath=None, name='', metadata='',
                 dst_dir=sample_dir, verbose=True):
    """
    This function writes the samples to a pickle file.
    ... The samples are saved as a list of Sample objects.
    ... The name of the file is 'sample_list-<name>-<timestamp>.pkl'

    Parameters
    ----------
    samples: list or Sample object, list of Sample objects
    name: str, name of the sample list
        ... Used in the filename
    metadata: str, metadata of the samples
        ... Will updata the comments of the samples. sample.metadata = metadata

    Returns
    -------
    filepath: Path object, path of the saved file
    """

    # Name of the file
    ## Default filepath: '/samples/sample_list-<name>.pkl'
    if filepath is None:
        filename = 'sample_list-' + name + '.pkl' # Name of the file
        filepath = dst_dir / filename # Path to save the samples
    else:
        if not str(filepath).endswith('.pkl'):
            filepath = Path(str(filepath) + '.pkl')
        filename = Path(filepath).name

    # If only one sample is given, convert it to a list
    if isinstance(samples, Sample):
        samples = [samples]
    if isinstance(dst_dir, str):
        dst_dir = Path(dst_dir)

    for i, sample in enumerate(samples):
        if name != '':
            sample.name = f'{name}_{i:03d}'
        sample.timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        sample.metadata = metadata
        sample.filepath = filepath

    # Save the samples to a pickle file
    dst_dir = Path(dst_dir)
    if not os.path.exists(dst_dir):
        os.makedirs('samples')
    with open(filepath, 'wb') as f:
        # pickle.dump([metadata, sample.relative_viscosity, samples], f) # Old version
        pickle.dump(samples, f)
    if verbose:
        print(f'save_samples:\n...Sample list was saved as {filename}')
    return filepath

def read_samples(sampleDir=sample_dir, flatten=True, verbose=True, sortByWell=True):
    """
    This function reads all samples in the given directory.

    Parameters
    ----------
    sampleDir: str or Path object, directory containing the samples
    verbose: bool, whether to print the filepath or not

    Returns
    -------
    all_samples: list, list of Sample objects
    """

    def alphanum_key(sample):
        """ Turn a string into a list of string and number chunks.
            "A10" -> ("A", 10) """
        if isinstance(sample.well, list):
            return tuple(int(text) if text.isdigit() else text for text in re.split('([0-9]+)', sample.well[0]))
        else:
            return tuple(int(text) if text.isdigit() else text for text in re.split('([0-9]+)', sample.well))


    if not isinstance(sampleDir, Path):
        sampleDir = Path(sampleDir)
    filepaths = list(sampleDir.glob('*.pkl'))
    print(f'Reading samples from {sampleDir}... {len(filepaths)} files found')
    all_samples = []
    for filepath in filepaths:
        with open(filepath, 'rb') as f:
            samples = pickle.load(f)
        all_samples.append(samples)
        if verbose:
            print(f'... Returning {filepath.name}')
    if flatten:
        all_samples = [item for sublist in all_samples for item in sublist]
        if sortByWell:
            all_samples = sorted(all_samples, key=alphanum_key)
        # Finally, sort samples by plate_deck_loc
        all_samples = sorted(all_samples, key=lambda x: x.plate_deck_loc)
    else:
        if sortByWell:
            # First sort samples in each list by well
            for i, samples in enumerate(all_samples):
                all_samples[i] = sorted(samples, key=alphanum_key)
            # Then sort the lists by the well of the first sample in each list
            all_samples = sorted(all_samples, key=lambda x: alphanum_key(x[0]))
            # Finally, sort samples by plate_deck_loc
            all_samples = sorted(all_samples, key=lambda x: x.plate_deck_loc)
    return all_samples



# STEP 3: DEFINE LABWARE CONFIGURATION
## STEP 3-1: ESTIMATE REQUIRED LABWARE AND CHEMICALS

def estimate_requirements(samples, verbose=True, pipetteTypes=[None, None]):
    """
    This function estimates the labware and chemicals required for the run based on the sample list and stock concentrations.

    Parameters
    ----------
    samples: list, list of Sample objects
        ... e.g. [sample1, sample2, ...]
    verbose: bool, whether to print the estimated requirements or not, default is True

    Returns
    -------
    nPlates: int, number of plates required for the run
    pipetteTypes: list, list of pipette type (p20, p300, p1000) required for the run
    nTipRacks: list, list of tip boxes required for the run
    nTips: list, list of tips required for the run
    requiredStockVol: dict, dictionary containing the stock volumes required for the run
    """
    def isStockConsistent(samples):
        """
        This function checks if the given samples refer to the same stock solutions.

        Parameters
        ----------
        samples: list, list of Sample objects

        Returns
        -------
        bool: True if the samples refer to the same stock solutions, False otherwise
        """
        stock_dict = samples[0].stock_dict
        for sample in samples:
            if sample.stock_dict != stock_dict:
                return False
        return True

    def doesRecipeExist(samples):
        """
        This function checks if the recipe exists for the samples.

        Parameters
        ----------
        samples: list, list of Sample objects

        Returns
        -------
        bool: True if the recipe exists, False otherwise
        """
        recipeExists = all(sample.recipe!={} for sample in samples)
        return recipeExists
    def estimate_required_stock_volumes(samples, units='mL'):
        """
        This function estimates the stock volumes required for the run based on the sample list and stock concentrations.
        The stock volumes are calculated based on the recipe for the samples.

        Parameters
        ----------
        samples: list, list of Sample objects
        units: str, units of the stock volumes, default is 'mL'

        Returns
        -------
        required_stock_volumes: dict, dictionary containing the stock volumes (mL) required for the run
        """
        required_stock_volumes = dict.fromkeys(samples[0].stock_dict, 0)
        for sample in samples:
            for chemical, volume in sample.recipe.items():
                vol = sample.recipe[chemical]
                if units == 'mL':
                    vol = vol / 1000
                elif units == 'uL':
                    pass
                required_stock_volumes[chemical] += vol
        return required_stock_volumes

    def estimate_pipettes_and_number_of_tips(samples, pipetteTypes=[None, None]):
        """
        This function estimates the number of tips required for the run based on the pipetting volumes.

        Parameters
        ----------
        pipetting_volumes: list, list of volumes to be dispensed

        Returns
        -------
        tips: list, list of tips required for the run
        """
        # CHECK IF RECIPES EXIST IN THE SAMPLES
        if not doesRecipeExist(samples):
            raise ValueError('Recipe does not exist for the samples',
                             'Define the stock solutions and the solvent for the samples via sample.stock_dict and sample.solvent.',
                             'Then create the recipe using the sample.update_recipe() method.',
                             )
        else:
            # CHECK IF THE SAMPLES REFER TO THE SAME STOCK SOLUTIONS
            if not isStockConsistent(samples):
                print('Error: The samples refer to different stock solutions. '
                      'Please check the stock concentrations via samples[i].stock_dict for i in range(len(samples))')
                raise ValueError('Inconsistent stock solutions')

        # Scan the pipetting volumes to determine the pipettes required
        pipetting_volumes = []
        for sample in samples:
            volume_list = list(sample.recipe.values())
            if any(vol < 0 for vol in volume_list):
                raise ValueError(f'Error: Negative volume for sample {sample.name} \n sample recipe: {sample.recipe}')
            else:
                pipetting_volumes += volume_list

        # Determine if p20 is required
        isP20Required = any(vol < 20 for vol in pipetting_volumes)

        # If p20 is required, check if p20 is provided in pipetteTypes
        if isP20Required:
            if 'p20' not in pipetteTypes and pipetteTypes != [None, None]:
                raise ValueError('Error: p20 is required but not provided in pipetteTypes')

        if isP20Required:
            # Estimate the number of tips required for the run using p20 and p1000
            tips_20_1000 = estimate_required_tips([0, 20], [100, 1000], pipetting_volumes)
            # Estimate the number of tips required for the run using p20 and p300
            tips_20_300 = estimate_required_tips([0, 20], [20, 300], pipetting_volumes)

            nTipsUsingP20P1000 = tips_20_1000[0] + tips_20_1000[1]
            nTipsUsingP20P300 = tips_20_300[0] + tips_20_300[1]
            if nTipsUsingP20P1000 < nTipsUsingP20P300:
                pipettesTypes = ['p20', 'p1000']
                nTipRacks = [np.ceil(tips_20_1000[0] / 96), np.ceil(tips_20_1000[1] / 96)]
                nTips = [tips_20_1000[0], tips_20_1000[1]]
            else:
                pipettes_types = ['p20', 'p300']
                nTipRacks = [np.ceil(tips_20_300[0] / 96), np.ceil(tips_20_300[1] / 96)]
                nTips = [tips_20_300[0], tips_20_300[1]]
        else:
            tips_300_1000 = estimate_required_tips([20, 300], [300, 1000], pipetting_volumes)

            pipettesTypes = ['p300', 'p1000']
            nTipRacks = [np.ceil(tips_300_1000[0] / 96), np.ceil(tips_300_1000[1] / 96)]
            nTips = [tips_300_1000[0], tips_300_1000[1]]

        nTipRacks, nTips = [int(val) for val in nTipRacks], [int(val) for val in nTips]
        return pipettes_types, nTipRacks, nTips

    if isinstance(samples[0], list):
        print('estimate_requirements:\n...You have entered a list of samples (a list of Sample objects).',
              'Flattening the list of lists to a single list.')
        samples = [item for sublist in samples for item in sublist]


    # Number of samples
    nSamples = len(samples)

    # Number of plates required for the run
    nPlates = int(np.ceil(nSamples / 96))

    # Estimate the number of plates required for the run
    pipetteTypes, nTipRacks, nTips = estimate_pipettes_and_number_of_tips(samples, pipetteTypes=pipetteTypes)

    # Estimate the stock volumes required for the run
    requiredStockVol = estimate_required_stock_volumes(samples, units='mL')

    if verbose:
        print('-'*40)
        print(f'Number of Plates Required: {nPlates}')
        print(f'Pipettes Required: {pipetteTypes}')
        print(f'Number of Tip Boxes Required: {nTipRacks}')
        print(f'Number of Tips Required: {nTips}')
        print(f'Stock Volumes Required in mL: {requiredStockVol}')
        print('-'*40)
        
    return nPlates, pipetteTypes, nTipRacks, nTips, requiredStockVol



## DEFAULT LABWARE CONFIGURATION (EXPECTED TO BE EDITED BY THE USER)
def create_default_labware_configuration(samples,
                                         plateType='thermoscientificnunc_96_wellplate_2000ul',
                                         pipetteTypes=[None, None],
                                         filepath='./preparation/labware_configuration.csv',
                                         stock_filepath='./preparation/stocks.csv'):
    # nPlates, pipetteTypes, nTipRacks, nTips, requiredStockVol
    """
    This function writes the labware configuration to a csv file.
    ... Labware names can be found in the Opentrons Labware Library.
    ...... Link: https://labware.opentrons.com/
    ... The default labware configuration is as follows:
        ... Plates: thermoscientificnunc_96_wellplate_2000ul
        ... pipetteTypes: p20, p300, p1000
        ... Tip boxes: opentrons_96_tiprack_20ul, opentrons_96_tiprack_300ul, opentrons_96_tiprack_1000ul
        ... Stock tubes: opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical
    ... The user is expected to edit the csv file to change the labware configuration as required.

    Parameters
    ----------
    nPlates: int, number of plates required for the run
        ... Choose from ['biorad_96_wellplate_200ul_pcr', 'thermoscientificnunc_96_wellplate_2000ul']
    pipetteTypes: list, list of pipetteTypes required for the run
        ... Choose from ['p20', 'p300', 'p1000']
        ... Corresponding names are ['p20_single_gen2', 'p300_single_gen2', 'p1000_single_gen2']
    tipracks: list, list of tip boxes required for the run, default is 'auto'
        ... If 'auto', the tip boxes are chosen based on the pipetteTypes. (Recommended)
        ... For reference, available tickrack names are
            ['opentrons_96_tiprack_20ul', 'opentrons_96_tiprack_300ul', 'opentrons_96_tiprack_1000ul']
    filepath: str, path to save the labware configuration, default is './preparation/labware_configuration.csv'
    stock_filepath: str, path to the stock solutions file, default is './preparation/stocks.csv'

    Returns
    -------

    """
    def isStockConsistent(samples, stock_dict):
        """
        This function checks if the given samples refer to the same stock solutions.

        Parameters
        ----------
        samples: list, list of Sample objects

        Returns
        -------
        bool: True if the samples refer to the same stock solutions, False otherwise
        """
        if isinstance(samples[0], list):
            ans = all([isStockConsistent(samples_, stock_dict) for samples_ in samples])
            return ans
        else:
            for sample in samples:
                required_chemicals = sample.stock_dict.keys()
                chemicals_in_stock = stock_dict.keys()
                if not set(required_chemicals).issubset(chemicals_in_stock):
                    return False
            return True
    def get_required_chemicals(samples):
        """
        This function returns the chemicals required for the run based on the samples.

        Parameters
        ----------
        samples: list, list of Sample objects

        Returns
        -------
        required_chemicals: list, list of chemicals required for the run
        """

        required_chemicals = set()
        # If the samples are in a list of lists- i.e.-[samples1, samples2, ...]
        if isinstance(samples[0], list):
            for samples_ in samples:
                required_chemicals.update(get_required_chemicals(samples_))
            return list(required_chemicals)
        # If the samples are in a single list- samples = [sample1, sample2, ...]
        else:
            for sample in samples:
                required_chemicals.update(sample.recipe.keys())
        return list(required_chemicals)

    pipetteNames = {'p20': 'p20_single_gen2',
                     'p300': 'p300_single_gen2',
                     'p1000': 'p1000_single_gen2'}
    ## Pipette - Tiprack correspondence: {pipette: tiprack}
    tiprackNames = {'p20': 'opentrons_96_tiprack_20ul',
                'p300': 'opentrons_96_tiprack_300ul',
                'p1000': 'opentrons_96_tiprack_1000ul'}
    ## Available pipette positions
    pipettePos = ['left', 'right']

    # Labware Name: opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical
    falcon15Wells = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
    falcon50Wells = ['A3', 'A4', 'B3', 'B4']


    # INITIALIZATION
    dec_loc = 11

    # LOAD STOCKS FROM FILE
    stock_dict = read_stocks(filepath=stock_filepath) # {'chemical species': concentration}
    if not isStockConsistent(samples, stock_dict):
        raise ValueError('Error: The samples refer to inconsistent stock solutions '
                         f'with the stock solutions in {stock_filepath}')

    # ESTIMATE REQUIREMENTS
    estimates = estimate_requirements(samples, verbose=False, pipetteTypes=pipetteTypes)
    nPlates, pipetteTypes, nTipRacks, nTips, requiredStockVol = estimates

    # FIND REQUIRED CHEMICALS
    requiredChemicals = get_required_chemicals(samples)

    with open(filepath, 'w', newline='') as f:
        fwriter = csv.writer(f, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)

        # WELL PLATES
        fwriter.writerow(['type', 'name', 'deck_location', 'starting_well_position'])
        for i in range(nPlates):
            dec_loc -= 1
            fwriter.writerow(['plate', plateType, str(dec_loc), 'A1'])

        fwriter.writerow([])
        fwriter.writerow([])

        # PIPETTES AND CORRESPONDING TIP BOXES
        # PIPETTE AND TIP BOXES (['pipette_list', 'name', 'location (right/left)'])
        for i, pipetteType in enumerate(pipetteTypes):
            fwriter.writerow(['type', 'name', 'deck_location', 'starting_tip'])
            # PIPETTES
            fwriter.writerow(['pipette', pipetteNames[pipetteType], pipettePos[i], 'n/a'])
            # fwriter.writerow(['tip box', 'name', 'deck location', 'starting tip position'])
            # TIP RACKS
            for j in range(nTipRacks[i]):
                dec_loc -= 1
                starting_tip = 'A1' if j == 0 else 'A1(cannot be changed)'
                fwriter.writerow(['tiprack', tiprackNames[pipetteType], str(dec_loc), 'A1'])
            fwriter.writerow([])
            fwriter.writerow([])

        # STOCK TUBES
        fwriter.writerow(
            ['type', 'chemical_name', 'concentration', 'volume_mL', 'labware_name', 'deck_location',
             'rack_location', 'tube_type'])
        dec_loc -= 1
        def get_default_stock_settings(required_volume):
            try:
                if required_volume < 8.0:
                    vol_stock = 10
                    tube = 'falcon15'
                    loc = falcon15Wells.pop(0)
                elif required_volume < 30.0:
                    vol_stock = 30
                    tube = 'falcon50'
                    loc = falcon50Wells.pop(0)
                else:
                    vol_stock = 50
                    tube = 'falcon50'
                    loc = falcon50Wells.pop(0)
            except:
                vol_stock, tube, loc = 'n/a', 'n/a', 'n/a'
            return vol_stock, tube, loc

        for chemical in requiredChemicals:
            required_vol_stock = requiredStockVol[chemical]
            if required_vol_stock > 0:
                vol_stock, tube, loc = get_default_stock_settings(required_vol_stock)
                fwriter.writerow(
                    ['stock', chemical, stock_dict[chemical], vol_stock, 'opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical',
                     dec_loc, loc, tube])
    print(f'create_default_labware_configuration:\n...Labware configuration saved as {filepath}')

def read_labware_configuration(filepath='./preparation/labware_configuration.csv'):
    """
    This function reads the labware configuration from a csv file, and returns the data as a list of DataFrames.

    Parameters
    ----------
    filepath: str, path to the labware configuration file, default is './preparation/labware_configuration.csv'

    Returns
    -------
    dfs: list, list of DataFrames containing the labware configuration data
        ... Typically, there are 4 dataframes: [df_plate, df_pipetteLeft, df_pipetteRight, df_stocks]
    """

    # Function to parse a section of the data into a DataFrame
    def parse_section(data_section):
        """
        This function parses a section of the data into a DataFrame.

        Parameters
        ----------
        data_section: list, section of the data containing the labware configuration

        Returns
        -------
        df: DataFrame, labware configuration dataframe
        """
        # Convert the section of data into a StringIO object to use with read_csv
        data_string = StringIO("\n".join(data_section))
        return pd.read_csv(data_string, header=None)

    def expand_string_to_list(s):
        """
        This function expands a string to a list of elements
        if the string contains ',' or '-'.

        Parameters
        ----------
        s: str, string to be expanded
        ... Examples: 'A1', 'A1, D2', 'A1-A12', 'A1, B5-F9'

        Returns
        -------
        elements or s: list or str, list of elements or the string itself
        """
        if ',' in s or '-' in s:
            elements = []
            for part in s.split(','):
                part = part.strip()
                if '-' in part:
                    start, end = part.split('-')
                    start_prefix = re.search(r'[A-Za-z]+', start).group()
                    start_num = int(re.search(r'\d+', start).group())
                    end_prefix = re.search(r'[A-Za-z]+', end).group()
                    end_num = int(re.search(r'\d+', end).group())

                    if start_prefix == end_prefix:
                        for num in range(start_num, end_num + 1):
                            elements.append(f'{start_prefix}{num}')
                    else:
                        # Handle if the prefix changes
                        for prefix in range(ord(start_prefix), ord(end_prefix) + 1):
                            # Handle the start prefix
                            if prefix == ord(start_prefix):
                                for num in range(start_num, 13):
                                    elements.append(f'{chr(prefix)}{num}')
                            # Handle middle prefixes
                            elif prefix > ord(start_prefix) and prefix < ord(end_prefix):
                                for num in range(1, 13):
                                    elements.append(f'{chr(prefix)}{num}')
                            # Handle the end prefix
                            elif prefix == ord(end_prefix):
                                for num in range(1, end_num + 1):
                                    elements.append(f'{chr(prefix)}{num}')
                else:
                    elements.append(part)
            return elements
        else:
            return s

    def format_volume_in_df(df_stock):
        """
        This function formats the volume_mL column in the DataFrame
        to match the data type of the `rack_location` column.

        Parameters
        ----------
        df_stock: DataFrame, labware configuration dataframe (stock solutions)

        Returns
        -------
        df_stock: DataFrame, labware configuration dataframe (stock solutions)
                  with formatted volume_mL column

        """
        # If rack_location is a list, then convert the volume into a list of volumes
        # Check for each row if rack_location is a list
        for i, row in df_stock.iterrows():
            if isinstance(row['rack_location'], list):
                rack_location = row['rack_location']
                nStockSolutions = len(rack_location)
                vol = row['volume_mL']
                # Update the volume_mL column
                df_stock['volume_mL'] = df_stock['volume_mL'].astype(object)
                df_stock.at[i, 'volume_mL'] = [vol] * nStockSolutions
        return df_stock

    def organize_labware_dataframes(df):
        """
        This function organizes the labware dataframes by setting the correct headers and cleaning them.

        Parameters
        ----------
        df: DataFrame, labware configuration dataframe

        Returns
        -------
        df: DataFrame, cleaned and organized labware configuration dataframe
        """
        # Set the correct headers for each dataframe and clean them
        ncols = df.shape[1] - df.iloc[0].isnull().sum() # Number of valid columns (name, deck_location, starting_well_position)
        nrows = df.shape[0] - df.iloc[:, 0].isnull().sum() # Number of valid rows
        df.columns = df.iloc[0]
        df = df.iloc[1:nrows, 0:ncols]  # Select only valid rows and columns

        # if df contains chemical_name, then use that as the index
        if 'chemical_name' in df.columns:
            df = df.set_index('chemical_name')
        # Convert the volume_mL column to float
        for col in df.columns:
            if col in ['volume_mL', 'concentration',
                       'dilution_factor', 'volume_source_mL',
                       'volume_dest_mL', 'repeat'
                       ]:
                df[col] = df[col].astype(float)
        # Convert rack_location to a list if it contains ',' or '-'
        ## This function converts 'A1-B12' to ['A1', 'A2', 'A12', 'B1', ..., 'B12']
        if 'rack_location' in df.columns:
            df['rack_location'] = df['rack_location'].apply(expand_string_to_list)
            if 'volume_mL' in df.columns:
                df = format_volume_in_df(df)
        return df


    # Load the file
    data = pd.read_csv(filepath, delimiter='\t', header=None)

    # Splitting the data into sections based on the 'type' header and parsing them into DataFrames
    sections = []
    current_section = []

    try:
        for row in data[0]:
            if row.startswith('type'):
                if current_section:
                    sections.append(current_section)
                    current_section = []
            current_section.append(row)

        # Add the last section if it exists
        if current_section:
            sections.append(current_section)

        # Parse each section into a DataFrame
        dataframes = [parse_section(section) for section in sections]

        # Organize the dataframes
        dfs = [organize_labware_dataframes(df) for df in dataframes]
        return dfs
    except:
        raise ValueError('Error: The labware configuration file is not formatted correctly. '
                         'Please check the file and try again.')



def create_sample_info(samples, filepath='./preparation/sample_info.csv'):
    # Output plate_number, well, sample name, and volume, composition, recipe, order_mixing
    # Sort them by the plate_number and well
    # Output the data to a csv file
    """
    This function creates a csv file containing the sample information.
    Parameters
    ----------
    samples: list, list of Sample objects
    filepath: str, path to save the sample information, default is './preparation/sample_info.csv'

    Returns
    -------

    """

    def natural_sort(df):
        """
        This function sorts the DataFrame by column 'plate_deck_loc' numerically and 'well_sorted' naturally.

        Parameters
        ----------
        df: DataFrame

        Returns
        -------
        df: DataFrame, sorted DataFrame
        """

        def alphanum_key(s):
            """ Turn a string into a list of string and number chunks.
                "A10" -> ("A", 10) """
            return tuple(int(text) if text.isdigit() else text for text in re.split('([0-9]+)', s))

        # Sort the DataFrame by column 'A' numerically and 'B' naturally
        df['well_sorted'] = df['well'].apply(alphanum_key)
        df.sort_values(by=['plate_deck_loc', 'well_sorted'], inplace=True)
        df.drop('well_sorted', axis=1, inplace=True)  # Clean up the temporary column
        return df
    # if samples is a list of list, flatten
    if isinstance(samples[0], list):
        samples = [item for sublist in samples for item in sublist]

    data = []
    for sample in samples:
        if isinstance(sample.well, list):
            for i, well in enumerate(sample.well):
                if i == len(sample.well)-1: footer = ''
                else: footer = '_prep'
                data.append([sample.plate_deck_loc, well, sample.id, sample.name + footer, sample.volume, sample.composition, sample.recipe, sample.order_mixing])
        else:
            data.append([sample.plate_deck_loc, sample.well, sample.id, sample.name, sample.volume, sample.composition, sample.recipe, sample.order_mixing])
    df = pd.DataFrame(data, columns=['plate_deck_loc', 'well', 'sample_id','name', 'volume', 'composition', 'recipe', 'order_mixing'])

    # natural sort df by plate_deck_loc and well
    df = natural_sort(df)

    df.to_csv(filepath, index=False)
    print(f'create_sample_info:\n...Sample information saved as {filepath}')
def update_dilution_info(filepath='./preparation/dilution_info.csv',
                       info=None, verbose=False, overwrite=False):
    # if the file does not exist, create a new file
    if info is None:
        info = {'plate_deck_loc': 10,
                'well': 'A1',
                'sample_id': None,
                'name': 'Sample_Dilution_x10',
                'volume': 100,
                'composition': None,
                'recipe': None,
                'order_mixing': ['sample', 'diluent'],
                }

    if not os.path.exists(filepath):
        df = pd.DataFrame(columns=['plate_deck_loc',
                                   'well',
                                   'sample_id',
                                   'name',
                                   'volume',
                                   'composition',
                                   'recipe',
                                   'order_mixing'
                                   ])
    else:
        df = pd.read_csv(filepath)  # update

    # read the last value of df[sample_id]
    if info['sample_id'] is None:
        try:
            info['sample_id'] = int(df.loc[:, 'sample_id'].values[-1]) + 1
        except:
            info['sample_id'] = 1

    # add a new row to the df
    df = df.append(info, ignore_index=True)
    # save the df
    df.to_csv(filepath, index=False)
    if verbose:
        print(f'... Updated the sample_info.csv file with the dilution info.')
    return df


# CREATING RUN.PKL FILE WHICH CONTAINS THE LABWARE CONFIGURATION AND SAMPLE LISTS ### COME BACK
def create_run_file(volume, sample_list_nos, stock_dict, solvent='water'):
    '''
    This function creates a run file for the OT2 based on the sample list and stock concentrations.
    The run file includes the pipetteTypes, tip boxes, stock concentrations, and sample lists that are required for the run.
    '''
    nPlates, pipetteTypes, tip_boxes, stock_vols = estimate_requirements(volume, sample_list_nos, stock_dict,
                                                                      show=False, solvent=solvent)

    create_default_labware_configuration(nPlates, pipetteTypes, tip_boxes, stock_dict, stock_vols)

    # save the run file
    filename = 'run_file-' + datetime.datetime.now().strftime("%Y-%m-%d__%H-%M-%S") + '.pkl'
    sample_lists = [sorted([f for f in os.listdir(sample_dir) if f.startswith('sample_list')])[i] for i in
                    sample_list_nos]

    with open(sample_dir / sample_lists[0], 'rb') as f:
        metadata, chemical_viscosity, samples = pickle.load(f)
        if chemical_viscosity is None:
            chemical_viscosity = samples[0].relative_viscosity

    while True:
        print('List of chemicals: {}'.format(chemical_viscosity.keys()))
        unmixed_chem_list = input('Enter the chemicals that need not be mixed immediately?').split()
        if set(unmixed_chem_list).issubset(stock_dict.keys()):
            break
        else:
            print('invalid input')
            continue

    if not os.path.exists(master_dir / 'runfiles'):
        os.makedirs('runfiles')
    with open(master_dir / 'runfiles' / filename, 'wb') as f:
        pickle.dump([chemical_viscosity, unmixed_chem_list, sample_lists], f)
    print('run file saved as {} in folder \'runfiles\''.format(filename))


# ESTIMATING REQUIRED LABWARE AND CHEMICALS
## TM: WORK ON THIS LATER

def estimate_required_tips(pipette1_vol_range, pipette2_vol_range, vols):
    """
    This function estimates the number of tips required for the run based on the volumes of the samples.

    Parameters
    ----------
    pipette1_vol_range: list, range of volumes that can be handled by pipette1
        ... e.g. [2, 20] for p20
    pipette2_vol_range: list, range of volumes that can be handled by pipette2
        ... e.g. [20, 300] for p300
    vols: list, list of volumes of the samples
        ... e.g. [10, 30, 350] for samples with volumes 10, 30, and 50 ul
    Returns
    -------
    tips: list, number of tips required for each pipette
    """

    tips = np.zeros(2)
    for vol in vols:
        if vol < pipette2_vol_range[0]:
            tips[0] += np.ceil(vol / pipette1_vol_range[1])
        elif vol > pipette2_vol_range[1]:
            tips[1] += np.ceil(vol / pipette2_vol_range[1])
        else:
            tips[1] += 1

    return tips

def estimate_required_volumes(conc_dict, stock_dict, vol_sample):
    """
    Function to calculate the volumes of the chemicals and water that need to be added to the sample to make it.

    Parameters
    ----------
    conc_dict: dictionary, concentrations of the chemicals in the sample
    stock_dict: dictionary, stock concentrations of the chemicals
    vol_sample: float, volume of the sample in uL

    Returns
    -------
    vol_dict: dictionary, volumes of the chemicals and water that need to be added to the sample to make it
    """

    vol_dict = {}
    for key in conc_dict:
        vol_dict[key] = [float(conc_dict[key][0]) / float(stock_dict[key]) * vol_sample,
                         conc_dict[key][1]]  # Volumes in uL
    # make up remaining volume with water
    vol_dict['water'] = [vol_sample - sum(vol_dict.values()), 0]
    return vol_dict


# Managing Samples (A list of Sample objects)
def areSampleIDsValid(samples):
    """
    This function checks if the sample IDs are valid.
    ... This function checks if the sample IDs are unique and start from 1.
    Parameters
    ----------
    samples: list, list of Sample objects

    Returns
    -------
    bool: True if the sample IDs are valid, False otherwise
    """
    ans = True
    # Flatten the list if necessary
    if isinstance(samples[0], list):
        samples = [item for sublist in samples for item in sublist]
    # Grab sample IDs
    ids = [sample.id for sample in samples]
    if len(ids) != len(set(ids)):
        ans = False
        print('Error: Duplicate sample IDs found in the samples')
    if ids != list(range(1, len(ids) + 1)):
        ans = False
        print('Error: Invalid sample IDs found in the samples')
    return ans
def assignSampleIDs(samples, verbose=True):
    """
    This function assigns sample IDs to the samples.
    ... The sample IDs start from 1 and are assigned in the order of the samples.

    Parameters
    ----------
    samples: list, list of Sample objects

    Returns
    -------
    None
    """
    # Flatten the list if necessary
    if isinstance(samples[0], list):
        samples = [item for sublist in samples for item in sublist]
    for i in range(len(samples)):
        samples[i].id = 1 + i
    if verbose:
        print('assignSampleIDs:\n... Sample IDs were assigned')

def assignWellPositions(samples, method='normal', nWellsPerPlate=96,
                        initialize_well_volume=True,
                        filepath_lab_config='./preparation/labware_configuration.csv',
                        save=True):
    """
    This function assigns well positions to the samples.

    Parameters
    ----------
    samples: list, list of Sample objects
    method: str, method to assign well positions, default is 'normal'
        ... Choose from ['normal', 'alternate']
    nWellsPerPlate: int, number of wells per plate, default is 96
    starts_with: str, starting well position, default is 'A1'
    filepath: str, path to the labware configuration file,
                   default is './preparation/labware_configuration.csv'
        ... Used to get the deck locations of the well plates
    save: bool, whether to save the samples or not, default is True
        ... If True, this will UPDATE the well positions in the sample files.

    Returns
    -------
    samples: list, list of Sample objects with well positions assigned
    """
    def getAvailableWellPlateDeckLocations(df_plate):
        """
        This function returns the deck locations (int) of the available well plates.
        ... This is read from the labware configuration file.

        Parameters
        ----------
        df_plate: DataFrame, labware configuration dataframe for the plates
            ... Typically, df_plate = read_labware_configuration(filepath)[0]

        Returns
        -------
        plateDeckLocations: list, list of deck locations of the available well plates
            ... e.g. [10, 9] # deck locations of the first and second well plates
        """
        plateDeckLocations = df_plate['deck_location'].values.astype(int)
        return plateDeckLocations

    def getReqNumWells(sample):
        """
        This function calculates the number of wells required to prepare a sample.

        Parameters
        ----------
        sample: Sample object
        nWellsPerPlate: int, number of wells per plate

        Returns
        -------
        nWellsPerSample: int, number of wells required to prepare a sample
        """

        def count_lists(lst):
            """ Count the number of lists in a nested list. """
            count = 1  # Start with 1 to count the current list
            for item in lst:
                if isinstance(item, list):
                    count += count_lists(item)  # Recursively count sublists
            return count
        # Number of required wells to prepare a sample is equal to the number of lists in the order_mixing
        # e.g. ['fus', ['fus_buffer', 'proline', 'water']] requires 2 wells to prepare the sample
        # ... This reads 1. Transfer 'fus' to Well A
        #                2. Transfer 'fus_buffer', 'proline', 'water' to Well B.
        #                3. Transfer liquid in Well B to Well A.
        nReqWells = count_lists(sample.order_mixing)
        return nReqWells

    # Get information about the wellplates from the labware configuration file
    df_plate = read_labware_configuration(filepath=filepath_lab_config)[0]
    # Get the deck locations of the well plates
    plateDeckLocations = getAvailableWellPlateDeckLocations(df_plate)
    # Get the starting well position
    starts_with = df_plate['starting_well_position'].values[0]

    # Flatten the list if necessary
    if isinstance(samples[0], list):
        # Keep track the number of samples in each list for later
        nSamplesInEachList = [len(sublist) for sublist in samples]
        # Flatten the list
        samples = [item for sublist in samples for item in sublist]
    else:
        nSamplesInEachList = [len(samples)]
    nSamples = len(samples) # Total number of samples to be prepared

    # Check if the sample IDs are valid (no duplicates and must start from 1)
    if not areSampleIDsValid(samples):
        print('... Assigning new sample IDs')
        assignSampleIDs(samples, verbose=False)
        print('...... Done')

    # Assign the well plate and well positions to the samples
    for i, sample in enumerate(samples):
        # Deck location of the well plate for the sample
        sample.plate_deck_loc = plateDeckLocations[int(i // nWellsPerPlate)]
        # Number of wells required to prepare a sample
        nWellsPerSample = getReqNumWells(sample)
        # Assign the well position(s) to the sample
        # Some samples require more than 1 well to prepare
        sample.well = getWellPosition(1, nWellsPerSample,
                                      method=method,
                                      starts_with=starts_with)
        starts_with = get_next_well(sample.well)

        # For samples with complex order_mixing, assign sub-well positions.
        if isinstance(sample.well, list):
            assignSubWells(sample)

        if initialize_well_volume:
            if isinstance(sample.well, str):
                sample.well_volume = 0
            elif isinstance(sample.well, list):
                sample.well_volume = dict(zip(sample.well, [0]*len(sample.well)))

    if save:
        # Figure out the filepaths to save the samples with updated well positions
        all_samples, filepaths, m = [], [], 0
        for i in range(len(nSamplesInEachList)):
            n = sum(nSamplesInEachList[:i+1])
            all_samples.append(samples[m:n])
            filepaths.append(samples[m].filepath)
            m = n
        # Save the samples with updated well positions
        for i, samples in enumerate(all_samples):
            if filepaths[i] == '':
                raise ValueError('Error: Filepath is not provided in the samples. samples[i].filepath is empty. '
                                 'Update it or set save=False')
            else:
                save_samples(samples, verbose=False, filepath=filepaths[i])
    print('assignWellPositions:\n... Well positions were assigned to the samples successfully')


def assignSubWells(sample):
    """
    This function assigns well positions when the sample requires more than 1 well to prepare.
    ... This occurs when Sample.order_mixing contains nested lists.
    ... For example, ['fus', ['fus_buffer', 'proline', 'water']] requires 2 wells to prepare the sample
    ...... This reads 1. Transfer 'fus' to Well A
    ......            2. Transfer 'fus_buffer', 'proline', 'water' to Well B
    ...... In such cases, assignWellPositions assigns multiple well positions (list) to the sample.
    ...... This function assigns the sub-well positions to the sample.

    Parameters
    ----------
    sample: Sample object

    Returns
    -------
    chem2well: dict, dictionary containing the sub-well positions assigned to the sample
        ... subwell = {chemical1: well1, chemical2: well2, ...}
        ... Example: {'fus': 'A1', 'fus_buffer': 'A2', 'proline': 'A2', 'water': 'A2'}
    """
    wells = sample.well
    chem2well = {}
    order = sample.order_mixing

    count = 0
    while count < len(wells):
        deepest_lists = find_deepest_lists(order)
        for i, lst in enumerate(deepest_lists):
            for chem in lst:
                if chem not in chem2well.keys():
                    chem2well[chem] = wells[count]
            count += 1

        order_flt_in = flatten_innermost(order)
        order = order_flt_in

    sample.subwell = chem2well # This dictionary assigns the sub-well positions to the sample

    return chem2well

def get_available_well_positions(samples, nWellsPerPlate=96, ignorePriorWells=True):
    """
    Function to generate a list of available well positions for 96-well plates based on the samples.

    Parameters
    ----------
    nWellsPerPlate: int, number of wells per plate, default is 96

    Returns
    -------
    available_wells: list, list of available well positions
    """

    # Flatten the list if necessary
    if isinstance(samples[0], list):
        samples = [item for sublist in samples for item in sublist]

    well_positions = natural_sort([sample.well for sample in samples])
    available_wells = [f'{chr(65 + i // 12)}{1 + i % 12}' for i in range(nWellsPerPlate)]

    if ignorePriorWells:
        start = available_wells.index(well_positions[0])
        available_wells = available_wells[start:]

    for well in well_positions:
        available_wells.remove(well)

    return available_wells


## HELPERS
def isBeforeOrEqualTo(well1, well2):
    """
    This function checks if well1 is before well2.

    Parameters
    ----------
    well1: str, well position e.g. 'A1'
    well2: str, well position e.g. 'A2'

    Returns
    -------
    bool: True if well1 is before well2, False otherwise
    """
    return ord(well1[0])+int(well1[1:]) <= ord(well2[0])+int(well2[1:])

def get_recipe(volume, sample, stock_dict, solvent='water'):
    """
    This function calculates the recipe for the sample based on the stock concentrations and the composition of the sample.

    Parameters
    ----------
    volume: float, volume of the sample to be prepared
    sample: Sample object
    stock_dict: dictionary containing the stock concentrations of the chemicals
    solvent: string, solvent used in the sample, default is water

    Returns
    -------
    recipe: dictionary containing the recipe for the sample
        ... recipe = {chemical1: volume in uL, chemical2:, volume in uL, ... , solvent: volume in uL}
    """
    recipe = {}
    for key in sample.composition:
        if key == solvent:
            continue
        else:
            recipe[key] = volume * sample.composition[key] / stock_dict[key]
    recipe[solvent] = volume - sum(recipe.values())
    return recipe


def calculate_max_concentration(sample, varying_chemical_species):
    """
    This function calculates the maximum concentration of the varying chemical species that can be added to the sample.
    ... The maximum concentration is calculated based on the stock concentrations of the chemicals in the sample.

    Parameters
    ----------
    sample: Sample, a sample object
    varying_chemical_species: str or list, chemical species that are varying in the sample
    ... e.g. 'proline' or ['proline', 'glycine']

    Returns
    -------
    max_concs: dict, dictionary containing the maximum concentrations of the varying chemical species
    """
    if isinstance(varying_chemical_species, str):
        varying_chemical_species = [varying_chemical_species]
    stock_dict, composition = sample.stock_dict, sample.composition
    chemicals, solvent = sample.chemicals, sample.solvent
    fixed_chems = [species for species in chemicals if species not in varying_chemical_species + [solvent]]
    max_concs = {}
    V = 1000
    vol = 0
    for species in fixed_chems:
        if species in stock_dict.keys() and species in composition.keys():
            vol += composition[species] * V / stock_dict[species]
        else:
            raise KeyError(f'{species} is not in the sample.stock_dict, {list(stock_dict.keys())} '
                           f'or sample.composition, {list(composition.keys())}')
    remainder_chems = [chem for chem in chemicals if chem not in fixed_chems + [solvent]]
    for chem in remainder_chems:
        max_concs[chem] = (V - vol) * stock_dict[chem] / V
    return max_concs

def get_relative_viscosity(chemical_list, viscous=['peg'], special=['bsa']):
    """
    This function determines the relative viscosity of the chemicals in the stock_dict.
    ... The relative viscosity is classified as 'thin' or 'thick' except special species.
    ... Special species are classified as 'bsa' for example.
    ...... This information is useful for the OT2 to determine the pipetting strategy.

    Parameters
    ----------
    chemical_dict: dict, dictionary containing the stock concentrations of the chemicals
    viscous: list, list of chemicals that are viscous. MUST BE LOWER CASE e.g. 'peg'
    special: list, list of special chemicals that need to be treated differently. MUST BE LOWER CASE e.g. 'bsa'
        ... Special species are classified by their chemical name. e.g. 'bsa'

    Returns
    -------
    relative_viscosity: dict, dictionary containing the relative viscosity of the chemicals
        ... e.g. {'bsa': 'bsa', 'peg': 'thick', 'water': 'thin'}
    """
    relative_viscosity = {}
    for species in chemical_list:
        # If it is viscous as PEG4k stock (60%w/v), it is considered 'viscous'.
        if species.lower() in viscous:
            relative_viscosity[species] = 'thick'  # e.g. 'peg'
        elif species.lower() in special:
            relative_viscosity[species] = species  # e.g. 'bsa'
        else:
            relative_viscosity[species] = 'thin'  # e.g. 'water'
    return relative_viscosity

def getWellPosition(sample_id, nWellsPerID=1, method='normal', starts_with='A1', nWellsPerPlate=96):
    """
    Function to generate a unique well position for a given sample number
    ... this is useful when you want to generate a unique well position for a sample in a 96-well plate
    ... Sample number 1 -> A1, Sample number 2 -> A2, ... Sample number 13 -> B1, Sample number 14 -> B2, ...

    Parameters
    ----------
    sample_number: int, sample number
    method: str, method to generate the well position, 'all' or 'diagonal'
        ... 'normal':    A1, A2, ..., A12,  B1,  B2, ..., B12, C1, C2, ..., C12, D1,   D2, ..., D12, ...
        ... 'alternate': A1, A2, ..., A12, B12, B11, ...,  B1, C1, C2, ..., C12, D12, D11, ...,  D1, ...
    starts_with: str, starting well position

    Returns
    -------
    well_position: str or list, well position for the given sample number
        ... e.g. 'A1', 'B2', or ['A1', 'A2', 'A3']
    """
    positions = []
    # Convert starts_with to a numerical index
    start_row = ord(starts_with[0]) - ord('A')
    start_col = int(starts_with[1:]) - 1
    start_offset = start_row * 12 + start_col

    # Calculate the start position taking into account the number of wells per ID and the starting offset
    start_index = start_offset + (sample_id - 1) * nWellsPerID

    for i in range(nWellsPerID):
        current_index = (start_index + i) % nWellsPerPlate
        row = current_index // 12
        col = current_index % 12 + 1

        if method == 'normal':
            positions.append(chr(65 + row) + str(col))
        elif method == 'alternate':
            if row % 2 == 0:
                positions.append(chr(65 + row) + str(col))
            else:
                positions.append(chr(65 + row) + str(13 - col))
    if len(positions) == 1:
        return positions[0]
    else:
        return positions


def get_next_well(current_well, wells=_wells):
    """
    Function to get the next well position in a 96-well plate.

    Parameters
    ----------
    current_well: str or list, current well position
        ... If it is a list, it will return the next well position for the last element in the list
    wells: itertools.Cycle, cycle object containing the well positions in a 96-well plate
        ... [A1, A2, ..., A12, B1, B2, ..., B12, ..., H12]

    Returns
    -------
    next_well: str, next well position
    """
    if isinstance(current_well, list):
        current_well = current_well[-1]
    for well in wells:
        if well == current_well:
            return next(wells)

def get_well_cycler(starts_with='A1'):
    """
    Function to get a cycle object for the well positions in a 96-well plate.

    Parameters
    ----------
    starts_with: str, starting well position, default is 'A1'

    Returns
    -------
    well_cycler: itertools.Cycle, cycle object containing the well positions in a 96-well plate
        ... [A1, A2, ..., A12, B1, B2, ..., B12, ..., H12]
    """
    well_cycler = itertools.cycle([f"{letter}{number}" for letter in 'ABCDEFGH' for number in range(1, 13)])  # 96-well plate

    # Convert starts_with to a numerical index
    start_row = ord(starts_with[0]) - ord('A')
    start_col = int(starts_with[1:]) - 1
    start_offset = start_row * 12 + start_col

    # Calculate the start position taking into account the number of wells per ID and the starting offset
    start_index = start_offset
    for i in range(start_index):
        next(well_cycler)
    return well_cycler

def get_wells(n=10, start='A1', end=None, inc=1, nWellsPerPlate=96):
    """
    Function to get a list of well positions in a 96-well plate.
    Parameters
    ----------
    n: int, number of well positions to get
    start: str, starting well position, default is 'A1'
    end: str, ending well position, default is None.
         If not None, it will return well positions from start to end
    inc: int, increment, default is 1
    nWellsPerPlate: int, number of wells per plate, default is 96

    Returns
    -------

    """
    if nWellsPerPlate != 96:
        raise ValueError('This function is only implemented for 96-well plates')

    well_cycler = get_well_cycler(starts_with=start)
    wells = []
    if end is None:
        # get n wells from a cycler starting from start
        for i in range(0, n*inc):
            well = next(well_cycler)
            if i % inc == 0:
                wells.append(well)
    else:
        # get wells from a cycler starting from start to end
        count = 0
        while True and count < nWellsPerPlate:
            well = next(well_cycler)
            if count % inc == 0:
                wells.append(well)
            if well == end:
                break
            count += 1
    return wells

# Height-to-volume conversion functions
# Volume-to-Height conversion functions
## Height from the BOTTOM of the tube
def v2z_deepwellplate(vol):
    '''
    Parameters: vol in uL
    Returns: height
    Description: Converts volume to height in deep well plates
    '''
    area = np.pi * (7.85 / 2) ** 2
    height = vol / area
    return max(height - 1, 1)  # height from the bottom of the well

def v2z_200uLwellplate(vol):
    '''
    Parameters: vol in uL
    Returns: height
    Description: Converts volume to height in deep well plates
    '''
    area = np.pi * 2.5 ** 2
    height = vol / area * (12 / 10.2) * 1.25
    return max(height, 1)  # height from the bottom of the well

def v2z_200uLeppendorf(vol):
    """
    Convert volume to height in 200uL eppendorf tubes
    ... Used for stock solutions of purified proteins stored in 200uL eppendorf tubes (e.g. FUS)

    Parameters
    ----------
    vol: float, volume in uL

    Returns
    -------
    float, height in mm
    """
    area = np.pi * 2.5 ** 2
    height = vol / area
    return max(height, 1)  # height from the bottom of the well

def v2z_falcon15ml(vol):
    '''
    Parameters: vol in ml for now should be more than 2ml
    Returns: height
    Description: Converts volume to height in falcon 15ml tubes
    '''
    vol = max(vol, 2)

    angle = 0.315  # radians
    vol = vol - 1.8
    height = vol * 1000 / (np.pi * 7.33 ** 2)
    height += 26.48

    h_from_bottom = max(height, 1)
    return h_from_bottom  # height from the bottom of the tube

def v2z_falcon50ml(vol):
    '''
    Parameters: vol in ml does not work for less than 5ml
    Returns: height
    Description: Converts volume to height in falcon 50ml tubes
    '''
    tanTheta = 1
    if vol < 5:
        height = (3 * vol * 1000 / tanTheta ** 2 / np.pi) ** (1 / 3) - 2.4
    else:
        vol = vol - 5
        height = vol * 1000 / (np.pi * (13.5) ** 2)
        height += 16
    h_from_bottom = max(height, 1)
    return h_from_bottom  # height from the bottom of the tube

## DEPRICATED due to the naming convention. VtoH is replaced with v2z.
def falcon15ml_VtoH(vol):
    '''
    Parameters: vol in ml for now should be more than 2ml
    Returns: height
    Description: Converts volume to height in falcon 15ml tubes
    '''
    print('DEPRICATED: Use v2z_falcon15ml')
    angle = 0.315  # radians
    if vol < 1:
        height = np.power(3 * vol / (np.pi * np.tan(angle) ** 2), 1 / 3) * 10  # height in mm
        raise ValueError('Volume too low V = {} should be more than 2ml'.format(vol))
    else:
        vol = vol - 1.8
        height = vol * 1000 / (np.pi * 7.33 ** 2)
        height += 26.48
    return max(height, 1)  # height from the bottom of the tube

def falcon50ml_VtoH(vol):
    '''
    Parameters: vol in ml does not work for less than 5ml
    Returns: height
    Description: Converts volume to height in falcon 50ml tubes
    '''
    print('DEPRICATED: Use v2z_falcon50ml')
    tanTheta = 1
    if vol < 5:
        height = (3 * vol * 1000 / tanTheta ** 2 / np.pi) ** (1 / 3) - 2.4

    else:
        vol = vol - 5
        height = vol * 1000 / (np.pi * (13.5) ** 2)
        height += 16
    return max(height, 1)  # height from the bottom of the tube

# GENERAL HELPERS
def extract_str(s, start, end):
    """
    Extract a string from a string like '_lalala_start5p3end_lalala' -> 5.3 (float)
    Parameters
    ----------
    s: string
    start: string, start of the string to extract
    end: string, end of the string to extract

    Returns
    -------
    str_extracted: string extracted from the input string

    """
    if s.find(start) < 0:
        print('ERROR: ' + start + ' is not in ' + s)
        raise NameError('issues with ' + start)
    if s.find(end) < 0:
        print('ERROR: ' + end + ' is not in ' + s)
        raise NameError('issues with ' + end)
    start_ind, end_ind = s.find(start), s.find(end)
    str_extracted = s[start_ind + len(start):end_ind]
    return str_extracted

def natural_sort(arr):
    """
    Natural-sorts elements in a given array
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)

    e.g.-  arr = ['a28', 'a01', 'a100', 'a5']
    ... WITHOUT natural sorting,
     -> ['a01', 'a100', 'a28', 'a5']
    ... WITH natural sorting,
     -> ['a01', 'a5', 'a28', 'a100']


    Parameters
    ----------
    arr: list or numpy array of strings

    Returns
    -------
    sorted_array: natural-sorted
    """

    def atoi(text):
        'natural sorting'
        return int(text) if text.isdigit() else text

    def natural_keys(text):
        return [atoi(c) for c in re.split('(\d+)', text)]

    return sorted(arr, key=natural_keys)

## NESTED LIST HELPERS
def flatten_list(nested_list):
    """
    This function flattens a nested list.

    Parameters
    ----------
    nested_list: list, nested list, e.g. [[1, 2], [[3,4,5], 6], 7]

    Returns
    -------
    flat_list: list, flattened list, e.g. [1, 2, 3, 4, 5, 6, 7]

    """
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten_list(item))
        else:
            flat_list.append(item)
    return flat_list

def flatten_innermost(nested_list):
    """
    This function flattens the innermost lists in a nested list.

    Parameters
    ----------
    nested_list: list, nested list

    Returns
    -------
    flattened_list: list, list with innermost lists flattened

    Examples
    --------
    order = ['1', ['2', ['3', '4'], ['5', '6'], '7'], '8', ['9', '10'], '11']
    -> flatten_innermost(order) returns ['1', ['2', '3', '4', '5', '6', '7'], '8', ['9', '10'], '11']
    """
    def max_depth(lst, depth=1):
        # This function determines the maximum depth of the nested list.
        if isinstance(lst, list):
            return max([depth] + [max_depth(item, depth + 1) for item in lst])
        else:
            return depth - 1

    def flatten_at_max_depth(lst, depth=1, target_depth=None):
        if target_depth is None:
            target_depth = max_depth(lst)  # Calculate maximum depth
        new_list = []
        for item in lst:
            if isinstance(item, list):
                if depth == target_depth - 1:
                    # Flatten this level
                    new_item = '_'.join(item) # Rename the innermost list as a proof of flattening
                    new_list.append(new_item) # [1, [2, 3, 4], 5] -> [1, '2_3_4', 5]
                    # Alternatively, you can extend the list
#                     new_list.extend(new_item) # [1, [2, 3, 4], 5] -> [1, 2, 3, 4, 5]
                else:
                    # Continue processing deeper lists
                    new_list.append(flatten_at_max_depth(item, depth + 1, target_depth))
            else:
                new_list.append(item)
        return new_list

    return flatten_at_max_depth(nested_list)

def find_deepest_lists(nested_list):
    """
    This function returns the deepest (innermost) lists in the nested list.
    ... If the given list is not nested, it returns the [original_list] (a list of a list)

    Parameters
    ----------
    nested_list: list, nested list

    Returns
    -------
    deepest_lists: list, list of deepest lists

    Examples
    --------
             order = ['1', ['2', ['3', '4'], ['5', '6'], '7'], '8', ['9', '10'], '11']
             -> find_deepest_lists(order) returns [['3', '4'], ['5', '6']]
    """
    def max_depth(lst, current_depth=0):
        """Finds the maximum depth of the nested list."""
        if isinstance(lst, list):
            return max((max_depth(item, current_depth + 1) for item in lst), default=current_depth)
        return current_depth - 1

    def collect_deepest_lists(lst, current_depth=0, target_depth=None):
        """Collects lists that exist at the maximum depth."""
        if target_depth is None:
            target_depth = max_depth(nested_list)
        if target_depth == 0:
            return [lst]
        else:
            deepest_lists = []
            for item in lst:
                if isinstance(item, list):
                    if current_depth == target_depth - 1:
                        deepest_lists.append(item)  # This list is at the deepest level
                    else:
                        deepest_lists.extend(collect_deepest_lists(item, current_depth + 1, target_depth))
            return deepest_lists

    return collect_deepest_lists(nested_list)


# DEPRICATED FUNCTIONS
# FUNCTIONS THAT ARE NEEDS TO BE CLEANED UP
## Currently, these functions are too long and need to be broken down into smaller functions
## These need to go. They are too long and hard to read. - Takumi


## Experient helpers
## These are too specific, and should be improved to keep them in this module.

def make_sample_from_csv(file, chemical_viscosity):
    '''
    This function reads the csv file and creates a sample object
    The first row of the csv file should be the header which contains the name of the chemical used.
    The first column should contain the sample ID
    The columns should contain the different concentrations of the chemicals in the sample
    The order from left to right should be the order in which the chemicals will be mixed

    Input:
    file: The path to the csv file

    Output:
    sample_list: A list of sample objects to be prepared by OT2
    '''
    sample_list = []
    with open(file, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        chemicals = header[1:]
        print(chemicals)
        for row in reader:
            composition = {chemicals[i]: float(row[i + 1]) for i in range(len(chemicals))}
            order_mixing = chemicals
            sample = Sample(composition, order_mixing)
            sample.id = int(row[0])
            sample_list.append(sample)
    name_list = input('what do you want to call the sample list?')
    comments = input('any comments on the samples?')
    time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    metadata = [name_list, time, comments]
    filename = 'sample_list-' + datetime.datetime.now().strftime("%Y-%m-%d__%H-%M-%S") + '.pkl'
    if not os.path.exists(sample_dir):
        os.makedirs('samples')
    with open(sample_dir / filename, 'wb') as f:
        pickle.dump([metadata, chemical_viscosity, sample_list], f)

    print('sample list saved as {}'.format(filename))

def make_sample_list_manually(chemical_viscosity, stock_dict, n=6, dst_dir='./samples'):
    '''
    This function makes a list of samples with variation of a given chemical.
    This is mainly if the sample list is pretty simple to make and involves only variation of one chemical while keeping the rest constant.
    The user has to input various details of the samples such as the chemicals, their composition, etc and the sample list is saved as a pickle file.

    Inputs:
    User inputs the details of the sample list such as the chemicals, their composition, etc.

    '''


    sample_list = []
    print('List of chemicals: {}'.format(chemical_viscosity.keys()))
    while True:
        str_list = input(
            'Enter list of set chemicals from chemical list {}'.format(chemical_viscosity.keys())).split()
        sample_chem_list = [i.lower().replace(' ', '') for i in str_list]
        sample_conc_list = []

        for chem in sample_chem_list:
            sample_conc_list.append(float(input('Set conc of {} in mM or g/L (for bsa and peg)'.format(chem))))

        var = input('Enter the chemical for variation from {}').format(chemical_viscosity.keys()).lower().replace(' ',
                                                                                                             '')
        var_max_conc = list(estimate_max_concentration(stock_dict, chemical_list=sample_chem_list + [var],
                                      composition=dict(zip(sample_chem_list, sample_conc_list)),
                                      fixed_chems=sample_chem_list).values())
        print(
            f'Suggested set of concentrations for {var} is {[math.floor(num) for num in np.linspace(0, var_max_conc[0], n)]}'.replace(
                ',', ''))
        while True:
            var_max_conc = list(estimate_max_concentration(stock_dict, chemical_list=sample_chem_list + [var],
                                          composition=dict(zip(sample_chem_list, sample_conc_list)),
                                          fixed_chems=sample_chem_list).values())
            str_list = input(
                'List of concs for variation of {} in mM or g/L (for bsa and peg)'.format(var)).split()
            conc_list = [float(i) for i in str_list]
            max_conc = max(conc_list)
            if all(max_conc <= x for x in var_max_conc):
                break
            else:
                print('Max conc allowed is {}'.format(var_max_conc))
        sample_chem_list.append(var)

        while True:
            str_list = input('list the order of mixing the chemicals in your sample - {}'.format(
                sample_chem_list + ['water'])).split()
            order_mixing = [i.lower().replace(' ', '') for i in str_list]
            if set(order_mixing).issubset(sample_chem_list + ['water']) and len(order_mixing) == len(
                    sample_chem_list) + 1:
                break
            else:
                print('invalid order of mixing')
                continue

        for conc in conc_list:
            composition = dict(zip(sample_chem_list, sample_conc_list))
            sample = Sample(composition, order_mixing)
            # sample.composition = dict(zip(sample_chem_list, sample_conc_list))
            # sample.composition[var] = conc
            sample.composition['water'] = 0
            sample.order_mixing = order_mixing
            sample_list.append(sample)

        if_bckgr = input('do you want to make background samples? (y/n)')
        if if_bckgr == 'y':
            chem_bckgr = input('which chemical do you want to make background samples for?').lower().replace(' ',
                                                                                                                  '')
            for conc in conc_list:
                sample = Sample(conc, order_mixing)
                sample.composition = dict(zip(sample_chem_list, sample_conc_list))
                sample.composition[chem_bckgr] = 0
                sample.composition['water'] = 0
                sample.order_mixing = order_mixing
                sample_list.append(sample)

        for sample in sample_list:
            print(sample.composition)
        print('order of mixing: {}'.format(sample_list[0].order_mixing))

        if input('confirm sample list? (y/n)') == 'y':
            break
        else:
            if input('restart sample list? (y/n)') == 'y':
                print('\n\n******* restarting sample list *******\n\n')
                sample_list = []
            else:
                return

    name_list = input('what do you want to call the sample list?')
    comments = input('any comments on the samples?')
    time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    metadata = [name_list, time, comments]
    filename = 'sample_list-' + name_list + '-' + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '.pkl'
    if not os.path.exists(dst_dir):
        os.makedirs('samples')
    with open(dst_dir / filename, 'wb') as f:
        pickle.dump([metadata, chemical_viscosity, sample_list], f)

    print('sample list saved as {}'.format(filename))

def make_sample_list_manually_bsa_peg(chemical_viscosity, stock_dict, n=6):
    '''
    This function makes a list of samples with variation of a given chemical.
    This is mainly if the sample list is pretty simple to make and involves only variation of one chemical while keeping the rest constant.
    The user has to input various details of the samples such as the chemicals, their composition, etc and the sample list is saved as a pickle file.

    Inputs:
    User inputs the details of the sample list such as the chemicals, their composition, etc.

    '''
    sample_list = []
    print('list of chemicals: {}'.format(chemical_viscosity.keys()))
    while True:
        ans = input('make default sample containing bsa,peg,kp,kcl? (y/n)')
        if ans == 'y':
            sample_chem_list = ['bsa', 'peg', 'kp', 'kcl']
            sample_conc_list = []
            var = input('chemical for variation {}').format(chemical_viscosity.keys()).lower().replace(' ',
                                                                                                            '')  # all chemicals are in lowercases and without spaces

            if var in sample_chem_list:
                sample_chem_list.remove(var)

            for chem in sample_chem_list:
                sample_conc_list.append(float(input('set conc of {} in mM or g/L (for bsa and peg)'.format(chem))))

            while True:
                var_max_conc = list(estimate_max_concentration(stock_dict, chemical_list=sample_chem_list + [var],
                                              composition=dict(zip(sample_chem_list, sample_conc_list)),
                                              fixed_chems=sample_chem_list).values())
                str_list = input(
                    'list of concs for variation of {} in mM or g/L (for bsa and peg)'.format(var)).split()
                conc_list = [float(i) for i in str_list]
                max_conc = max(conc_list)
                if all(max_conc <= x for x in var_max_conc):
                    break
                else:
                    print('max conc allowed is {}'.format(var_max_conc))
                    continue
            sample_chem_list.append(var)

            while True:
                str_list = input('list the order of mixing the chemicals in your sample - {}'.format(
                    sample_chem_list + ['water'])).split()
                order_mixing = [i.lower().replace(' ', '') for i in str_list]
                if set(order_mixing).issubset(sample_chem_list + ['water']) and len(order_mixing) == len(
                        sample_chem_list) + 1:
                    break
                else:
                    print('invalid order of mixing')
                    continue

        elif ans == 'n':
            str_list = input(
                'enter list of set chemicals from chemical list {}'.format(chemical_viscosity.keys())).split()
            sample_chem_list = [i.lower().replace(' ', '') for i in str_list]
            sample_conc_list = []
            for chem in sample_chem_list:
                sample_conc_list.append(float(input('set conc of {} in mM or g/L (for bsa and peg)'.format(chem))))

            var = input('chemical for variation from {}').format(chemical_viscosity.keys()).lower().replace(' ',
                                                                                                                 '')
            var_max_conc = list(estimate_max_concentration(stock_dict, chemical_list=sample_chem_list + [var],
                                          composition=dict(zip(sample_chem_list, sample_conc_list)),
                                          fixed_chems=sample_chem_list).values())
            print(
                f'Suggested set of concentrations for {var} is {[math.floor(num) for num in np.linspace(0, var_max_conc[0], n)]}'.replace(
                    ',', ''))
            str_list = input(
                'list of concs for variation of {} in mM or g/L (for bsa and peg)'.format(var)).split()
            conc_list = [float(i) for i in str_list]
            max_conc = max(conc_list)
            if all(max_conc <= x for x in var_max_conc):
                break
            else:
                print('max conc allowed is {}'.format(var_max_conc))
                continue
        else:
            print('invalid input')
            continue

        for conc in conc_list:
            sample = Sample()
            sample.composition = dict(zip(sample_chem_list, sample_conc_list))
            sample.composition[var] = conc
            sample.composition['water'] = 0
            sample.order_mixing = order_mixing
            sample_list.append(sample)

        if_bckgr = input('do you want to make background samples? (y/n)')
        if if_bckgr == 'y':
            chem_bckgr = input('which chemical do you want to make background samples for?').lower().replace(' ',
                                                                                                                  '')
            for conc in conc_list:
                sample = Sample()
                sample.composition = dict(zip(sample_chem_list, sample_conc_list))
                sample.composition[chem_bckgr] = 0
                sample.composition['water'] = 0
                sample.order_mixing = order_mixing
                sample_list.append(sample)

        for sample in sample_list:
            print(sample.composition)
        print('order of mixing: {}'.format(sample_list[0].order_mixing))

        if input('confirm sample list? (y/n)') == 'y':
            break
        else:
            if input('restart sample list? (y/n)') == 'y':
                print('\n\n******* restarting sample list *******\n\n')
                sample_list = []
            else:
                return

    name_list = input('what do you want to call the sample list?')
    comments = input('any comments on the samples?')
    time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    metadata = [name_list, time, comments]
    filename = 'sample_list-' + name_list + '-' + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '.pkl'
    if not os.path.exists(sample_dir):
        os.makedirs('samples')
    with open(sample_dir / filename, 'wb') as f:
        pickle.dump([metadata, chemical_viscosity, sample_list], f)

    print('sample list saved as {}'.format(filename))


## CUSTOM
### PURPOSE: TO MAKE A SAMPLE LIST WITH VARIATION OF A GIVEN CHEMICAL (BSA/PEG SUSCEPTIBILITY EXP)
def make_sample_list_bsa_peg_susceptibility(chemical_viscosity, stock_dict, solute, n=6,
                                        composition = {'bsa': 46.6, 'peg': 226., 'kp': 100., 'kcl': 200., 'water': 0.}):
    '''
    This function makes a list of samples with variation of a given chemical.
    This is mainly if the sample list is pretty simple to make and involves only variation of one chemical while keeping the rest constant.
    The user has to input various details of the samples such as the chemicals, their composition, etc and the sample list is saved as a pickle file.

    Inputs:
    User inputs the details of the sample list such as the chemicals, their composition, etc.

    '''


    sample_list = []
    print('list of chemicals: {}'.format(chemical_viscosity.keys()))
    while True:
        sample_conc_list = list(composition.items()) # list of concentrations of the chemicals in the sample
        sample_chem_list = list(composition.keys()) # list of chemicals in the sample
        var = input('chemical for variation {}').format(chemical_viscosity.keys()).lower().replace(' ', '')

        if var in sample_chem_list:
            sample_chem_list.remove(var)

        for chem in sample_chem_list:
            if chem not in ['bsa', 'peg']:
                sample_conc_list.append(float(input('set conc of {} in mM'.format(chem))))

        while True:
            sample = Sample(composition = {'fus_buffer': 50., solute: 0},
                         vary=solute, concentrations=[0, 25, 50, 100, 200, 300, 400, 500, 600, 700], solvent='water',
                         order_mixing=[solute, 'fus_buffer'],
                         stock_dict=stock_dict)
            var_max_conc = list(estimate_max_concentration(stock_dict, chemical_list=sample_chem_list + [var],
                                          composition=dict(zip(sample_chem_list, sample_conc_list)),
                                          fixed_chems=sample_chem_list).values())
            print(
                f'Suggested set of concentrations for {var} is {[math.floor(num) for num in np.linspace(0, var_max_conc[0], n)]}'.replace(
                    ',', ''))
            str_list = input(
                'list of concs for variation of {} in mM or g/L (for bsa and peg)'.format(var)).split()
            conc_list = [float(i) for i in str_list]
            max_conc = max(conc_list)
            if all(max_conc <= x for x in var_max_conc):
                break
            else:
                print('max conc allowed is {}'.format(var_max_conc))
                continue
        sample_chem_list.append(var)

        while True:
            # print(f'Suggested order of mixing is kp kcl water {var} bsa peg')
            # str_list = input('list the order of mixing the chemicals in your sample - {}'.format(sample_chem_list+['water'])).split()
            str_list = f'kp kcl water {var} bsa peg'.split()
            order_mixing = [i.lower().replace(' ', '') for i in str_list]
            if set(order_mixing).issubset(sample_chem_list + ['water']) and len(order_mixing) == len(
                    sample_chem_list) + 1:
                break
            else:
                print('invalid order of mixing')
                continue

        for conc in conc_list:
            sample = Sample()
            sample.composition = dict(zip(sample_chem_list, sample_conc_list))
            sample.composition[var] = conc
            sample.composition['water'] = 0
            sample.order_mixing = order_mixing
            sample_list.append(sample)

        for sample in sample_list:
            print(sample.composition)
        print('order of mixing: {}'.format(sample_list[0].order_mixing))

        if input('confirm sample list? (y/n)') == 'y':
            break
        else:
            if input('restart sample list? (y/n)') == 'y':
                print('\n\n******* restarting sample list *******\n\n')
                sample_list = []
            else:
                return

    # name_list = input('what do you want to call the sample list?')
    name_list = f'{var}_variation'
    # comments = input('any comments on the samples?')
    comments = ''
    time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    metadata = [name_list, time, comments]
    filename = 'sample_list-' + name_list + '-' + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S") + '.pkl'
    if not os.path.exists(sample_dir):
        os.makedirs('samples')
    with open(sample_dir / filename, 'wb') as f:
        pickle.dump([metadata, chemical_viscosity, sample_list], f)

    print('sample list saved as {}'.format(filename))




## add_chemical: Bad. water_transfer_and_mix, bsa_transfer_and_mix, PEG_transfer_and_mix ?
### Aren't they all the same except the mixing procedures? Solution- Call one function and pass different parameters.
def add_chemical(protocol, pipettes, stock_dict, out_well, chem_name, chem_type, vol_added, vol_well, basic_mix,
                 drop_tip=True):
    '''
    Function to add a particular chemical in the recipe based on its viscosity (chem_type)

    protocol: protocol_api.ProtocolContext object

    '''
    if vol_added > 0:
        stock_tube_rack = stock_dict[chem_name]['labware']
        if chem_type in ['water', 'thin']:
            water_transfer_and_mix(protocol, pipettes, vol_added, vol_well, stock_tube_rack, stock_dict[chem_name],
                                   out_well, basic_mix=basic_mix, drop_tip=drop_tip)
        elif chem_type in ['bsa']:
            bsa_transfer_and_mix(protocol, pipettes, vol_added, stock_tube_rack, stock_dict[chem_name],
                                 out_well, vol_well, drop_tip=drop_tip)
        elif chem_type == ['peg', 'thick', 'viscous']:
            PEG_transfer_and_mix(protocol, pipettes, vol_added, stock_tube_rack, stock_dict[chem_name],
                                 out_well, vol_well, drop_tip=drop_tip)
        else:
            raise ValueError('Chemical type not recognized')



def PEG_transfer_and_mix(protocol, pipettes, vol, stocks_tuberack, stock_dict, out_well, vol_well, mix=True,
                         drop_tip=True):
    '''
    Description: Transfers a volume of liquid with viscosity similar to PEG from the stock tube to the well and mixes the solution in the well.
    Input:
    pipettes = list of pipettes in the protocol
    vol = volume of the liquid to be transferred
    stocks_tuberack = labware where the stock solutions are stored
    stock_dict = dictionary containing the stock solutions
    out_well = well to be mixed
    vol_well = volume of the sample currently in the well
    '''

    pipette = pipettes[vol > pipettes[0].max_volume]
    if pipette.has_tip == False:
        pipette.pick_up_tip()
    peg_blowout_rate = 7
    default_pipette = pipette.flow_rate.blow_out
    pipette.flow_rate.blow_out = peg_blowout_rate

    height_disp = v2z_deepwellplate(vol_well)
    # height_disp = 2
    for vol_disp in control.sub_volumes(pipette, vol):
        #       pipette.aspirate(vol_disp,in_well.bottom(z=height_asp),rate=0.3)
        aspirate(protocol, pipette, vol_disp, rate=0.3, stocks_tuberack=stocks_tuberack, stock_dict_entry=stock_dict,
                 delay_seconds=7)

        protocol.delay(seconds=1)  # just checking

        pipette.dispense(vol_disp, out_well.bottom(z=height_disp), rate=0.2)
        protocol.delay(seconds=3)

        pipette.blow_out(out_well.bottom(z=v2z_deepwellplate(vol_well) + 2))

    pipette.flow_rate.blow_out = default_pipette
    # mix the sample
    if mix:
        # CriticalPoint_mix(protocol,pipettes,out_well,vol_well,N_mix=3)
        MixedPhase_mix(protocol, pipettes, out_well, vol_well)

    protocol.delay(seconds=1)  # just checking

    # make sure all the tips are dropped
    if pipette.has_tip == True and drop_tip:
        pipette.drop_tip()


def bsa_transfer_and_mix(protocol, pipettes, vol, stocks_tuberack, stock_dict, out_well, vol_well, mix=True,
                         drop_tip=True):
    '''
    Description: Transfers a volume of liquid with viscosity similar to BSA from the stock tube to the well and mixes the solution in the well.
    Input:
    pipettes = list of pipettes in the protocol
    vol = volume of the liquid to be transferred
    stocks_tuberack = labware where the stock solutions are stored
    stock_dict = dictionary containing the stock solutions
    out_well = well to be mixed
    vol_well = volume of the sample currently in the well
    '''
    pipette = pipettes[vol > pipettes[0].max_volume]
    if pipette.has_tip == False:
        pipette.pick_up_tip()

    bsa_blowout_rate = 30
    default_pipette = pipette.flow_rate.blow_out
    pipette.flow_rate.blow_out = bsa_blowout_rate

    for vol_disp in control.sub_volumes(pipette, vol):
        #        pipette.aspirate(vol_disp,in_well.bottom(z=3),rate=0.8)
        aspirate(protocol, pipette, vol_disp, rate=0.6, stocks_tuberack=stocks_tuberack, stock_dict_entry=stock_dict,
                 delay_seconds=2)
        protocol.delay(seconds=1)  # just checking

        height_disp = max(0.25 * v2z_deepwellplate(vol_well), 2)
        pipette.dispense(vol_disp, out_well.bottom(z=height_disp), rate=0.8)

        pipette.blow_out(out_well.bottom(z=v2z_deepwellplate(vol_well) + 2))

    pipette.flow_rate.blow_out = default_pipette  # set to default blow out rate

    if mix:
        thorough_mix(protocol, pipettes, out_well, vol_well)

    if pipette.has_tip == True and drop_tip:
        pipette.drop_tip()

def transfer_well_to_well(protocol, pipettes, vol, source_well, vol_src_well, dest_well, rate=1, mix=False,
                          drop_tip=True, sep=False):
    '''
    Function to transfer liquid from one well to another
    Input:
        pipettes = list of pipettes in the protocol
        vol = volume of the liquid to be transferred
        source_well = well from which the liquid is to be transferred
        vol_src_well = volume of the liquid currently in the source well
        height_asp = height from the bottom of the source well from which the liquid is to be aspirated
        dest_well = well to which the liquid is to be transferred
        rate = rate of aspiration and dispensing
        mix = whether to mix the liquid in the destination well
        drop_tip = whether to drop the tip after the transfer
        sep = whether to separate the dense and dilute phase mixtures
    '''
    pipette = pipettes[vol > pipettes[0].max_volume]
    if pipette.has_tip == False:
        pipette.pick_up_tip()

    bsa_blowout_rate = 30
    default_pipette = pipette.flow_rate.blow_out
    pipette.flow_rate.blow_out = bsa_blowout_rate

    if sep:
        for vol_disp in control.sub_volumes(pipette, vol):
            print(control.sub_volumes(pipette, vol))
            if vol_src_well + vol_disp > 2000:
                height_asp = v2z_deepwellplate(vol_src_well - 1.1 * vol_disp / 2)
                pipette.aspirate(vol_disp / 2, source_well.bottom(z=height_asp), rate=rate)
                protocol.delay(seconds=0.5)
                height_asp = v2z_deepwellplate(vol_src_well - 1.1 * vol_disp)
                pipette.aspirate(vol_disp / 2, source_well.bottom(z=height_asp), rate=rate)
                pipette.move_to(source_well.top(z=4))
                pipette.dispense(vol_disp, dest_well.bottom(z=2), rate=rate)
                pipette.blow_out(dest_well.bottom(z=v2z_deepwellplate(vol * 1.5)))
                vol_src_well -= vol_disp
            else:
                height_asp = v2z_deepwellplate(vol_src_well - 1.1 * vol_disp)
                pipette.aspirate(vol_disp, source_well.bottom(z=height_asp), rate=rate)
                pipette.move_to(source_well.top(z=4))
                pipette.dispense(vol_disp, dest_well.bottom(z=2), rate=rate)
                vol_src_well -= vol_disp
                pipette.blow_out(dest_well.bottom(z=v2z_deepwellplate(vol * 1.5)))
    else:
        height_asp = v2z_deepwellplate(vol_src_well - 1.1 * vol)
        pipette.aspirate(vol, source_well.bottom(z=height_asp), rate=rate)
        protocol.delay(seconds=0.5)
        pipette.move_to(source_well.top(z=4))
        pipette.dispense(vol, dest_well.bottom(z=2), rate=rate)
        pipette.blow_out(dest_well.bottom(z=v2z_deepwellplate(vol * 1.5)))

    if mix:
        thorough_mix(protocol, pipette, dest_well, vol)

    height_blow_out = v2z_deepwellplate(vol * 1.5)
    pipette.blow_out(dest_well.bottom(z=height_blow_out))
    pipette.flow_rate.blow_out = default_pipette  # set to default blow out rate

    if pipette.has_tip == True and drop_tip:
        pipette.drop_tip()


def transfer_well_to_trash(protocol, pipettes, vol, source_well, height_asp, stock_dict):
    '''
    Function to transfer liquid from one well to the trash
    This is mainly to get rid of any excess dilute phase liquid above the dense phase after separating the two phases
    Input:
        pipettes = list of pipettes in the protocol
        vol = volume of the liquid to be transferred
        source_well = well from which the liquid is to be transferred
        height_asp = height from the bottom of the source well from which the liquid is to be aspirated
        stock_dict = dictionary containing the location of the destination tube
    '''
    pipette = pipettes[vol > pipettes[0].max_volume]
    if pipette.has_tip == False:
        pipette.pick_up_tip()

    bsa_blowout_rate = 30
    default_pipette = pipette.flow_rate.blow_out
    pipette.flow_rate.blow_out = bsa_blowout_rate
    pipette.aspirate(vol, source_well.bottom(z=height_asp), rate=1)
    pipette.move_to(source_well.top(z=4))
    dest_tube = stock_dict['trash']['labware'][stock_dict['trash']['loc']]

    pipette.dispense(vol, dest_tube.bottom(z=falcon50ml_VtoH(stock_dict['trash']['vol'])), rate=1)
    pipette.blow_out(dest_tube.bottom(z=falcon50ml_VtoH(stock_dict['trash']['vol'])))
    pipette.flow_rate.blow_out = default_pipette  # set to default blow out rate
    pipette.drop_tip()


def add_mixing_vols(protocol, pipettes, plate, stock_dict, chemical_viscosity, unmixed_chem_list=[], sample_df=None):
    '''
    Function to add all the chemicals in the recipe to the sample
    '''

    output_dir = Path('/var/lib/jupyter/notebooks/')

    for index, row in sample_df.iterrows():
        sample_name = row['Name of sample list']
        sample_number = row['id']
        sample_well = row['Well']
        try:
            with open(output_dir / 'logfile.txt', 'a') as f:  # Open the file for writing
                f.write(f'Mixing sample - {sample_name}, sample number - {sample_number} in well - {sample_well} \n')
        except:
            pass

        for ind, chem in enumerate(row['order_mixing']):
            if chem not in unmixed_chem_list:
                if ind == 0:
                    basic_mix = False
                else:
                    basic_mix = True
                vol_well = row['Well Volume'] + row[chem]
                out_well = plate[row['Plate']][row['Well']]
                row['Well Volume'] = vol_well
                add_chemical(protocol, pipettes, stock_dict, out_well, chem, chemical_viscosity[chem], row[chem],
                             vol_well, basic_mix)
        sample_df.at[index, 'Well Volume'] = vol_well

    return sample_df


def add_unmixing_vols(protocol, pipettes, plate, stock_dict, unmixed_chem_list=[], sample_df=None):
    '''
    Function to add the volumes of the sample that do not require mixing to the wells
    '''

    for chem in unmixed_chem_list:
        chem_vols = np.column_stack((sample_df.index.values, sample_df[chem].values))
        mask = chem_vols[:, 1] > pipettes[0].max_volume
        for_pipette_2 = chem_vols[mask]
        for_pipette_1 = chem_vols[~mask]

        if for_pipette_1.size > 0:
            pipettes[0].pick_up_tip()
            for index, vol in for_pipette_1:
                out_well = plate[sample_df.loc[index]['Plate']][sample_df.loc[index]['Well']]
                add_chemical(protocol, pipettes, stock_dict, out_well, chem, 'water', vol,
                             sample_df.loc[index]['Well Volume'], basic_mix=False, drop_tip=False)
                sample_df.at[index, 'Well Volume'] = sample_df.loc[index]['Well Volume'] + vol
            pipettes[0].drop_tip()
        if for_pipette_2.size > 0:
            pipettes[1].pick_up_tip()
            for index, vol in for_pipette_2:
                out_well = plate[sample_df.loc[index]['Plate']][sample_df.loc[index]['Well']]
                add_chemical(protocol, pipettes, stock_dict, out_well, chem, 'water', vol,
                             vol_well=sample_df.loc[index]['Well Volume'], basic_mix=False, drop_tip=False)
                sample_df.at[index, 'Well Volume'] = sample_df.loc[index]['Well Volume'] + vol
            pipettes[1].drop_tip()

    return sample_df

# Just... don't use this.
def mix_Everything(protocol, pipettes, sample, out_well, stock_dict, chemical_viscosity):
    '''
    Function to add all the chemicals in the recipe to the sample
    '''
    vol_well = 0
    for ind, chem in enumerate(sample.order_mixing):
        # This if statement is confusing and needs to be fixed. - Takumi
        if ind == 0:
            basic_mix = False
        else:
            basic_mix = True
        vol_well += sample.recipe[chem]

        print('adding chemical {} and vol {}'.format(chem, sample.recipe[chem]))

        add_chemical(protocol, pipettes, stock_dict, out_well, chem, chemical_viscosity[chem], sample.recipe[chem],
                     vol_well, basic_mix)

# CriticalPoint_mix: NEVER use a capital letter for the first letter of a function name. - Takumi
def CriticalPoint_mix(protocol, pipettes, out_well, vol_well, N_mix=1):
    '''
    This is for mixing phase separating and one-phase samples close to the critical point
    The difference here is that both the solution added (PEG) and the solution in well (BSA+salts) are viscous and mixing of two viscous liquids takes more effort in mixing.
    Input:
    pipettes = list of pipettes in the protocol
    out_well = well to be mixed
    vol_well = volume of the sample currently in the well
    N_mix = number of times to mix the sample typically 3
    '''

    vol_mix = 0.2 * vol_well  # always mix 20% of the volume of solution in well
    vol_mix_static = 0.45 * vol_well
    if vol_mix < pipettes[0].min_volume:
        vol_mix = pipettes[0].min_volume

    for p in pipettes:
        if p.has_tip:
            if vol_mix <= p.max_volume and vol_mix >= p.min_volume:
                pipette = p
                break
            else:
                p.drop_tip()
        else:
            pipette = pipettes[vol_mix > pipettes[1].min_volume]

    if pipette.has_tip == False:
        pipette.pick_up_tip()
    #    N_mix_init = 20
    N_mix_init = 10
    N_mix_rot = 15
    N_mix_out = 5

    height_sample = v2z_deepwellplate(vol_well)
    h_mix_up = max(0.9 * height_sample, 2)
    h_mix_down = max(0.75 * height_sample, 2)

    for i in range(N_mix_init):
        pipette.aspirate(vol_mix, location=out_well.bottom(z=h_mix_up), rate=4)
        protocol.delay(seconds=0.8)
        pipette.dispense(vol_mix, location=out_well.bottom(z=h_mix_up), rate=10)

    h_mixrot_up_1 = max(0.75 * height_sample, 2)
    h_mixrot_up_2 = max(0.50 * height_sample, 2)
    h_mixrot_down_1 = max(0.25 * height_sample, 2)
    h_mixrot_down_2 = max(0.80 * height_sample, 2)

    n_mix = 0
    while n_mix < N_mix:

        for i in range(N_mix_rot):
            pipette.aspirate(vol_mix, location=out_well.bottom(z=h_mixrot_up_1), rate=5)
            protocol.delay(seconds=0.8)
            pipette.aspirate(vol_mix, location=out_well.bottom(z=h_mixrot_up_2), rate=5)
            protocol.delay(seconds=0.8)
            pipette.dispense(vol_mix, location=out_well.bottom(z=h_mixrot_down_1), rate=10)
            protocol.delay(seconds=0.8)
            pipette.dispense(vol_mix, location=out_well.bottom(z=h_mixrot_down_2), rate=10)
            protocol.delay(seconds=1)

        for i in range(N_mix_init):
            pipette.aspirate(vol_mix_static, location=out_well.bottom(z=h_mix_up), rate=4)
            protocol.delay(seconds=0.3)
            h_mix_down = max(np.random.uniform(0.25, 0.9) * height_sample, 2)
            pipette.dispense(vol_mix_static, location=out_well.bottom(z=h_mix_down), rate=10)
        n_mix += 1

    peg_blowout_rate = 7
    default_pipette = pipette.flow_rate.blow_out
    pipette.flow_rate.blow_out = peg_blowout_rate
    pipette.blow_out(out_well.bottom(z=h_mix_up))
    pipette.flow_rate.blow_out = default_pipette  # set to default blow out rate

    pipette.drop_tip()

# MixedPhase_mix: NEVER use a capital letter for the first letter of a function name. - Takumi
def MixedPhase_mix(protocol, pipettes, out_well, vol_well):
    '''
    For mixing samples which contain PEG - mixed phase samples

    Input:
    pipettes = list of pipettes in the protocol
    out_well = well to be mixed
    vol_well = volume of the sample currently in the well
    '''

    vol_mix = 0.2 * vol_well  # always mix 20% of the volume of solution in well

    if vol_mix < pipettes[0].min_volume:
        vol_mix = pipettes[0].min_volume

    for p in pipettes:
        if p.has_tip:
            if 2 * vol_mix <= p.max_volume and 2 * vol_mix >= p.min_volume:
                pipette = p
                break
            else:
                p.drop_tip()
        else:
            pipette = pipettes[vol_mix > pipettes[1].min_volume]

    if pipette.has_tip == False:
        pipette.pick_up_tip()
    #    N_mix_init = 20
    N_mix_init = 10
    N_mix_rot = 15
    N_mix_out = 5

    height_sample = v2z_deepwellplate(vol_well)
    h_mix_up = max(0.85 * height_sample, 2)
    h_mix_down = max(0.75 * height_sample, 2)

    for i in range(N_mix_init):
        pipette.aspirate(vol_mix, location=out_well.bottom(z=h_mix_up), rate=2)
        protocol.delay(seconds=0.8)
        pipette.dispense(vol_mix, location=out_well.bottom(z=h_mix_down), rate=10)

    h_mixrot_up_1 = max(0.75 * height_sample, 2)
    h_mixrot_up_2 = max(0.50 * height_sample, 2)
    h_mixrot_down_1 = max(0.25 * height_sample, 2)
    h_mixrot_down_2 = max(0.80 * height_sample, 2)

    print(pipette, vol_mix)
    for i in range(N_mix_rot):
        pipette.aspirate(vol_mix, location=out_well.bottom(z=h_mixrot_up_1), rate=2)
        protocol.delay(seconds=0.8)
        pipette.aspirate(vol_mix, location=out_well.bottom(z=h_mixrot_up_2), rate=2)
        protocol.delay(seconds=0.8)
        pipette.dispense(vol_mix, location=out_well.bottom(z=h_mixrot_down_1), rate=10)
        protocol.delay(seconds=0.8)
        pipette.dispense(vol_mix, location=out_well.bottom(z=h_mixrot_down_2), rate=10)
        protocol.delay(seconds=1)

    peg_blowout_rate = 7
    default_pipette = pipette.flow_rate.blow_out
    pipette.flow_rate.blow_out = peg_blowout_rate
    pipette.blow_out(out_well.bottom(z=h_mix_up))
    pipette.flow_rate.blow_out = default_pipette  # set to default blow out rate

    pipette.drop_tip()

def water_transfer_and_mix(protocol, pipettes, vol, vol_well, stocks_tuberack, stock_dict, out_well, basic_mix=False,
                           thorough_mix=False, drop_tip=True):
    '''
    Parameters: pipette,vol,in_well,height_asp,out_well,height_disp,n_mix,vol_mix,height_mix,mix=False
    Description: Transfers a volume of liquid with viscosity similar to water from the stock tube to the well
    '''

    pipette = pipettes[vol > pipettes[0].max_volume]
    if pipette.has_tip == False:
        pipette.pick_up_tip()

    water_blowout_rate = 50
    default_pipette = pipette.flow_rate.blow_out
    pipette.flow_rate.blow_out = water_blowout_rate

    if basic_mix:
        height_disp = v2z_deepwellplate(vol_well * 1.5)
    else:
        height_disp = 2

    for vol_disp in control.sub_volumes(pipette, vol):
        # aspirate solution from stock tube
        aspirate(protocol, pipette, vol_disp, rate=1, stocks_tuberack=stocks_tuberack, stock_dict_entry=stock_dict)

        pipette.dispense(vol_disp, out_well.bottom(z=height_disp))

        pipette.blow_out(out_well.bottom(z=v2z_deepwellplate(vol_well * 1.5)))

    pipette.flow_rate.blow_out = default_pipette  # set to default blow out rate

    if thorough_mix:
        thorough_mix(protocol, pipette, out_well, vol_well)
    elif basic_mix:
        n_mix = 2
        vol_mix = max(0.2 * vol_well, pipette.min_volume)
        pipette.mix(repetitions=n_mix, volume=vol_mix, location=out_well.bottom(z=2), rate=5)
        pipette.blow_out(out_well.bottom(z=20))

    if pipette.has_tip == True and drop_tip:
        pipette.drop_tip()

# FUNCTIONS THAT ARE TOO SPECIFIC
def create_labware_file_for_dilutions(nPlates, pipetteTypes, tip_boxes, stock_dict):
    '''
    Writes a csv file that contains the run details for the OT2.
    The csv file includes all sorts of details such as the pipetteTypes, tip boxes, stock concentrations, etc. that have to be specified before the run.
    '''

    tip_names = {'p20': 'opentrons_96_tiprack_20ul', 'p300': 'opentrons_96_tiprack_300ul',
                 'p1000': 'opentrons_96_tiprack_1000ul'}
    pipette_pos = ['left', 'right']
    dec_loc = 11
    falcon15Wells = ['A1', 'A2', 'B1', 'B2', 'C1', 'C2']
    falcon50Wells = ['A3', 'A4', 'B3', 'B4']
    with open('dilution_labware.csv', 'w', newline='') as f:
        fwriter = csv.writer(f, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        fwriter.writerow(['Sample_plate', 'plate_name', 'deck_location'])
        dec_loc -= 1
        fwriter.writerow(['plate', 'thermoscientificnunc_96_wellplate_2000ul', str(dec_loc)])
        fwriter.writerow([])
        fwriter.writerow(['Tranfer_plate', 'plate_name', 'deck_location'])
        dec_loc -= 1
        fwriter.writerow(['plate', 'thermoscientificnunc_96_wellplate_2000ul', str(dec_loc)])
        fwriter.writerow([])
        fwriter.writerow(['DilutionPlates_list', 'plate_name', 'deck_location'])
        for i in range(nPlates):
            dec_loc -= 1
            fwriter.writerow(['plate', 'thermoscientificnunc_96_wellplate_2000ul', str(dec_loc)])
        fwriter.writerow([])

        fwriter.writerow(['pipette_list', 'name', 'location (right/left)'])
        for ind, pipette in enumerate(pipetteTypes):
            if pipette == 'p20':
                fwriter.writerow(['pipette', 'p20_single_gen2', pipette_pos[ind]])
            elif pipette == 'p300':
                fwriter.writerow(['pipette', 'p300_single_gen2', pipette_pos[ind]])
            elif pipette == 'p1000':
                fwriter.writerow(['pipette', 'p1000_single_gen2', pipette_pos[ind]])
            fwriter.writerow(['tip box', 'name', 'deck location', 'starting tip position'])

            for i in range(int(tip_boxes[ind])):
                dec_loc -= 1
                fwriter.writerow(['tiprack', tip_names[pipette], str(dec_loc), 'A1'])
            fwriter.writerow([])

        fwriter.writerow([])
        fwriter.writerow(
            ['stock_list', 'chemical name', 'concentration', 'volume (mL)', 'labware_name', 'deck_location',
             'location_on_labware', 'tube_type'])

        vol_water = stock_dict['water'] / 1000
        while vol_water > 0:
            dec_loc -= 1
            loc = falcon50Wells.pop(0)
            fwriter.writerow(
                ['stock', 'water', 0, 45, 'opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical', dec_loc, loc,
                 'falcon50'])
            vol_water -= 45
        dec_loc -= 1
        loc = falcon50Wells.pop(0)
        fwriter.writerow(
            ['stock', 'trash', 0, 10, 'opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical', dec_loc, loc, 'falcon50'])

def create_run_file_for_dilutions(sample_list_nos=[], vol_sample_per_dil=20, dilution_factor=10, Num_dilutions=1,
                                  stock_dict={}, vol_fraction_dense=0.2):
    all_sample_lists = sorted([f for f in os.listdir(master_dir / 'prepared_samples/') if f.startswith('sample_list')])
    dilution_samples = []
    vols = []
    sample_number = 0
    data_list = []

    sample_lists = [sorted([f for f in os.listdir(master_dir / 'prepared_samples/') if f.startswith('sample_list')])[i]
                    for i in sample_list_nos]

    filename = 'dilutions_run_file-' + datetime.datetime.now().strftime("%Y-%m-%d__%H-%M-%S") + '.pkl'
    with open(master_dir / filename, 'wb') as f:
        pickle.dump([vol_sample_per_dil, Num_dilutions, dilution_factor, sample_lists], f)

    vol_water = 0

    for i in sample_list_nos:
        with open(master_dir / 'prepared_samples' / all_sample_lists[i], 'rb') as f:
            metadata, chemical_viscosity, samples = pickle.load(f)

        for ind, sample in enumerate(samples):
            sample.id = ind + 1  # remove this later
            for i in range(Num_dilutions):
                sample_number += 1
                vols.append(vol_sample_per_dil)
                vols.append(vol_sample_per_dil * (dilution_factor - 1))
                vol_water += vol_sample_per_dil * (dilution_factor - 1)

            vol_dense_phase = (vol_fraction_dense) * sample.volume
            vol_dilute_phase = sample.volume - vol_dense_phase
            transfer_vol = 0.6 * vol_dilute_phase
            vols.append(transfer_vol)
    stock_dict['water'] = vol_water

    if any(vol < 0 for vol in vols):
        print('error: negative volume')
        return
    if any(vol > 0 and vol < 20 for vol in vols):
        p20_req = 1
        tips_20_1000 = estimate_required_tips([0, 20], [100, 1000], vols)
        tips_20_300 = estimate_required_tips([0, 20], [20, 300], vols)
        if tips_20_1000[0] + tips_20_1000[1] < tips_20_300[0] + tips_20_300[1]:
            pipetteTypes = ['p20', 'p1000']
            tip_boxes = [np.ceil(tips_20_1000[0] / 96), np.ceil(tips_20_1000[1] / 96)]
            Num_tips = [tips_20_1000[0], tips_20_1000[1]]
        else:
            pipetteTypes = ['p20', 'p300']
            tip_boxes = [np.ceil(tips_20_300[0] / 96), np.ceil(tips_20_300[1] / 96)]
            Num_tips = [tips_20_300[0], tips_20_300[1]]
    else:
        tips = [estimate_required_tips([0, 20], [20, 300], vols), estimate_required_tips([0, 20], [100, 1000], vols),
                estimate_required_tips([20, 300], [300, 1000], vols)]
        total_tips = np.array([tips[0][0] + tips[0][1], tips[1][0] + tips[1][1], tips[2][0] + tips[2][1]])
        min_tip = np.argmin(np.flip(total_tips))
        if min_tip == 2:
            pipetteTypes = ['p20', 'p300']
            tip_boxes = [np.ceil(tips[0][0] / 96), np.ceil(tips[0][1] / 96)]
            Num_tips = [tips[0][0], tips[0][1]]
        elif min_tip == 2:
            pipetteTypes = ['p20', 'p1000']
            tip_boxes = [np.ceil(tips[1][0] / 96), np.ceil(tips[1][1] / 96)]
            Num_tips = [tips[1][0], tips[1][1]]
        else:
            pipetteTypes = ['p300', 'p1000']
            tip_boxes = [np.ceil(tips[2][0] / 96), np.ceil(tips[2][1] / 96)]
            Num_tips = [tips[2][0], tips[2][1]]

    print('summary of the minimum requirements:')
    print('pipetteTypes: {}'.format(pipetteTypes))
    print('tip boxes: {}'.format(tip_boxes))
    print('Number of tips used: {}'.format(Num_tips))

    nPlates_dils = int(np.ceil(sample_number / 96))
    create_labware_file_for_dilutions(nPlates_dils, pipetteTypes, tip_boxes, stock_dict)

    print('run file saved as {}'.format(filename))




# Good ideas, bad implementation. It is too hard to read. - Takumi
## Calling an arg 'purpose' is not a good idea. It is too subjective.
def thorough_mix(protocol, pipettes, out_well, vol_well, purpose='SamplePrep'):
    '''
    This is for thorough mixing of the sample in the well - Not to be used for phase separated samples
    Typeically used for mixing the BSA and other salts in the sample just before the PEG is added to it.

    Parameters:
    pipette: pipette to be used for mixing
    out_well: well to be mixed
    vol_well: volume of the sample currently in the well
    purpose: purpose of the mixing - 'SamplePrep' or 'Dilutions' or 'CriticalPt' - I use different mixing parameters for these two purposes
            for the dilutions, the number of mixing steps is greater and use a larger volume for mixing
    '''

    if purpose == 'SamplePrep':
        vol_mix = 0.20 * vol_well  # always mix 20% of the volume of solution in well
        n_mix_down = 10
        n_mix_up = 10
    elif purpose == 'Dilutions':
        vol_mix = 0.50 * vol_well  # always mix 50% of the volume of solution in well
        n_mix_down = 15
        n_mix_up = 15
    elif purpose == 'CriticalPt':
        vol_mix = 0.20 * vol_well
        n_mix_down = 15
        n_mix_up = 15

    # vol_mix = 200    # for the CP samples
    if vol_mix < pipettes[0].min_volume:
        vol_mix = pipettes[0].min_volume

    for p in pipettes:
        if p.has_tip:
            if vol_mix <= p.max_volume and vol_mix >= p.min_volume:
                pipette = p
                break
            else:
                p.drop_tip()
        else:
            pipette = pipettes[vol_mix > pipettes[1].min_volume]

    if pipette.has_tip == False:
        pipette.pick_up_tip()

    height_sample = v2z_deepwellplate(vol_well)

    h_mix_down_1 = max(0.50 * height_sample, 2)
    h_mix_down_2 = max(0.35 * height_sample, 2)

    for i in range(n_mix_down):
        pipette.aspirate(vol_mix, location=out_well.bottom(z=h_mix_down_1), rate=3)
        pipette.dispense(vol_mix, location=out_well.bottom(z=h_mix_down_2), rate=5)

    h_mix_up_1 = max(0.50 * height_sample, 2)
    h_mix_up_2 = max(1 * height_sample, 2)
    for i in range(n_mix_up):
        pipette.aspirate(vol_mix, location=out_well.bottom(z=h_mix_up_1), rate=3)
        pipette.dispense(vol_mix, location=out_well.bottom(z=h_mix_up_2), rate=5)

    # set blow out rates and blow out
    water_blowout_rate = 30
    default_pipette = pipette.flow_rate.blow_out
    pipette.flow_rate.blow_out = water_blowout_rate
    pipette.blow_out(out_well.bottom(z=h_mix_up_2 + 2))
    # set blow out rates back to default:
    pipette.flow_rate.blow_out = default_pipette  # set to default blow out rate
    pipette.move_to(out_well.top(z=4))

    pipette.drop_tip()


def aspirate(protocol, pipette, vol, rate, stocks_tuberack, stock_dict_entry, delay_seconds=0.5):
    '''
    Parameters: pipette,vol,rate,stocks_tuberack,stock_dict_entry for the stock solution to be added.
    Description: Aspirates a volume of liquid from the stock solution - when stock is in falcon tubes
    '''
    height_asp = 2
    if stock_dict_entry['Tube'] == 'falcon15':
        height_asp = falcon15ml_VtoH(stock_dict_entry['vol'] - vol / 1000)
    elif stock_dict_entry['Tube'] == 'falcon50':
        height_asp = falcon50ml_VtoH(stock_dict_entry['vol'] - vol / 1000)
    stock_dict_entry['vol'] -= vol / 1000
    pipette.aspirate(vol, stocks_tuberack[stock_dict_entry['loc']].bottom(z=height_asp), rate=rate)
    protocol.delay(seconds=delay_seconds)
    pipette.move_to(stocks_tuberack[stock_dict_entry['loc']].bottom(z=height_asp + 10), speed=15)



# FUNCTIONS THAT ARE NOT NECESSARY
## list_labware- Who needs this? You already have a csv. - Takumi
def list_labware(filepath=master_dir / 'labware_configuration.csv'):
    """
    Function to list the labware configuration from the csv file

    Parameters
    ----------
    filepath: path to the csv file containing the labware configuration

    Returns
    -------
    pipette_names: list of pipette names
    plate_names: list of plate names
    starting_well: starting well for the samples
    starting_tips_1: starting tip for the first pipette
    starting_tips_2: starting tip for the second pipette
    stocks: dictionary containing the stock solution location,
    plate_names, starting_well, starting_tips_1, starting_tips_2, stocks
    """
    # AVAIlABLE PIPEETTES AND TIPRACKS
    pipette_list = ['p20_single_gen2', 'p300_single_gen2', 'p1000_single_gen2']
    tiprack_list = ['opentrons_96_tiprack_20ul', 'opentrons_96_tiprack_300ul', 'opentrons_96_tiprack_1000ul']

    # INITIALIZAITON
    ## PIPEETTES
    pipette_cnt = -1
    pipette_names = []
    # pipettes = []

    # TIPRACKS
    tipracks_1 = []
    tipracks_2 = []
    starting_tips_1 = []
    starting_tips_2 = []
    tiprack_names = []

    # WELL PLATES
    plate_names = []

    # STOCKS
    stocks = {}
    stock_tube_racks = []
    stock_names = []
    stock_index = -1

    with open(filepath, newline='') as csvfile:
        labware_reader = csv.reader(csvfile, delimiter=',', quotechar='|')

        for row in labware_reader:
            if any(x.strip() for x in row):
                if row[0] == "pipette":
                    pipette_cnt += 1
                    name_pipette = row[1]
                    loc_pipette = row[2]
                    pipette_names.append(name_pipette)
                    print(pipette_cnt)
                if row[0] == "tiprack":
                    name_tiprack = row[1]
                    tiprack_deck = int(row[2])
                    starting_tip = row[3]
                    tiprack_names.append(name_tiprack)
                    # There is a bug that needs to be fixed.
                    # starting_tips_1 is associated with the LEFT TIP.
                    # starting_tips_2 is associated with the RIGHT TIP.
                    if name_tiprack in tiprack_list and tiprack_deck in np.arange(1, 12):
                        if pipette_cnt == 1:
                            starting_tips_2.append(starting_tip)
                        else:
                            starting_tips_1.append(starting_tip)
                    else:
                        raise ValueError('Tiprack not specified correctly')
                if row[0] == "plate":
                    name_plate = row[1]
                    plate_deck = int(row[2])
                    starting_well = row[3]
                    plate_names.append(name_plate)
                if row[0] == "stock":
                    key = row[1]
                    conc = float(row[2])
                    vol = float(row[3])
                    labware_name = row[4]
                    stock_deck = int(row[5])
                    loc = row[6]
                    tube = row[7]

                    stocks[key] = {}
                    stocks[key]['conc'] = conc
                    stocks[key]['vol'] = vol
                    stocks[key]['deck'] = stock_deck
                    stocks[key]['loc'] = loc
                    stocks[key]['Tube'] = tube
        labware_info = pipette_names, plate_names, starting_well, starting_tips_1, starting_tips_2, stocks
        description = ['Pipette names', 'Plate names', 'Starting pos of well plate',
                       f'Starting Pos of Tip rack {tiprack_names[0]}',
                       f'Starting Pos of Tip rack {tiprack_names[1] if len(tiprack_names) > 1 else "None"}',
                       'Stock solutions']

        for desc, info in zip(description, labware_info):
            print(desc + ':', info)

        return labware_info

# 1. I would not use these functions. Messing with a file architecture hinders readability and reproducibility. - Takumi
# 2. Confusing names! - TM
## Do not use "add" for two different operations.
## Currently, you "add" to the sample list and "add" (chemicals) to the prepared samples.
def add_to_prepared_samples(sample_list_nos):
    '''
    This function adds the sample lists which are prepared by otto to the prepared_samples directory.
    Inputs:
        sample_list_nos: list of sample list numbers that need to be added to the prepared_samples directory
    Function:
    Sample lists corresponding to the input will be deleted from the samples folder and put in the prepared_samples folder.
    '''

    if not os.path.exists(master_dir / 'prepared_samples'):
        os.makedirs(master_dir / 'prepared_samples')

    sample_lists = sorted([f for f in os.listdir(sample_dir) if f.startswith('sample_list')])
    for ind in sample_list_nos:
        file = sample_lists[ind]
        shutil.move(sample_dir / file, master_dir / 'prepared_samples' / file)
        print('moved {}'.format(file))

# This function is okay but not essential. - Takumi
## Why would you want to remove a sample list this way while you ask users to edit a csv file directly?
def view_sample_lists(sample_dir=master_dir / 'samples', show_recipe=False):
    """
    Returns a list of Sample objects saved in the samples directory

    Parameters
    ----------
    sample_dir: str, directory where the sample lists are saved, default is 'samples'
    show_recipe: bool, whether to show the recipe for the samples or not, default is False

    Returns
    -------
    samples: list of Sample objects
    """
    sample_list = []
    cnt = 0
    for file in sorted(os.listdir(dir)):
        if file.startswith('sample_list'):
            with open(dir / file, 'rb') as f:
                metadata, chemical_viscosity, samples = pickle.load(f) # This is abnocious- TM
                print('Sample list {}:'.format(cnt))
                print('Name: {}'.format(metadata[0]))
                print('Time when this object was created: {}'.format(metadata[1]))
                print('Comments: {}'.format(metadata[2]))
                for sample in samples:
                    print(sample.composition)
                    if show_recipe:
                        print(sample.recipe)
                print('Order of mixing: {}'.format(samples[0].order_mixing))
                # print('recipe for sample 0: {}'.format(samples[0].recipe))
                print('\n\n')
                cnt += 1
    return samples
def remove_sample_list(sample_list_nos=None, sample_dir=sample_dir):
    '''
    This function removes the sample lists that are no longer required.
    Input:
        sample_list_nos: list of sample list numbers that need to be removed
        dir: directory where the sample lists are saved - by default this is going to be master_dir/samples but can be changed if the sample lists are saved in a different directory.
    '''
    if os.path.exists(sample_dir):
        sample_lists = sorted([f for f in os.listdir(sample_dir) if f.startswith('sample_list')])
        if sample_list_nos is None:
            sample_list_nos = list(range(len(sample_lists)))
        if len(sample_lists) > 0:
            for ind in sample_list_nos:
                file = sample_lists[ind]
                os.remove(sample_dir / file)
                print('removed {}'.format(file))


def estimate_requirements_DEPRICATED(volume, sample_list_nos, stock_dict, show=True,
                          file_header='sample_list', solvent='water'):
    """
    This function estimates the requirements for the run based on the sample list and stock concentrations.
    The requirements include the number of plates, pipettes, tip boxes, and stock volumes required for the run.

    Parameters
    ----------
    volume: int, volume of the sample in ul
    sample_list_nos: list, list of sample list numbers that need to be considered for the run
    stock_dict: dict, dictionary containing the stock concentrations of the chemicals
    show: bool, whether to show the requirements or not, default is True
    file_header: str, header of the sample list files, default is 'sample_list'
    solvent: str, solvent used in the samples, default is 'water'

    Returns
    -------
    nPlates: int, number of plates required for the run
    pipetteTypes: list, list of pipetteTypes required for the run
    tip_boxes: list, list of tip boxes required for the run
    stock_vols: dict, dictionary containing the stock volumes required for the run
    """
    print('-' * 40)
    print('ESTIMATING REQUIREMENTS FOR THE RUN')
    p20_req = 0
    tips_required = []
    all_sample_lists = sorted([f for f in os.listdir(sample_dir) if f.startswith(file_header)])
    vols = []
    n_samples = 0
    stock_vols = dict.fromkeys(stock_dict, 0)
    for ind in sample_list_nos:
        with open(sample_dir / all_sample_lists[ind], 'rb') as f:
            metadata, chemical_viscosity, samples = pickle.load(f)
            for sample in samples:
                n_samples += 1
                recipe = get_recipe(volume, sample, stock_dict, solvent=solvent)
                vols += list(recipe.values())
                print(recipe, solvent)
                for key in recipe:
                    stock_vols[key] += recipe[key]
                sample.recipe = recipe
                sample.volume = volume
                sample.stock_concs = stock_dict
                if any(vol < 0 for vol in list(recipe.values())):
                    print('error: negative volume for sample {}'.format(sample.composition))
                    print('sample recipe: {}'.format(sample.recipe))
                    return

        with open(sample_dir / all_sample_lists[ind], 'wb') as f:
            pickle.dump([metadata, chemical_viscosity, samples], f)

    nPlates = np.ceil(n_samples / 96)

    if any(vol < 0 for vol in vols):
        print('error: negative volume')

        return

    if any(vol > 0 and vol < 20 for vol in vols):
        p20_req = 1

    if p20_req:
        tips_20_1000 = estimate_required_tips([0, 20], [100, 1000], vols)
        tips_20_300 = estimate_required_tips([0, 20], [20, 300], vols)
        if tips_20_1000[0] + tips_20_1000[1] < tips_20_300[0] + tips_20_300[1]:
            pipetteTypes = ['p20', 'p1000']
            tip_boxes = [np.ceil(tips_20_1000[0] / 96), np.ceil(tips_20_1000[1] / 96)]
            Num_tips = [tips_20_1000[0], tips_20_1000[1]]
        else:
            pipetteTypes = ['p20', 'p300']
            tip_boxes = [np.ceil(tips_20_300[0] / 96), np.ceil(tips_20_300[1] / 96)]
            Num_tips = [tips_20_300[0], tips_20_300[1]]
    else:
        tips = [estimate_required_tips([0, 20], [20, 300], vols), estimate_required_tips([0, 20], [100, 1000], vols),
                estimate_required_tips([20, 300], [300, 1000], vols)]
        total_tips = np.array([tips[0][0] + tips[0][1], tips[1][0] + tips[1][1], tips[2][0] + tips[2][1]])
        min_tip = np.argmin(total_tips)
        if min_tip == 0:
            pipetteTypes = ['p20', 'p300']
            tip_boxes = [np.ceil(tips[0][0] / 96), np.ceil(tips[0][1] / 96)]
            Num_tips = [tips[0][0], tips[0][1]]
        elif min_tip == 1:
            pipetteTypes = ['p20', 'p1000']
            tip_boxes = [np.ceil(tips[1][0] / 96), np.ceil(tips[1][1] / 96)]
            Num_tips = [tips[1][0], tips[1][1]]
        else:
            pipetteTypes = ['p300', 'p1000']
            tip_boxes = [np.ceil(tips[2][0] / 96), np.ceil(tips[2][1] / 96)]
            Num_tips = [tips[2][0], tips[2][1]]
    estimates = int(nPlates), pipetteTypes, tip_boxes, stock_vols
    if show:
        print('summary of the minimum requirements:')
        print('number of plates: {}'.format(nPlates))
        print('pipetteTypes: {}'.format(pipetteTypes))
        print('tip boxes: {}'.format(tip_boxes))
        print('Number of tips used: {}'.format(Num_tips))
        print('minimum stock volumes required:')
        for key in stock_vols:
            print('{}: {} mL'.format(key, stock_vols[key] / 1000))
    nPlates = int(nPlates)
    print('-' * 40)
    return nPlates, pipetteTypes, tip_boxes, stock_vols

def estimate_max_concentration(stock_dict,chemical_list,composition,fixed_chems):
    '''
    This function calculates the maximum concentrations of the chemicals in the sample given the stock concentrations and the composition of the sample.
    Inputs:
        stock_dict: dictionary containing the stock concentrations of the chemicals
        chemical_list: list of chemicals in the sample
        composition: dictionary containing the composition of the sample
        fixed_chems: list of chemicals that are fixed in the sample and cannot be varied
    Output:
        max_concs: dictionary containing the maximum concentrations of the chemicals in the sample
    '''
    max_concs = {}
    V = 1000
    vol = 0
    for chem in fixed_chems:
        vol += composition[chem]*V/stock_dict[chem]
    remainder_chems = [chem for chem in chemical_list if chem not in fixed_chems]
    for chem in remainder_chems:
        max_concs[chem] = (V - vol)*stock_dict[chem]/V
    return max_concs
