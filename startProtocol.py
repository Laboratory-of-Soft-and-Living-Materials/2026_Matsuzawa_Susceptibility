"""
Author: Takumi Matsuzawa
Date: 2026/04/27

Code for sample preparation using OT-2 (Opentrons Labworks Inc.)
"""
import copy
import glob
import os
import sys
# Add the path to the sys.path
sys.path.append(os.getcwd())

import pandas as pd
import sample_handling as sh # sample preparation
import rotbot_control as control # mechanical control
import pickle
import datetime
import shutil
import subprocess
import platform
import pathlib
from pathlib import Path
from opentrons import protocol_api

# Necessary for Path library to work on Windows and Linux
plt = platform.system()
if plt == 'Linux': pathlib.WindowsPath = pathlib.PosixPath

metadata = {
    'apiLevel': '2.13',
    'protocolName': 'General sample prep',
    'description': 'A script to prepare samples using OT-2',
    'author': 'Takumi Matsuzawa',
}

master_dir = sh.master_dir # Path to the master directory (Pathlib Path)
sample_dir = sh.sample_dir # Path to the samples directory (Pathlib Path)
out_dir = sh.out_dir # Path to the output directory (Pathlib Path)

# DEFAULT TRANSFER SETTINGS
# Thin fluids (e.g. water)
transfer_thin_kwargs = {'flow_rate': 1.0, 'aspirate_wait_time': 0.5, 'hang_time': 0.5,
                        'blowout_rate': None, 'max_blowout_height': 200,
                        'v_src': 0, 'v_dst': 0, 'v2z_src': None, 'v2z_dst': None, 'z_src': None, 'z_dst': None,
                        'z_offset_src': 0, 'z_offset_dst': 0, 'drop_tip': False,
                        'default_z': 5, 'mix': False,
                        'mix_kwargs': {'repetitions': 4, 'volume': 20, 'location': None, 'rate': 5}
                        }
# Thick fluids (e.g. PEG, glycerol, etc.)
transfer_thick_kwargs = {'flow_rate': 0.3, 'aspirate_wait_time': 7, 'hang_time': 8,
                         'blowout_rate': 7, 'max_blowout_height': 200,
                         'v_src': 0, 'v_dst': 0, 'v2z_src': None, 'v2z_dst': None, 'z_src': None, 'z_dst': None,
                         'z_offset_src': 0, 'z_offset_dst': 0, 'drop_tip': False,
                         'default_z': 5, 'mix': True,
                         'mix_kwargs': {'repetitions': 4, 'volume': 20, 'location': None, 'rate': 5}}
# SPECIAL INSTRUCTIONS PER CHEMICAL
# BSA
transfer_BSA_kwargs =  {'flow_rate': 0.6, 'aspirate_wait_time': 2, 'hang_time': 5,
                        'blowout_rate': 30, 'max_blowout_height': 200,
                        'v_src': 0, 'v_dst': 0, 'v2z_src': None, 'v2z_dst': None, 'z_src': None, 'z_dst': None,
                        'z_offset_src': 0, 'z_offset_dst': 0,  'drop_tip': False,
                        'default_z': 5, 'mix': True,
                        'mix_kwargs': {'repetitions': 4, 'volume': 20, 'location': None, 'rate': 5}}
# FUS (z_offset_src is NON-ZERO because one must use 300uL tip racks to hold Eppendorf tubes)
transfer_Bik1_kwargs = transfer_FUS_kwargs =  {'flow_rate': 1.0, 'aspirate_wait_time': 2, 'hang_time': 5,
                        'blowout_rate': None, 'max_blowout_height': 200,
                        'v_src': 0, 'v_dst': 0, 'v2z_src': None, 'v2z_dst': None, 'z_src': None, 'z_dst': None,
                        'z_offset_src': 41.32, 'z_offset_dst': 0, 'drop_tip': False,
                        'default_z': 100, 'mix': False,
                        'mix_kwargs': {'repetitions': 4, 'volume': 20, 'location': None, 'rate': 5}
                        }
transfer_default_kwargs = transfer_thin_kwargs

def find_samples(sample_dir=sample_dir):
    """
    Find the samples in the sample directory
    Parameters
    ----------
    sample_dir : pathlib.Path
        Path to the sample directory

    Returns
    -------
    sample_paths : list
        List of paths to the samples
    """
    sample_paths = sorted(glob.glob(str(sample_dir / '*.pkl')))
    return sample_paths

# HELPER FUNCTIONS
def get_transfer_kwargs(rel_viscosity):
    """
    Returns the transfer kwargs for the chemical based on the relative viscosity to water
    
    Parameters
    ----------
    rel_viscosity: str
        Relative viscosity of the chemical (thick, thin, name of the chemical (e.g. BSA))

    Returns
    -------
    kwargs : dict
        Transfer kwargs for the chemical
    """
    print('...... transfer mode:', rel_viscosity)
    if rel_viscosity.lower() == 'thick':
        kwargs = transfer_thick_kwargs
    elif rel_viscosity.lower() == 'thin':
        kwargs = transfer_thin_kwargs
    elif rel_viscosity.lower() == 'bsa':
        kwargs = transfer_BSA_kwargs
    elif rel_viscosity.lower() == 'fus':
        kwargs = transfer_FUS_kwargs
    else:
        print(f'... No settings found to transfer {rel_viscosity}. Using default transfer kwargs')
        kwargs = transfer_default_kwargs
    return kwargs

def get_stock_locator(tuberacks, df_stocks):
    """
    Returns the location of the stock tube
    Parameters
    ----------
    tuberacks : list, List of tuberacks
    df_stocks : pd.DataFrame, Dataframe of the stocks

    Returns
    -------
    stock_location : list
        Location of the stock tube
    """
    stock_locator = {}
    chemicals = df_stocks.index.values

    for chem in chemicals:
        loc = df_stocks.loc[chem, 'rack_location']
        for i, tuberack in enumerate(tuberacks):
            tuberack_name = tuberack.name
            if tuberack_name == df_stocks.loc[chem, 'labware_name']:
                if isinstance(loc, str):
                    stock_locator[chem] = tuberack[loc]
                elif isinstance(loc, list):
                    stock_locator[chem] = [tuberack[wellname] for wellname in loc]
                else:
                    print(f'... Incorrect location format for {chem}. '
                          f'Please check the labware configuration file.')
    return stock_locator

def get_deck2plate(well_plates):
    """
    Returns a dictionary that maps the plate object from e deck location
    Parameters
    ----------
    plates: list, List of well plate objects
    deck_loc: int, Deck location, e.g. 10, 11

    Returns
    -------
    deck2plate: dict, Dictionary that maps the deck location to the plate object
    """
    deck_loc_registered = [int(plate.parent) for plate in well_plates]
    deck2plate = dict(zip(deck_loc_registered, well_plates))
    return deck2plate

# HELPER FUNCTIONS FOR PREPARING SAMPLES
# PREPARE SAMPLES
def prepare_sample(protocol, labware, dfs_labware, sample):
    """
    Prepare the sample according to the recipe
    Parameters
    ----------
    protocol : protocol_api.ProtocolContext
        The protocol context
    sample : sh.Sample
        The sample to be prepared
    dfs_labware : dict
        Dictionary of labware dataframes

    Returns
    -------
    df_stocks : pd.DataFrame, Updated stock dataframe
    """
    def max_depth(lst, depth=1):
        """
        Returns the maximum depth of the nested list

        Parameters
        ----------
        lst: list, The list to be analyzed
        depth: int, The depth level of the list if it is not a nested list. ['a', 'b', 'c'] has a depth of 1.

        Returns
        -------
        int: The maximum depth of the nested list

        Examples
        --------
            >> max_depth([1, [2, [3, 4]],5])
            >> 3
        """
        # This function determines the maximum depth of the nested list.
        if isinstance(lst, list):
            return max([depth] + [max_depth(item, depth + 1) for item in lst])
        else:
            return depth - 1

    pipettes, plates, tuberacks = labware
    if len(dfs_labware) == 4:
        df_plate, df_pipetteLeft, df_pipetteRight, df_stocks = dfs_labware
    else:
        df_stocks = dfs_labware[-1]
        print(f'... More than 4 labware dataframes found. Using the last one as the stock dataframe. '
              f'If this assumption were incorrect, revise the labware configuration file.')
    # TODO: look for a better way to handle the misassignment of the deck location
    # sys.exit()
    # sample.plate_deck_loc = int(df_plate.loc[1, 'deck_location'])
    deck2plate = get_deck2plate(plates)

    # For a sample with a simple order of mixing
    if isinstance(sample.well, str):
        order_mixing = sh.flatten_list(sample.order_mixing)  # ORDERED chemicals from the recipe
        stock_locator = get_stock_locator(tuberacks, df_stocks)  # Location of the stock tube

        # Location of the well to dispense the chemicals
        # sample.plate_deck_loc: Deck location (e.g. '10', '11')
        # opentrons plate object: plates[0].parent stores the deck location
        dst = deck2plate[sample.plate_deck_loc][sample.well]
        print(f'... Preparing sample {sample.id}, {sample.name} in well {dst.well_name} in Deck {sample.plate_deck_loc}')
        for j, chem in enumerate(order_mixing):
            if chem in sample.skip_species:
                print(f'... Skipping {chem} for {sample.name}')
            else:
                print(f'... STEP {j + 1}: Adding {chem} in well {sample.well}')
                src = stock_locator[chem]  # Location of the stock tube
                vol = sample.recipe[chem]  # Volume of the chemical to be dispensed
                kwargs = get_transfer_kwargs(sample.relative_viscosity[chem])  # Transfer kwargs for the chemicalkwargs
                kwargs['v_src'] = max(0, df_stocks.loc[chem, 'volume_mL'] * 0.9)  # Aspirate from 80% of the stock tube
                kwargs['v_dst'] = sample.volume  # Dispense at the total volume of the sample in the well
                if j == len(order_mixing) - 1: # If this is the last chemical in the list, mix the solution
                    kwargs['mix'] = True
                    kwargs['mix_kwargs'] = {'repetitions': 4, 'volume':  (sample.well_volume + vol) / 2, 'location': None, 'rate': 5}
                # Transfer the liquids
                pipette, last_location = control.transfer(protocol, pipettes, src, dst, vol, **kwargs)

                ## Kwarg options
                # pipette, last_location = control.transfer(protocol, pipettes, src, dst, vol,
                #                  flow_rate=1., aspirate_wait_time=0.5,
                #                  blowout_rate=None, max_blowout_height=200,
                #                  v_src=v_src, v_dst=v_dst, v2z_src=None, v2z_dst=None,
                #                  # default: Pass volume not height to determine the aspirating/dispensing heights
                #                  z_src=None, z_dst=None,  # alternative: pass height from the bshm
                #                  z_offset_src=0, z_offset_dst=0, # offset from the default height
                #                  # starting_tip=None, # Unused. Already configured
                #                  drop_tip=True, # Drop the tip after dispensing to avoid contamination
                #                  home=False,
                #                  default_z=5,
                #                  mix=mix,
                #                  mix_kwargs={'repetitions': 4, 'volume': 10, 'location': None, 'rate': 5},)

                # Update df_stocks
                df_stocks.loc[chem, 'volume_mL'] -= vol / 1000.  # Update the volume of the stock tube
                # Update the volume in the well
                sample.well_volume += vol

                control.drop_tip(pipette)
    # For a sample with a complex order of mixing
    elif isinstance(sample.well, list):
        order_mixing = sample.order_mixing # ORDERED chemicals from the recipe. e.g.- ['fus', ['fus_buffer','water', 'Proline']]
        stock_locator = get_stock_locator(tuberacks, df_stocks) # Location of the stock tube

        chem2well = sample.subwell # Dictionary to keep track of the chemicals mixed in the same well
        maxDepth = max_depth(order_mixing)

        print(f'... This sample requires {maxDepth} wells to prepare.\n')

        for j in range(max_depth(order_mixing)):
            prioritized_lists = sh.find_deepest_lists(order_mixing)
            for chem_list in prioritized_lists:
                print(f'... STEP{j+1}: Adding {chem_list} in well {sample.well[j]}')
                added_vol_in_well = 0
                for i, chem in enumerate(chem_list):
                    if chem in sample.skip_species:
                        print(f'... Skipping {chem} for {sample.name}')
                    else:
                        # Source well
                        src = stock_locator[chem]  # Labware['Tuberack1']['A1'] or [[Tuberack1['A1'], Tuberack2['A1'], ...]
                        # DEFAULT
                        if not isinstance(src, list):
                            # Destination well
                            if i == 0 and chem2well[chem] is None:
                                # If premixed chemicals are already in the first cell, move onto the next chemical
                                continue
                            else:
                                # Transfer everything to the first well of the list
                                chem2well[chem] = chem2well[chem_list[0]] # Set the destination well to the same well of the first list
                                dst = deck2plate[sample.plate_deck_loc][chem2well[chem]]
                            # Volume of the chemical to be dispensed
                            vol = sample.recipe[chem]
                            # Transfer kwargs for the chemical
                            kwargs = get_transfer_kwargs(sample.relative_viscosity[chem])
                            print(chem, df_stocks.loc[chem, 'volume_mL'])
                            kwargs['v_src'] = max(0, df_stocks.loc[chem, 'volume_mL'] * 0.9) # Aspirate from 80% of the stock tube
                            kwargs['v_dst'] = sample.volume # Dispense at the total volume of the sample in the well
                            # if i == len(chem_list) - 1:
                            #     kwargs['mix'] = True
                            #     kwargs['mix_kwargs'] = {'repetitions': 4, 'volume':  (sample.well_volume[chem2well[chem]] + vol) / 2, 'location': None, 'rate': 5}

                            # Transfer the liquids
                            pipette, last_location = control.transfer(protocol, pipettes, src, dst, vol, **kwargs)

                            # Update df_stocks
                            df_stocks.loc[chem, 'volume_mL'] -= vol / 1000.  # Update the volume of the stock tube

                            # Update the volume in the well
                            sample.well_volume[chem2well[chem]] += vol

                            # Drop the tip(s) to avoid contamination
                            control.drop_tips(pipettes)

                            # Update the volume in the well
                            added_vol_in_well += vol
                        # SPECIAL CASES (e.g. FUS- stock solution is in the well plates)
                        else:
                            # Volume of the chemical to be dispensed
                            vol = sample.recipe[chem]
                            stock_vols = df_stocks.loc[chem, 'volume_mL'] # list of stock volumes
                            # Find which stock tube to use
                            idx = [idx for idx, v in enumerate(stock_vols) if v >= vol/1000][0]

                            # Source well
                            src = stock_locator[chem][idx]

                            # Destination well
                            if i == 0 and chem2well[chem] is None:
                                # Set the destination well the same as the well where this chemical is placed
                                chem2well[chem] = stock_locator[chem].well
                                # If premixed chemicals are in the first cell, move onto the next chemical
                                continue
                            else:
                                # Volume of the chemical to be dispensed
                                vol = sample.recipe[chem]
                                # Transfer everything to the first well of the list
                                chem2well[chem] = chem2well[chem_list[0]] # Set the destination well to the same well of the first list
                                dst = deck2plate[sample.plate_deck_loc][chem2well[chem]]
                            # Transfer kwargs for the chemical
                            kwargs = get_transfer_kwargs(sample.relative_viscosity[chem])
                            kwargs['v_src'] = max(0, df_stocks.loc[chem, 'volume_mL'][idx] * 0.9) # Aspirate from 80% of the stock tube
                            kwargs['v_dst'] = sample.volume # Dispense at the total volume of the sample in the well
                            # if i == len(chem_list) - 1:
                            #     kwargs['mix'] = True
                            #     kwargs['mix_kwargs'] = {'repetitions': 4, 'volume':  (sample.well_volume[chem2well[chem]] + vol) / 2, 'location': None, 'rate': 5}


                            # Transfer the liquids
                            pipette, last_location = control.transfer(protocol, pipettes, src, dst, vol, **kwargs)

                            # Update df_stocks
                            df_stocks.loc[chem, 'volume_mL'][idx] -= vol / 1000.  # Update the volume of the stock tube

                            # Update the volume in the well
                            sample.well_volume[chem2well[chem]] += vol

                            # Drop the tip(s) to avoid contamination
                            control.drop_tips(pipettes)

                            # Update the volume in the well
                            added_vol_in_well += vol
                            print('\n')
                        print(f'...... Added {vol:.1f}uL of {chem} \n'
                              f'............ from {src}\n'
                              f'............   to {dst}\n')

                # Now update chem2well to indicate that the chemicals are mixed to this well
                # chem2well has a new entry: {'fus_buffer_water': 'A2'} for example to indicate that fus_buffer and water are mixed in well A2
                transferred_chem = '_'.join(chem_list)
                chem2well[transferred_chem] = dst.well_name
                # Update the location of the stock tube
                stock_locator[transferred_chem] = dst
                # Update the recipe of the sample. recipe['fus_buffer_water'] = vol_fus_buffer + vol_water
                sample.recipe[transferred_chem] = added_vol_in_well
                # Assume that the mixture is thin # TODO: Update this based on the added chemicals
                sample.relative_viscosity[transferred_chem] = 'thin'

                new_stock = {'type': 'temporary',
                             'concentration': None,
                             'volume_mL': added_vol_in_well,
                             'labware_name': dst._parent.name,
                             'deck_location': dst._parent.parent,
                             'rack_location': dst.well_name,
                             'tube_type': 'well_plate'
                             }
                df_stocks.loc[transferred_chem] = new_stock
                print(f'... {transferred_chem} was successfully prepared in well {dst.well_name} in Deck {dst._parent.parent}\n')
                print('-' * 50)
            # Prioritized chemicals are now mixed. Update the order_mixing by flattening the innermost list
            # The example is given below.
            # ['fus', ['fus_buffer', 'water'], ['Proline', 'NaCl']] -> ['fus', 'fus_buffer_'water', 'Proline_'NaCl']
            # This reads ['fus_buffer', 'water'] is prepared in the same well ('A2' for example)
            # and ['Proline', 'NaCl'] are prepared in the same well ('A3' for example).
            order_mixing = sh.flatten_innermost(order_mixing)
            # ['fus', ['fus_buffer', 'water'], ['Proline', 'NaCl']] -> ['fus', 'fus_buffer_'water', 'Proline_'NaCl']
        control.drop_tips(pipettes)

        sample.prepared = True
        sample.well_obj = dst
    else:
        print(f'... Incorrect well format for {sample.name}. It must be a string or a list of strings.')
    print(f'... Sample {sample.id}, {sample.name} was successfully prepared in well {sample.well} in Deck {sample.plate_deck_loc}')
    print('-' * 50)
    return df_stocks

def prepare_samples(protocol, labware, dfs_labware, list_of_samples):
    """
    Prepare the samples according to the recipe.
    ...
    ... With a simple order of mixing (e.g. ['fus', 'fus_buffer', 'water'])

    Parameters
    ----------
    protocol : protocol_api.ProtocolContext, The protocol context
    labware : tuple, Tuple of labware
        ... pipettes, plates, tuberacks = labware
    dfs_labware : list, List of labware dataframes
        ... df_plate, df_pipetteLeft, df_pipetteRight, df_stocks = dfs_labware
    list_of_samples : list, List of samples to be prepared
        ... samples = [sample_list1, sample_list2, ...]
                    = [[sample1, sample2, ..., sampleN], [sampleN+1, sampleN+2, ...], ..., [...]]

    Returns
    -------
    dfs_labware : list
    """
    def max_depth(lst, depth=1):
        """
        Returns the maximum depth of the nested list

        Parameters
        ----------
        lst: list, The list to be analyzed
        depth: int, The depth level of the list if it is not a nested list. ['a', 'b', 'c'] has a depth of 1.

        Returns
        -------
        int: The maximum depth of the nested list

        Examples
        --------
        >>> max_depth([1, [2, [3, 4]],5])
        3
        """
        # This function determines the maximum depth of the nested list.
        if isinstance(lst, list):
            return max([depth] + [max_depth(item, depth + 1) for item in lst])
        else:
            return depth - 1

    def _check_order_of_mixing(samples):
        """
        Check if all the samples have the same order of mixing
        Parameters
        ----------
        samples : list, list of samples

        Returns
        -------
        bool: True if all the samples have the same order of mixing
        """
        order_mixing = samples[0].order_mixing
        for sample in samples:
            if sample.order_mixing != order_mixing:
                return False
        return True

    def update_stock_locators(samples, tuberacks, df_stocks):
        """
        Add the stock locators to the given Sample objects
        Parameters
        ----------
        samples: list, List of Sample objects
            ... [sample1, sample2, ...,]
        stock_locator
            ... dict, Dictionary that maps a location of a chemical in the stock tube or well in the well plate
            ... e.g.- {'fus': Well('A1'), 'fus_buffer': Tuberack('A2'), ...}
            ... The values are opentrons Well objects not str!
        Returns
        -------

        """
        for sample in samples:
            sample.stock_locator = get_stock_locator(tuberacks, df_stocks)
        return samples

    def handle_simple_mixing(samples, protocol, pipettes, plates, df_stocks, chemicals):
        step = 1
        deck2plate = get_deck2plate(plates)
        for j, chem in enumerate(chemicals, start=1):
            print('-' * 50)
            print(f'... STEP {step}-{step+len(samples)}: Adding {chem} to Well {[s.well for s in samples]}')
            print('-' * 50)
            if j == len(chemicals): mix = True
            else: mix = False
            for i, sample in enumerate(samples):
                print(f'... Step {step}')
                if chem in sample.skip_species:
                    print(f'... Skipping {chem} for {sample.name}')
                else:
                    if i % 5 == 0:
                        print(f'... Distributing {chem} to sample {i + 1}/{len(samples)}')
                    print(f'...... Relative viscosity of {chem}: {sample.relative_viscosity[chem]}')

                    src = sample.stock_locator[chem]
                    dst = deck2plate[sample.plate_deck_loc][sample.well]
                    vol = sample.recipe[chem]

                    # Handle different types of `src`
                    if not isinstance(src, list):
                        df_stocks = transfer_single_source(protocol, pipettes, dst, vol, chem, sample, df_stocks, mix=mix)
                    else:
                        df_stocks = transfer_multiple_sources(protocol, pipettes, dst, vol, chem, sample, df_stocks, mix=mix)

                    print(f'...... Added {vol:.1f}uL of {chem} \n'
                          f'............ from {src}\n'
                          f'............   to {dst}\n')
                    sample.prepared = True
                    sample.well_obj = dst
                    step += 1
            print(f'... {chem} was successfully distributed to all samples\n')
            drop_tips(protocol)

        return df_stocks

    def transfer_single_source(protocol, pipettes, dst, vol, chem, sample, df_stocks, mix=False):
        """
        Transfer the chemical from a single source to the destination well
        ... From a single source (e.g. Tuberack1['A1']) to the destination well (e.g. Plate['A1'])

        Parameters
        ----------
        protocol: protocol_api.ProtocolContext
        pipettes: tuple, tuple of  pipette objects
        dst: opentrons.legacy_api.containers.placeable.Well, destination well
        vol: float, volume of the chemical to be transferred
        chem: str, name of the chemical
        sample: sh.Sample, sample object
        df_stocks: pd.DataFrame, stock dataframe
        mix: bool, True if the solution needs to be mixed

        Returns
        -------
        df_stocks: pd.DataFrame, updated stock dataframe
        """
        src = sample.stock_locator[chem]

        kwargs = get_transfer_kwargs(sample.relative_viscosity[chem])
        kwargs['v_src'] = max(0, df_stocks.loc[chem, 'volume_mL'] * 0.9)
        kwargs['v_dst'] = sample.volume

        if mix:
            kwargs['mix'] = True
            kwargs['mix_kwargs'] = {'repetitions': 4, 'volume': (sample.well_volume + vol) / 2, 'location': None,
                                    'rate': 5}

        pipette, last_location = control.transfer(protocol, pipettes, src, dst, vol, **kwargs)

        df_stocks.loc[chem, 'volume_mL'] -= vol / 1000.
        sample.well_volume += vol

        return df_stocks

    def transfer_multiple_sources(protocol, pipettes, dst, vol, chem, sample, df_stocks, mix=False):
        """
        Transfer the chemical from multiple sources to the destination well
        ... From a single source (e.g. Tuberack1['A1']) to the destination well (e.g. Plate['A1'])

        Parameters
        ----------
        protocol: protocol_api.ProtocolContext
        pipettes: tuple, tuple of  pipette objects
        dst: opentrons.legacy_api.containers.placeable.Well, destination well
        vol: float, volume of the chemical to be transferred
        chem: str, name of the chemical
        sample: sh.Sample, sample object
        df_stocks: pd.DataFrame, stock dataframe
        mix: bool, True if the solution needs to be mixed

        Returns
        -------
        df_stocks: pd.DataFrame, updated stock dataframe
        """
        stock_vols = df_stocks.loc[chem, 'volume_mL']
        idx = [i for i, v in enumerate(stock_vols) if v * 1000 >= max(vol, 1)][0]
        src = sample.stock_locator[chem][idx]

        kwargs = get_transfer_kwargs(sample.relative_viscosity[chem])
        kwargs['v_src'] = max(0, df_stocks.loc[chem, 'volume_mL'][idx] * 0.9)
        kwargs['v_dst'] = sample.volume

        if mix:
            kwargs['mix'] = True
            kwargs['mix_kwargs'] = {'repetitions': 4, 'volume': (sample.well_volume + vol) / 2, 'location': None,
                                    'rate': 5}

        pipette, last_location = control.transfer(protocol, pipettes, src, dst, vol, **kwargs)

        df_stocks.loc[chem, 'volume_mL'][idx] -= vol / 1000.
        sample.well_volume += vol

        return df_stocks

    def handle_complex_mixing(samples, protocol, pipettes, plates, df_stocks, order_mixing):
        step = 1
        deck2plate = get_deck2plate(plates)
        for j in range(max_depth(order_mixing)):
            # Identify the chemicals to be mixed
            prioritized_lists = sh.find_deepest_lists(order_mixing)
            for chem_list in prioritized_lists:
                # Loop over the chemicals to be mixed
                for i, chem in enumerate(chem_list):
                    print('-' * 50)
                    print(f'... STEP{step}: Adding {chem} to {len(samples)} wells')
                    print('-' * 50)

                    # Loop over samples
                    for k, sample in enumerate(samples):
                        print(f'...... Progress: {k + 1}/{len(samples)} for Sample {sample.name}')
                        if chem in sample.skip_species:
                            print(f'... Skipping {chem} for {sample.name}')
                        else:
                            # LOAD THE SAMPLE WELLS HERE
                            chem2well = sample.subwell  # Dictionary to keep track of the chemicals mixed in the same well
                            # Source well
                            src = sample.stock_locator[chem] # Labware['Tuberack1']['A1'] or [[Tuberack1['A1'], Tuberack2['A1'], ...]

                            # DEFAULT: SOURCE IS A TUBE RACK
                            if not isinstance(src, list):
                                # Destination well
                                if i == 0 and chem2well[chem] is None:
                                    # If premixed chemicals are already in the first cell, move onto the next chemical
                                    continue
                                else:
                                    # Set the destination well to the same well of the first list
                                    chem2well[chem] = chem2well[chem_list[0]]
                                    dst = deck2plate[sample.plate_deck_loc][chem2well[chem]]
                                # Volume of the chemical to be dispensed
                                vol = sample.recipe[chem]
                                # Transfer kwargs for the chemical
                                kwargs = get_transfer_kwargs(sample.relative_viscosity[chem])
                                kwargs['v_src'] = max(0, df_stocks.loc[
                                    chem, 'volume_mL'] * 0.9)  # Aspirate from 95% of the stock tube height
                                kwargs['v_dst'] = sample.volume  # Dispense at the total volume of the sample in the well

                                if i == len(chem_list) - 1:
                                    kwargs['mix'] = True
                                    kwargs['mix_kwargs'] = {'repetitions': 4,
                                                            'volume': (sample.well_volume[chem2well[chem]] + vol) / 2,
                                                            'location': None, 'rate': 5}

                                # Transfer the liquids
                                pipette, last_location = control.transfer(protocol, pipettes, src, dst, vol, **kwargs)
                                # Update df_stocks
                                df_stocks.loc[chem, 'volume_mL'] -= vol / 1000.  # Update the volume of the stock tube

                                # Update the volume in the well
                                sample.well_volume[chem2well[chem]] += vol
                                print('\n')
                            # SPECIAL CASES: when stock solutions are not set in the tube racks
                            #                (e.g. FUS- stock solution is in the well plates)
                            else:
                                # Volume of the chemical to be dispensed
                                vol = sample.recipe[chem]
                                stock_vols = df_stocks.loc[chem, 'volume_mL']  # list of stock volumes
                                # Find which stock tube to use
                                try:
                                    idx = [i for i, v in enumerate(stock_vols) if v * 1000 >= max(vol, 1)][0]
                                except:
                                    print(
                                        f'... Not enough stock solution for {chem}. Please refill the stock tube. Aborting the protocol.')
                                    sys.exit()

                                # Source well
                                src = sample.stock_locator[chem][idx]

                                # Destination well
                                if i == 0 and chem2well[chem] is None:
                                    # Set the destination well the same as the well where this chemical is placed
                                    chem2well[chem] = sample.stock_locator[chem].well
                                    # If premixed chemicals are in the first cell, move onto the next chemical
                                    continue
                                else:
                                    # Volume of the chemical to be dispensed
                                    vol = sample.recipe[chem]
                                    # Transfer everything to the first well of the list
                                    chem2well[chem] = chem2well[
                                        chem_list[0]]  # Set the destination well to the same well of the first list
                                    dst = deck2plate[sample.plate_deck_loc][chem2well[chem]]
                                # Transfer kwargs for the chemical
                                kwargs = get_transfer_kwargs(sample.relative_viscosity[chem])
                                kwargs['v_src'] = max(0, df_stocks.loc[chem, 'volume_mL'][
                                    idx] * 0.9)  # Aspirate from 80% of the stock tube
                                kwargs['v_dst'] = sample.volume  # Dispense at the total volume of the sample in the well

                                if i == len(chem_list) - 1:
                                    kwargs['mix'] = True
                                    kwargs['mix_kwargs'] = {'repetitions': 4,
                                                            'volume': (sample.well_volume[chem2well[chem]] + vol) / 2,
                                                            'location': None, 'rate': 5}

                                # Transfer the liquids
                                pipette, last_location = control.transfer(protocol, pipettes, src, dst, vol, **kwargs)

                                # Update df_stocks
                                df_stocks.loc[chem, 'volume_mL'][idx] -= vol / 1000.  # Update the volume of the stock tube

                                # Update the volume in the well
                                sample.well_volume[chem2well[chem]] += vol

                                if src.well_name == chem2well.keys():
                                    sample.well_volume[src.well_name] -= vol
                            print(f'...... Added {vol:.1f}uL of {chem} \n'
                                  f'............ from {src}\n'
                                  f'............   to {dst}\n')
                            print(f'...... Volume in well {chem2well[chem]}: {sample.well_volume[chem2well[chem]]:.1f}uL\n')

                            if i == len(chem_list) - 1:
                                # Now update chem2well to indicate that the chemicals are mixed to this well
                                # chem2well has a new entry: {'fus_buffer_water': 'A2'} for example to indicate that fus_buffer and water are mixed in well A2
                                transferred_chem = '_'.join(chem_list)
                                chem2well[transferred_chem] = dst.well_name
                                sample.subwell = chem2well  # Update the subwell dictionary
                                # Update the location of the stock tube
                                sample.stock_locator[transferred_chem] = dst

                                # Update the recipe of the sample. recipe['fus_buffer_water'] = vol_fus_buffer + vol_water
                                sample.recipe[transferred_chem] = sample.well_volume[dst.well_name]

                                # Assume that the mixture is thin # TODO: Update this based on the added chemicals
                                sample.relative_viscosity[transferred_chem] = 'thin'

                                # Update the stock dataframe to keep track of the intermediate solutions
                                new_stock = {'type': 'temporary',
                                             'concentration': None,
                                             'volume_mL': sample.well_volume[dst.well_name],  # TODO THIS IS WRONG
                                             'labware_name': dst._parent.name,
                                             'deck_location': dst._parent.parent,
                                             'rack_location': dst.well_name,
                                             'tube_type': 'well_plate'
                                             }
                                df_stocks.loc[transferred_chem] = new_stock
                    # Keep track of the number of steps
                    step += 1
                    control.drop_tips(pipettes)

            # Prioritized chemicals are now mixed. Update the order_mixing by flattening the innermost list
            # The example is given below.
            # ['fus', ['fus_buffer', 'water'], ['Proline', 'NaCl']] -> ['fus', 'fus_buffer_'water', 'Proline_'NaCl']
            # This reads ['fus_buffer', 'water'] is prepared in the same well ('A2' for example)
            # and ['Proline', 'NaCl'] are prepared in the same well ('A3' for example).
            order_mixing = sh.flatten_innermost(order_mixing)
            # ['fus', ['fus_buffer', 'water'], ['Proline', 'NaCl']] -> ['fus', 'fus_buffer_'water', 'Proline_'NaCl']
        control.drop_tips(pipettes)

        for sample in samples:
            sample.prepared = True
            sample.well_obj = deck2plate[sample.plate_deck_loc][sample.well]
        return df_stocks

    # Main function
    pipettes, plates, tuberacks = labware
    if len(dfs_labware) == 4:
        df_plate, df_pipetteLeft, df_pipetteRight, df_stocks = dfs_labware
    else:
        df_stocks = dfs_labware[-1]
        print(f'... More than 4 labware dataframes found. Using the last one as the stock dataframe. '
              f'If this assumption was incorrect, revise the labware configuration file.')

    for j, samples in enumerate(list_of_samples):
        step = 1
        print('-' * 50)
        print(f'Making samples {1+j}/{len(list_of_samples)}:', samples[0].name, '-', samples[-1].name)
        # Special case: If the samples do not have the same order of mixing
        if not _check_order_of_mixing(samples):
            print(f'... Samples {samples[0].name}-{samples[-1].name} do not have the same order of mixing.'
                  f' Preparing each sample one by one using "prepare_sample()".'
                  )
            for sample in samples:
                df_stocks = prepare_sample(protocol, labware, dfs_labware, sample)
        # Default case: If all the samples have the same order of mixing
        else:
            # For all samples with a simple order of mixing (no nested lists)
            if isinstance(samples[0].well, str):
                update_stock_locators(samples, tuberacks, df_stocks)  # Add the stock locators to the sample
                # Chemicals to be dispensed
                chemicals = samples[0].order_mixing
                df_stocks = handle_simple_mixing(samples, protocol, pipettes, plates, df_stocks, chemicals)

                # Update the stock dataframe
                dfs_labware[-1] = df_stocks

            # For all samples with a complex order of mixing (nested lists)
            elif isinstance(samples[0].well, list):
                # Initialization steps
                order_mixing = samples[0].order_mixing  # ORDERED chemicals from the recipe. e.g.- ['fus', ['fus_buffer','water', 'Proline']]
                chem2well = samples[0].subwell  # Dictionary to keep track of the chemicals mixed in the same well
                maxDepth = max_depth(order_mixing)
                print(f'... This sample requires {maxDepth} wells to prepare.\n')

                update_stock_locators(samples, tuberacks, df_stocks)  # Add the stock locators
                df_stocks = handle_complex_mixing(samples, protocol, pipettes, plates, df_stocks, order_mixing)

                # Update the stock dataframe
                dfs_labware[-1] = df_stocks
    protocol.home()
    return dfs_labware

def dilute_samples(protocol, labware, dfs_labware):
    """
    Dilute the samples according to /samples/dilution_configuration.csv.

    Parameters
    ----------
    protocol: protocol_api.ProtocolContext
    labware: tuple, Tuple of opentrons labware objects
    ... pipettes, plates, tuberacks = labware
    dfs_labware: list, List of labware dataframes
    ... df_plate, df_pipetteLeft, df_pipetteRight, df_stocks, df_protocol = dfs_labware

    Returns
    -------
    None
    """
    # Main function
    pipettes, plates, tuberacks = labware
    if len(dfs_labware) == 5:
        df_plate, df_pipetteLeft, df_pipetteRight, df_stocks, df_protocol = dfs_labware
    else:
        raise ValueError(f'... More than 5 labware dataframes found. Please revise the dilution configuration file.')
    # deck2plate: Deck number (int) to Plate object
    deck2plate = get_deck2plate(plates)

    # Name of the diluent (str)
    diluent = df_stocks.index.values[0]
    # stock_locator: dict, maps the diluent (str) to the location (Well obj)
    stock_locator = get_stock_locator(tuberacks, df_stocks)
    # Location of diluent
    diluent_loc = stock_locator[diluent] # diluent_loc: Well Object
    # Volume of the diluent stock tube
    vol_diluent_stock = df_stocks.loc[diluent, 'volume_mL'] # mL to uL


    # Figure out the source wells and destination wells
    src_dn = int(df_protocol['source_deck_location'].values[0])
    src_start_well, src_end_well = df_protocol['source_wells'].values[0].split('-')
    src_inc = int(df_protocol['source_increment'].values[0])
    src_volume = float(df_protocol['source_volume_mL'].values[0]) # Volume in the source well
    dst_dn = [int(dn) for dn in df_protocol['destination_deck_location'].values[0].split(',')]
    dst_start_well = df_protocol['destination_starting_well'].values[0]
    dst_volume = float(df_protocol['destination_volume_mL'].values[0])
    repeat = int(df_protocol['repeat'].values[0])
    dilution_factor = float(df_protocol['dilution_factor'].values[0])
    vol = dst_volume / dilution_factor # Volume of the source to be transferred

    src_wells = sh.get_wells(start=src_start_well, end=src_end_well, inc=src_inc)
    dst_wells = sh.get_wells(n=len(src_wells)*repeat, start=dst_start_well, inc=1) # list of strings

    # Find the source wells (Well obj)
    src_Wells = [deck2plate[src_dn][well] for well in src_wells for _ in range(repeat)] # Well objects
    # Assign the destination wells
    dst_Wells = [] # Well objects
    for dn in dst_dn:
        for j, well in enumerate(dst_wells):
            dst_Wells.append(deck2plate[dn][well])
            if well == 'H12':
                break
        dst_wells = dst_wells[j+1:]
    # Initialize
    df = pd.DataFrame(columns=['plate_deck_loc',
                               'well',
                               'sample_id',
                               'name',
                               'volume',
                               'composition',
                               'recipe',
                               'order_mixing'
                               ])
    df.to_csv('./preparation/dilution_info.csv', index=False)

    # Transfer the diluent to the destination wells
    for src, dst in zip(src_Wells, dst_Wells):
        print(f'... Diluting {src.well_name}({src._parent.parent}) to {dst.well_name}({dst._parent.parent})\n'
              f'...... Dilution factor: {dilution_factor}\n'
              f'...... Solvent: {diluent}\n'
              f'...... Final volume: {dst_volume} uL\n')

        transfer_src_kwargs = copy.deepcopy(transfer_default_kwargs)
        transfer_src_kwargs['v_src'] = (src_volume - vol*repeat) # 0.9 is a user-adjusted parameter
        transfer_dil_kwargs = copy.deepcopy(transfer_default_kwargs)
        transfer_dil_kwargs['v_src'] = vol_diluent_stock - vol  # 0.9 is a user-adjusted parameter
        transfer_dil_kwargs['mix'] = True  # 0.9 is a user-adjusted parameter
        transfer_dil_kwargs['mix_kwargs'] = {'repetitions': 4, 'volume': dst_volume / 2, 'location': None, 'rate': 5}

        print(vol_diluent_stock, vol)
        control.dilute(protocol, pipettes, src, dst, vol, diluent_loc, dilution_factor,
                       transfer_src_kwargs=transfer_src_kwargs,
                       transfer_dil_kwargs=transfer_dil_kwargs,)

        # Update volume in the stock tube
        vol_diluent_stock -= vol/1000

        # Update the sample_info.csv file
        sh.update_dilution_info(info={'plate_deck_loc': dst._parent.parent,
                                'well': dst.well_name,
                                'sample_id': None,
                                'name': f'deck{src._parent.parent}_{src.well_name}_dilution_x{dilution_factor:.0f}',
                                'volume': dst_volume,
                                'composition': None,
                                'recipe': None,
                                'order_mixing': ['sample', 'diluent']},
                            )

# DESIGN EXPERIMENTS HERE
def execute_default_protocol(protocol: protocol_api.ProtocolContext):
    """
    Prepare the samples per list of samples
    ... Assumes that each list of samples uses the same chemicals
    ... Appropriate for variation of a single chemical in the samples

    Parameters
    ----------
    protocol: protocol_api.ProtocolContext

    Returns
    -------

    """

    ## DIFFERENT METHOD
    labware, dfs_labware = control.load_labware(protocol, filepath=master_dir / 'preparation/labware_configuration.csv')
    samples_list = sh.read_samples(sample_dir, flatten=False)  # List of a list of Sample objects
    samples_flattened = sh.flatten_list(samples_list)  # List of Sample objects
    nSamples = len(samples_flattened)
    print(f'... Found {nSamples} samples to prepare')

    dfs_labware = prepare_samples(protocol, labware, dfs_labware, samples_list)
    print('\nDone!')

    print('-' * 50)
    print('Final Stock Volumes:')
    print(dfs_labware[-1].loc[:, 'volume_mL'])
    print('-' * 50)
    return dfs_labware, samples_flattened
def execute_protocol_one_by_one(protocol: protocol_api.ProtocolContext, start=0):
    """
    Prepare the samples one by one
    ... This is the most robust method but is slow

    Parameters
    ----------
    protocol: protocol_api.ProtocolContext

    Returns
    -------
    dfs_labware: list, list of labware dataframes
    all_samples: list of Sample objects
    """
    labwares, dfs_labware = control.load_labware(protocol, filepath=master_dir / 'preparation/labware_configuration.csv')
    all_samples = sh.read_samples(sample_dir, flatten=True)  # List of Sample objects
    nSamples = len(all_samples)
    print(f'... Found {nSamples} samples to prepare')

    # 1. Prepare the sample one by one (MOST ROBUST but SLOW)
    for i, sample in enumerate(all_samples):
        if i >= start:
            if i % 10 == 0:
                print(f'Preparing sample {i + 1}/{nSamples}')
            df_stocks = prepare_sample(protocol, labwares, dfs_labware, sample)
            dfs_labware[-1] = df_stocks  # Update the stock dataframe

    print('\nDone!')

    print('-'*50)
    print('Final Stock Volumes:')
    print(dfs_labware[-1].loc[:, 'volume_mL'])
    print('-'*50)
    return dfs_labware, all_samples

def execute_dilution_protocol(protocol):
    try:
        labware, dfs_labware = control.load_labware(protocol, filepath=master_dir / 'preparation/dilution_configuration.csv')
    except:
        control.unload_labware(protocol)
        labware, dfs_labware = control.load_labware(protocol, filepath=master_dir / 'preparation/labware_configuration.csv')


    dilute_samples(protocol, labware, dfs_labware)

    # samples_list = sh.read_samples(sample_dir, flatten=False)  # List of a list of Sample objects
    # samples_flattened = sh.flatten_list(samples_list)  # List of Sample objects
    # nSamples = len(samples_flattened)
    # print(f'... Found {nSamples} samples to prepare')

    # dfs_labware = prepare_samples(protocol, labware, dfs_labware, samples_list)
    print('\nDone!')

    print('-'*50)
    print('Final Stock Volumes:')
    print('-'*50)
    return dfs_labware

#
def drop_tips(protocol):
    """
    Drop the tips given the protocol object
    ... control.drop_tips(pipettes) takeks the pipettes as an argument not the protocol object

    Parameters
    ----------
    protocol: protocol_api.ProtocolContext object

    Returns
    -------

    """
    # Specify the Pipette objects from protocol (must be preloaded)
    pipettes = [instrument for i, instrument in enumerate(protocol._instruments.values()) if i < 2]
    control.drop_tips(pipettes)

def save_run(out_dir=out_dir):
    """
    Save the files required for the run
    Parameters
    ----------
    out_dir: pathlib.Path

    Returns
    -------

    """
    if isinstance(out_dir, Path):
        out_dir = Path(out_dir)

    timestamp = (datetime.datetime.now() - datetime.timedelta(hours=4)).strftime('%Y-%m-%d-%H-%M-%S')
    print(f'Saved the files at {out_dir / timestamp}.')
    for dirname in ['preparation', 'samples']:
        src = master_dir / dirname
        dst = out_dir / timestamp / dirname
        shutil.copytree(src, dst)

    return out_dir / timestamp

def upload_exerimental_protocol(ip_address='10.49.62.184', ot2_ssh_key='~/.ssh/ot2_ssh_key',
                                src=os.getcwd(), dst='/var/lib/jupyter/notebooks/'):
    """
    This function uploads and updates relevant codes to run a job using sh.
    ... 1. It removes the existing folder /sh/.../sh_sample_prep
    ... 2. Then, it replaces it with the local folder /LOCAL/.../sh_sample_prep
    ... 2. Then, it replaces it with the local folder /LOCAL/.../sh_sample_prep

    Compatible with Windows and OS/Linux

    Parameters:
    ip_address: str, an IP address of sh
    ot2_ssh_key: str, location to the SSH private key to connect to sh
    src: str, a path to the directory where sh-related codes are stored
    dest: str, a path to the directory where the `src` folder will be copied
    """
    directory_name = 'sh_sample_prep'
    if os.name == 'nt': # Windows
        # Remove files from previous experiments on the remote server
        remove_command = f"ssh -i {ot2_ssh_key} root@{ip_address} rm -r {dst}{directory_name}"
        subprocess.run(remove_command, shell=True)

        # Upload files to remote server
        upload_command = f"scp -O -i {ot2_ssh_key} -r {src} root@{ip_address}:{dst}"
        result = subprocess.run(upload_command, shell=True)
    elif os.name == 'posix':  # Mac OS
        # Remove files from previous experiments on the remote server
        remove_command = f"ssh -i {ot2_ssh_key} {dst.split(':')[0]} 'rm -r {dst}{directory_name}'"
        subprocess.run(remove_command, shell=True)

        # Upload files to remote server
        upload_command = f"scp -O -r -i {ot2_ssh_key} {src} {dst}"
        result = subprocess.run(upload_command, shell=True)
    else:
        print('... Unkonwn OS!')
        sys.exit(1)
    print("Return code:", result.returncode)
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)
    print(upload_command)
    print('... Uploaded the files (run files, python scripts, etc.) to sh successfully.')


# EXPERIMENTAL DESIGNS
## OT-2 will prepare the samples according the run function below
def run(protocol: protocol_api.ProtocolContext, start=0):
    """
    Run the protocol to prepare the samples

    Parameters
    ----------
    protocol
    start: int, the index of the sample to start from (default is 0)

    Returns
    -------

    """
    print('-'*50)
    print(f'Protocol API Version: ', protocol._api_version)
    print(f'Run has started at {(datetime.datetime.now() - datetime.timedelta(hours=4)).strftime("%Y-%m-%d-%H-%M-%S")}.')

    # PROTOCOL
    # df_stocks, samples = execute_default_protocol(protocol) # Execute the protocol per list of samples (Recommended)
    df_stocks, samples = execute_protocol_one_by_one(protocol, start=start) # Execute the protocol one by one (Changes a tip every time, robust but slow)

    # UNDER CONSTRUCTION--------
    # wells_obj = [s.well_obj for s in samples if 'prep' not in s.name]]
    # execute_dilution_protocol(protocol)
    # --------------------------

    protocol.home()

    save_run(out_dir)
    print(f'Run has ended at {(datetime.datetime.now() - datetime.timedelta(hours=4)).strftime("%Y-%m-%d-%H-%M-%S")}.')
    print('... All samples were successfully prepared. Good job!')
    print('-'*50)
    print('-'*50)


if __name__ == '__main__':
    print("To submit a job on sh, use the command 'opentrons_execute startProtocol.py -n'")

    sh.remove_sample_list() # Remove the sample list

    # STEP 1: Create /preparation/stocks.csv
    ## Method 1: Create a dictionary of the stock solutions
    # stock_dict = {
    #     'fus': 53.9,  # g/L fus=1.07uM + 6M urea in 10 uL
    #     'fus_buffer': 1000 / 3.,  # 50 mM Hepes, 150mM NaCl,
    #     'Proline': 1000,
    #     'NaHCO3': 1000,
    #     'NH4OH': 1000,
    #     'Lysine': 1000,
    #     'water': 0,
    # }
    # sh.write_stocks(stock_dict)
    ## Method 2 Directly edit the csv file at /preparation/stocks.csv
    ### If you edit the csv file directly, you must load the file.

    ans = input('... Is there a ./preparation/stocks.csv file? (y/n)')
    if ans == 'y':
        stock_dict = sh.read_stocks(filepath='./preparation/stocks.csv')
    else:
        print('... Please create a ./preparation/stocks.csv file.')
        sys.exit(1)

    # STEP 2: Create Sample objects
    ## A Sample object contains `composition` and `order_mixing`
    ## If you pass, `stock_dict` and `volume`, it will create a recipe automatically. (Recommended)
    all_samples = []
    for solute in ['Proline', 'NaHCO3']: #  ['NH4OH', 'Proline', 'NaHCO3', 'Lysine']
        # `samples` is a list of Sample objects.
        samples = sh.create_samples(composition={'fus_buffer': 50., solute: 0, 'water': 0},
                                      vary=solute, concentrations=[0, 25, 50, 100, 200, 300, 400, 500, 600, 700],
                                      solvent='water',
                                      order_mixing=[solute, 'fus_buffer', 'water'],
                                      stock_dict=stock_dict, volume=118,
                                      # Provide with stock_dict and sample volume to create a recipe
                                      name='fus_' + solute,
                                      save=True,
                                      )
        all_samples.append(samples)
    #
    # STEP 3: Create a labware configuration file
    ans = input('... Is there a ./preparation/lab_configuration.csv file? (y/n)')
    if ans == 'n':
        sh.create_default_labware_configuration(all_samples)

    # STEP 4: USER MUST EDIT THE FILE at /preparation/lab_configuration.csv
    input('... Please edit the labware configuration file at /preparation/lab_configuration.csv '
          'and press Enter to continue.')

    # STEP 5: Read the labware configuration file
    dfs = sh.read_labware_configuration()

    # STEP 6: Assign a unique well position to each sample
    sh.assignWellPositions(all_samples)

    # STEP 7: (Optional) Create a sample sheet (Sample name, well number, composition, etc.)
    sh.create_sample_info(all_samples)

    # STEP 8: (Optional) Simulate the protocol
    ans = input('... Do you want to simulate the protocol? (y/n)')
    if ans == 'y':
        # Simulate the protocol
        command = ['opentrons_simulate', 'startProtocol.py', '>', './archive/out.txt']
        subprocess.run(' '.join(command), shell=True)
    print('... Protocol was successfully simulated.')
    # STEP 9: Upload the experimental protocol to the OT-2
    upload_exerimental_protocol() # Upload the experimental protocol to the OT-2

    # STEP 10: Run the protocol
    ans = input('... Do you want to run the protocol? (y/n)')
    if ans == 'y':
        command = [
            'ssh', '-i', '~/.ssh/ot2_ssh_key', 'root@10.49.35.97',
            'cd /var/lib/jupyter/notebooks/sh_sample_prep && opentrons_execute ./startProtocol.py &'
        ]

        # Run the command
        result = subprocess.run(' '.join(command), shell=True)
        output = result.stdout
        print(output)
    else:
        print('... Protocol was not executed.')

