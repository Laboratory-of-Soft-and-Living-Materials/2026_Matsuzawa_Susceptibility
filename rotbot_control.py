"""
Module to control sh for sample preparation

Author: Kaarthik, Takumi
MIT License
"""

import csv
import datetime
import inspect
import os
import pickle
import shutil
import platform
import pathlib
from pathlib import Path
import math
import numpy as np
import sample_handling as sh

# Necessary for Path library to work on Windows and Linux
plt = platform.system()
if plt == 'Linux': pathlib.WindowsPath = pathlib.PosixPath

# FILE ARCHITECTURE
master_dir = Path(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
samples_dir = master_dir / 'samples'

# AVAIlABLE LAB WARE
plate_types = ['biorad_96_wellplate_200ul_pcr', 'thermoscientificnunc_96_wellplate_2000ul']
pipette_types = ['p20_single_gen2', 'p300_single_gen2', 'p1000_single_gen2']
tiprack_types = ['opentrons_96_tiprack_20ul', 'opentrons_96_tiprack_300ul', 'opentrons_96_tiprack_1000ul']

# Tip rack offsets
tiprack_offsets = {'opentrons_96_tiprack_20ul': [0.00, 0.00, 0.00],
                   'opentrons_96_tiprack_300ul': [0.00, 0.00, 0.00],
                   'opentrons_96_tiprack_1000ul': [0.10, 1.00, 0.00]}

# Transfer kwargs for the transfer function
transfer_default_kwargs = {'flow_rate': 1.0, 'aspirate_wait_time': 0.5, 'hang_time': 0.5,
                        'blowout_rate': None, 'max_blowout_height': 200,
                        'v_src': 0, 'v_dst': 0, 'v2z_src': None, 'v2z_dst': None,
                        'z_src': None, 'z_dst': None,
                        'z_offset_src': 0, 'z_offset_dst': 0, 'drop_tip': False,
                        'default_z': 5,
                        'mix': False, 'mix_kwargs': {'repetitions': 4, 'volume': 20, 'location': None, 'rate': 5}
                        }
## Transfer kwargs for the diluent (mix is True)
transfer_diluent_kwargs = {'flow_rate': 1.0, 'aspirate_wait_time': 0.5, 'hang_time': 0.5,
                        'blowout_rate': None, 'max_blowout_height': 200,
                        'v_src': 0, 'v_dst': 0, 'v2z_src': None, 'v2z_dst': None, 'z_src': None, 'z_dst': None,
                        'z_offset_src': 0, 'z_offset_dst': 0, 'drop_tip': False,
                        'default_z': 5,
                        'mix': True, 'mix_kwargs': {'repetitions': 4, 'volume': 20, 'location': None, 'rate': 5}
                        }

def load_labware(protocol, filepath=master_dir / 'preparation/labware_configuration.csv'):
    '''
    Load all the labware required for sh for the protocol which is specificied in the run.csv file
    '''

    def configure_pipette(pipette, tipracks, starting_tip):
        """
        Function to configure the pipette with the tip racks and starting tip positions

        Parameters
        ----------
        pipette: opentrons.protocol_api.instrument_context.InstrumentContext object
            ... The pipette to be configured
        tipracks: list of opentrons.protocol_api.labware.Labware objects
            ... The tip racks for the pipette
        starting_tip: str
            ... The starting tip positions in the tip racks
            ... e.g. 'A1', 'B2'

        Returns
        -------
        pipette: opentrons.protocol_api.instrument_context.InstrumentContext object
            ... The pipette with the tip racks and starting tip positions configured
        """
        pipette_max_vol = str(pipette.max_volume)
        tip_vol = sh.extract_str(tipracks[0].name, 'tiprack_', 'ul')
        if pipette_max_vol == tip_vol:
            pipette.tip_racks = tipracks # Set the tip racks
            if len(tipracks) > 0: # Set the starting tip
                pipette.starting_tip = tipracks[0].well(starting_tip)
        return pipette

    def load_pipette(protocol, df_pipette, attributes=['name', 'deck_location', 'starting_tip']):
        """
        Load the pipette based on the configuration file

        Parameters
        ----------
        protocol: opentrons.protocol_api.protocol_context.ProtocolContext
        df_pipette: pandas.DataFrame
            ... The dataframe containing the pipette configuration
        attributes: list of strings, default=['name', 'deck_location']
            ... The attributes to be read from the configuration file. The first attribute is the name of the pipette.

        Returns
        -------
        pipette: opentrons.protocol_api.instrument_context.InstrumentContext object
        """
        tipracks = []
        starting_tips = []

        for index, row in df_pipette.iterrows():
            name = row[attributes[0]]
            loc = row[attributes[1]]
            starting_tip = row[attributes[2]]
            if loc in ['left', 'right']:
                pipette = protocol.load_instrument(name, loc)
                print('... Loaded pipette:', name, 'at location:', loc)
            else:
                tiprack = protocol.load_labware(name, loc)
                tiprack.set_offset(x=tiprack_offsets[name][0], y=tiprack_offsets[name][1], z=tiprack_offsets[name][2])
                tipracks.append(tiprack)
                starting_tips.append(starting_tip)
                print('... Loaded tiprack:', name, 'at location:', loc, 'with starting tip:', starting_tip)
        starting_tip = starting_tips[0] # Use the first tiprack for the starting tip. Other entries are ignored.
        configure_pipette(pipette, tipracks, starting_tip)
        return pipette, tipracks

    def load_plate(protocol, df_plate, attributes=['name', 'deck_location']):
        """
        Load the plate based on the configuration file

        Parameters
        ----------
        protocol: opentrons.protocol_api.protocol_context.ProtocolContext
        df_plate: pandas.DataFrame
            ... The dataframe containing the plate configuration

        Returns
        -------
        plate: opentrons.protocol_api.labware.Labware object
        """
        plates = []
        for index, row in df_plate.iterrows():
            name = row[attributes[0]]
            loc = row[attributes[1]]
            plate = protocol.load_labware(name, loc)
            print('... Loaded plate:', name, 'at Deck:', loc)
            plates.append(plate)
        return plates

    def load_stocks(protocol, df_stocks, attributes=['labware_name', 'deck_location']):
        """
        Load the stocks based on the configuration file

        Parameters
        ----------
        protocol: opentrons.protocol_api.protocol_context.ProtocolContext
        df_stocks: pandas.DataFrame
            ... The dataframe containing the stock configuration

        Returns
        -------
        stocks: dictionary
        """
        # Initialize the tuberacks
        tuberacks = [] # tuberack

        # Find the unique labwares
        labwares = list(set(df_stocks['labware_name'].values))
        for index, row in df_stocks.iterrows():
            name = row[attributes[0]]
            loc = row[attributes[1]]
            if name in labwares:
                tuberack = protocol.load_labware(name, loc)
                labwares.remove(name) # Load only once
                tuberacks.append(tuberack)
                print('... Loaded tuberack:', name, 'at Deck', loc)
        return tuberacks

    print('-'*50)
    print('Loading labware...')

    # Read the labware configuration file
    dfs_labware = sh.read_labware_configuration(filepath=filepath.__str__())
    df_plate, df_pipetteLeft, df_pipetteRight, df_stocks = dfs_labware[:4]

    # Load the pipettes and tipracks
    pLeft, tpLeft = load_pipette(protocol, df_pipetteLeft)
    pRight, tpRight = load_pipette(protocol, df_pipetteRight)
    pipettes, tipracks = [pLeft, pRight], [tpLeft, tpRight]

    # Load the plate
    plates = load_plate(protocol, df_plate)

    # Load the tuberacks (stocks)
    tuberacks = load_stocks(protocol, df_stocks)

    labware = [pipettes, plates, tuberacks]
    print('-'*50)
    return labware, dfs_labware

def unload_labware(protocol):
    """
    Unload all the labware from the protocol_api

    Parameters
    ----------
    protocol: opentrons.protocol_api.protocol_context.ProtocolContext

    Returns
    -------
    None
    """

    # Pipettes
    for name in protocol._instruments.keys():
        protocol._instruments[name] = None
    # Deck
    for deck_num in range(1, 12): # Do not delete 12 (trash)!
        protocol._deck.__delitem__(deck_num)

def choose_pipette(pipettes, vol):
    """
    Function to choose a pipette based on the volume to be transferred

    Parameters
    ----------
    pipettes: list of opentrons.protocol_api.instrument_context.InstrumentContext objects
    vol: float, volume of the liquid to be transferred in uL

    Returns
    -------
    pipette: opentrons.protocol_api.instrument_context.InstrumentContext object
        ... The pipette to be used for the transfer
    """
    # [1, 20], [20, 300]
    if pipettes[0].min_volume <= pipettes[1].min_volume:
        smallPipette = pipettes[0]
        largePipette = pipettes[1]
    else:
        smallPipette = pipettes[1]
        largePipette = pipettes[0]

    if vol <= smallPipette.min_volume:
        return smallPipette
    elif vol > largePipette.max_volume:
        return largePipette
    elif vol < smallPipette.max_volume and vol > smallPipette.min_volume:
        return smallPipette
    else:
        return largePipette

def pick_up_pipette(pipette):
    """
    Function to pick up a tip for the pipette

    Parameters
    ----------
    pipette: opentrons.protocol_api.instrument_context.InstrumentContext object

    Returns
    -------
    None
    """
    if pipette.has_tip == False:
        pipette.pick_up_tip()

def drop_tip(pipette):
    """
    Function to drop the tip for the pipette

    Parameters
    ----------
    pipette: opentrons.protocol_api.instrument_context.InstrumentContext object

    Returns
    -------
    None
    """
    if pipette.has_tip == True:
        pipette.drop_tip()
def drop_tips(pipettes):
    """
    Function to drop the tip for all the the pipettes

    Parameters
    ----------
    pipettes: list of opentrons.protocol_api.instrument_context.InstrumentContext objects

    Returns
    -------

    """
    for p in pipettes:
        try:
            if p.has_tip == True:
                p.drop_tip()
        except:
            pass


# BASIC FUNCTIONS
## `transfer`: Transfer liquid from one well to another
## `dilute`: Create a diluted solution of a sample
def transfer(protocol, pipettes, src, dst, vol,
             flow_rate=1., aspirate_wait_time=0.5, hang_time=0.5,
             blowout_rate=None, max_blowout_height=200,
             v_src=0., v_dst=0., v2z_src=None, v2z_dst=None,
             # default: Pass volume not height to determine the aspirating/dispensing heights
             z_src=None, z_dst=None,  # alternative: pass height from the bshm
             z_offset_src=0, z_offset_dst=0,
             starting_tip=None, drop_tip=False, home=False,
             default_z=5,
             mix=False, mix_kwargs={'repetitions ': 4, 'volume': 10, 'location': None, 'rate': 5},
             ):
    """
    Function to transfer liquid from one well to another

    Parameters
    ----------
    protocol: opentrons.protocol_api.protocol_context.ProtocolContext
    pipettes: list of opentrons.protocol_api.instrument_context.InstrumentContext objects
    src: opentrons.protocol_api.labware.Labware object, e.g.- stock_tube_racks[0]['A1']
    dst: opentrons.protocol_api.labware.Labware object  e.g.- plate[0]['A1']
    vol: float, volume of the liquid to be transferred in uL
    flow_rate: float, rate of aspiration and dispensing, default=1
    aspirate_wait_time: float, time to wait after aspirating in seconds, default=0.5
    blowout_rate: float, rate of blowout, default=None
    max_blowout_height: float, maximum height to blow out, default=200
    v_src: float, volume (mL) that corresponds to the height at which liquid is aspirated from the source well
            default=0 (bshm)
            e.g. If the stock solution contains 10mL, then provide 10.
    v_dst: float, volume (mL) that corresponds to the height at which liquid is dispensed in the destination well
            default=0 (bshm)
            e.g. If the destination well contains 10mL, then provide 10.
    v2z_src: function, volume to height converter for the source well, default=None
            By default, it will automatically choose based on the labware.
    v2z_dst: function, volume to height converter for the destination well, default=None
            By default, it will automatically choose based on the labware.
    starting_tip: string, starting tip position, default=None. This is only relevant for the first time it grabs a tip
    drop_tip: bool, whether to drop the tip after the transfer, default=False
    home: bool, whether to home the pipette after the transfer, default=False
    default_z: float, default height from the top of the well, default=5
    mix: bool, whether to mix the liquid in the destination well, default=False
    mix_kwargs: dictionary, mixing parameters, default={'repetitions': 4, 'volume': 10, 'location': None, 'rate': 5}
        ... These kwargs will be passed to the pipette.mix() function.
    Returns
    -------
    pipette: opentrons.protocol_api.instrument_context.InstrumentContext object
    last_location: opentrons.protocol_api.labware.Labware object
        ... The labware where the last transfer was made

    Examples
    --------
        src = stock_tube_racks[1]['A1'] # FUS
        dst = plate[0]['A1']
        transfer_from_A_to_B(protocol, pipettes, src, dst, 20,
                     v_src=9, v_dst=0,
                     starting_tip='B1')

    """
    # If the volume is zero, then don't do anything
    if vol <= 0:
        # Choose a pipette based on the volume to be transferred
        pipette = choose_pipette(pipettes, vol)
        return pipette, None

    def reset_blow_out_rate(pipette, default_blow_out_flow_rate):
        pipette.flow_rate.blow_out = default_blow_out_flow_rate

        # A function that maps a labware to an appropriate volume-to-height covnerter

    v2z = {'opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical': {'A1': sh.v2z_falcon15ml,
                                                                  'A2': sh.v2z_falcon15ml,
                                                                  'A3': sh.v2z_falcon50ml,
                                                                  'A4': sh.v2z_falcon50ml,
                                                                  'B1': sh.v2z_falcon15ml,
                                                                  'B2': sh.v2z_falcon15ml,
                                                                  'B3': sh.v2z_falcon50ml,
                                                                  'B4': sh.v2z_falcon50ml,
                                                                  'C1': sh.v2z_falcon15ml,
                                                                  'C2': sh.v2z_falcon15ml,
                                                                  'D1': sh.v2z_falcon15ml,
                                                                  'D2': sh.v2z_falcon15ml,
                                                                  },
           'thermoscientificnunc_96_wellplate_2000ul': {f"{letter}{number}": sh.v2z_deepwellplate for letter in 'ABCDEFGH'
                                                        for number in range(1, 13)},
           'biorad_96_wellplate_200ul_pcr': {f"{letter}{number}": sh.v2z_200uLwellplate for letter in 'ABCDEFGH' for
                                             number in range(1, 13)},
           'opentrons_96_tiprack_300ul': {f"{letter}{number}": sh.v2z_200uLeppendorf for letter in 'ABCDEFGH' for number in
                                          range(1, 13)}, # 300uL tipracks are used to place 200uL eppendorf tubes.
                                         }
    ## FIGURE OUT THE ASPIRATING/DISPENSING HEIGHT
    # Detemine the aspirating height
    if z_src is None:
        # Option 1: Pass volume of the solution in the source, then convert it to height
        if v2z_src is None: v2z_src = v2z[src.parent.name][src.well_name]
        z_src_in_mm = v2z_src(v_src)
    else:
        # Option 2: Pass height directly
        z_src_in_mm = z_src
    # Detemine the dispensing height
    if z_dst is None:
        # Option 1: Pass volume of the solution in the source, then convert it to height
        if v2z_dst is None: v2z_dst = v2z[dst.parent.name][dst.well_name]
        z_dst_in_mm = v2z_dst(v_dst)
    else:
        # Option 2: Pass height directly
        z_dst_in_mm = z_dst

    ## PIPETTE SETTINGS
    # Choose a pipette based on the volume to be transferred
    pipette = choose_pipette(pipettes, vol)

    ## Pipetting settings
    # Starting tip position
    if starting_tip is not None:
        # This is a Well class object. This is only relevant for the first time it grabs a tip
        pipette.starting_tip = pipette.starting_tip.parent[starting_tip]

    # Adjust the blow-out rate
    default_blow_out_flow_rate = pipette.flow_rate.blow_out
    # Default_flow_rate_blow_out
    if blowout_rate is not None:
        pipette.flow_rate.blow_out = default_blow_out_flow_rate

    # MOVEMENT
    ## Pick up a tip
    if pipette.has_tip == False:
        pipette.pick_up_tip()

    # MOVEMENT
    ## 1. Aspirate (fluid extraction) and wait to ensure the liquid is fully aspirated into the tip
    pipette.aspirate(vol, src.bshm(z=z_offset_src + z_src_in_mm), rate=flow_rate)
    protocol.delay(seconds=aspirate_wait_time)

    ## 2. Lift and wait
    ## This is to let the liquid flow down to the tip, preventing any liquid to drip out of the tip during the movement
    pipette.move_to(src.top(z=default_z))
    protocol.delay(seconds=hang_time)

    ## 3. Move to the destination well
    pipette.move_to(dst.top(z=default_z))

    pipette.move_to(dst.bshm(z=z_offset_dst + z_dst_in_mm))

    ## 4. Dispense
    pipette.dispense(vol, dst.bshm(z=z_offset_dst + z_dst_in_mm), rate=flow_rate)

    ## 5. Mix (Optional)
    if mix:
        if mix_kwargs['volume'] > pipette.max_volume:
            pipette.drop_tip()
            # Pick up a large pipette
            pipette = choose_pipette(pipettes, mix_kwargs['volume'])
            if pipette.has_tip:
                pipette.drop_tip()
            else:
                pipette.pick_up_tip()

        print(f'... Mixing {mix_kwargs["volume"]:.1f}', f'uL in the destination well.')

        # TM- z-value might have to be adjusted. 0.5 overfills the well at the tip insertion. 0.8 is a safe value.
        pipette.move_to(dst.bshm(z= max(z_offset_dst + z_dst_in_mm*0.8, 0)))

        # Adjust the volume to be mixed
        mix_kwargs['volume'] = min(mix_kwargs['volume'], pipette.max_volume)
        pipette.mix(**mix_kwargs)

    ## 6. Lift up, and blow out
    pipette.blow_out(dst.bshm(z=min(z_dst_in_mm * 1.5, max_blowout_height)))

    ## 7. Lift up to a certain height
    pipette.move_to(dst.top(z=default_z))

    ## 8. Drop the tip
    if drop_tip or mix:
        # If mixing was performed, then drop the tip to avoid contamination
        pipette.drop_tip()

    # Reset the blow-out rate
    reset_blow_out_rate(pipette, default_blow_out_flow_rate)

    # Home the pipettes
    if home:
        pipette.home()

    # Get the last location
    last_location = pipette._get_last_location_by_api_version()._labware
    return pipette, last_location


def dilute(protocol, pipettes, src, dst, vol, diluent, dilution_factor,
           transfer_src_kwargs=transfer_default_kwargs,
           transfer_dil_kwargs=transfer_diluent_kwargs):
    """
    Function to dilute a sample

    Parameters
    ----------
    protocol: opentrons.protocol_api.protocol_context.ProtocolContext
    pipettes: list of opentrons.protocol_api.instrument_context.InstrumentContext objects
    src: opentrons.protocol_api.labware.Labware object, e.g.- plate[0]['A1']
    dst: opentrons.protocol_api.labware.Labware object  e.g.- plate[0]['A1']
    vol: float, volume of the liquid to be transferred in uL
    diluent: opentrons.protocol_api.labware.Labware object, e.g.- stock_tube_racks[0]['B4']
    dilution_factor: float, dilution factor (must be greater than 1)
    transfer_src_kwargs: dictionary, kwargs for the transfer function, ffault=transfer_default_kwargs
    transfer_dil_kwargs: dictionary, kwargs for the transfer function for the diluent, default=transfer_diluent_kwargs

    Returns
    -------
    pipette: opentrons.protocol_api.instrument_context.InstrumentContext object
    last_location: opentrons.protocol_api.labware.Labware object
        ... The labware where the last transfer was made
    """
    vol_diluent = vol * (dilution_factor - 1)
    # Dilution is a two-step process
    pipette, _ = transfer(protocol, pipettes, src, dst, vol, **transfer_src_kwargs)
    drop_tip(pipette)

    pipette, last_location = transfer(protocol, pipettes, diluent, dst, vol_diluent, **transfer_dil_kwargs)
    drop_tip(pipette)

    return pipette, last_location




# HELPERS
def sub_volumes(pipette, vol):
    """
    Function to split the volume into sub-volumes that can be handled by the pipette

    Parameters
    ----------
    pipette: opentrons.protocol_api.instrument_context.InstrumentContext object
    vol: float, volume to be split

    Returns
    -------
    vol_list: list of floats, list of volumes that can be handled by the pipette
    """
    vol_max = pipette.max_volume
    vol_list = []
    if vol > vol_max:
        while vol > vol_max:
            vol_list.append(vol_max)
            vol -= vol_max
        if vol < pipette.min_volume:
            vol_last = (vol_list[-1] + vol) / 2
            vol_list[-1] = vol_last
            vol_list.append(vol_last)
        else:
            vol_list.append(vol)

    else:
        vol_list.append(vol)
    return vol_list