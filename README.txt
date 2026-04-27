Author: Takumi Matsuzawa
Date: 2026/04/27

1. Description:
- This repository stores codes and a notebook used in (Matsuzawa, et al., bioRxiv, 2026 "Metabolites Shift Equilibria of Biomolecular Condensates"). 
- The included codes were developed to automate sample preparation of phase-separated solutions  using OT-2 (Opentrons Labworks Inc.).
- Link to the paper: https://www.biorxiv.org/content/10.64898/2026.01.14.699531v2

2. Overview of autosample synthesis
2.1. Load labwares (pipettes, tips, etc.) to OT-2
2.2. Prepare concentrated solutions (we call them 'stock solutions').
2.3. Create a sample list in csv. (We recommend using our notebooks to generate a sample list.)
2.4. Copy the generated files to OT-2 through ssh and scp protocols. 
2.5. Follow 'master.ipynb' to generate a sample preparation protocol for OT-2.
2.5. SSH into OT-2.
2.6. Navigate yourself to '/var/lib/jupyter/notebooks/'.
2.7. Run '/var/lib/jupyter/notebooks/2026_Matsuzawa_Susceptibility/startProtocol.py'
2.8. Monitor OT-2 to make sure that nothing crashes. 


3.Spec:
3.1. OT-2, Opentrons
- Firmware Version: v1.1.0-25e5ceaa
- Supported Protocol API Versions: v2.0 - v2.17
- Robot Server Version: v7.2.2
- API version: 2.13
3.2. OS of the controlling PC: Windows 11
3.3. Connection between OT-2 and PC: A/B USB


4. How to auto-prepare using OT-2:
4.1. Setting up Hardware
4.1.1. Attach 2 uL or 20 uL pipettes to OT-2. Follow steps on the OT-2 Opentrons wizard
4.1.2. Run positional calibration through Opentrons native software for OT-2 (Instructions: https://support.opentrons.com/s/article/Get-started-Calibrate-the-deck?_gl=1*164of2d*_ga*MjAxOTE0OTk4My4xNzc3MzAwOTI1*_ga_66HK7MC5D7*czE3NzczMDA5MjQkbzEkZzAkdDE3NzczMDA5MjQkajYwJGwwJGgw*_gcl_aw*R0NMLjE3NzczMDA5MjcuQ2owS0NRandrcnpQQmhDcUFSSXNBSk40NjBsT0FUU3k0dVgtUXNZU2xnSGFCTURQdkFRMXE0M3h4c3RtU3hQdjFybmR0QmpiWFE4dzZ3RWFBalZFRUFMd193Y0I.*_gcl_au*MTEzMTcxNTgyNi4xNzc3MzAwOTI1*_ga_GNSMNLW4RY*czE3NzczMDA5MjQkbzEkZzAkdDE3NzczMDA5MjQkajYwJGwwJGgxMDE2Mjk0NzIy)
3.1.3. Load labwares in the setup (pipetting tips, well plates, and falcon tube holders)
3.1.4. Register loaded labwares to /preparation/labware_configuration.csv
... Type must be one of 'plate', 'pipette', 'tiprack', 'stock'
... 'name' must be registered in Opentrons. e.g. p300_single_gen2, thermoscientificnunc_96_wellplate_2000ul
... 'deck location' are the deck number in the machine. Range: 1 - 11.
... For a 'pipette', 'deck location' must be either 'left' or 'right'
... 'starting_tip' refers to the location of the first tip to be attached to the pipette. e.g. A1, A2, B1, etc.  
... For 'stock', enter the name of the chemical to 'chemical_name'. This chemical name must be identical to the ones in '/preparation/stock.csv'
... Enter volume of each stock solution.
... You must use a tuberack to place stock solutions in falcon tubes, such as 'opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical', enter the 'rack_location' (indicated on the rack you can purchase from Opentrons). Enter the 'tube_type' either as 'falcon15' or 'falcon50'.


4: Software Setup 
4.1 Set up an environment on your machine. The used conda environment for this study is stored in 'environment.yml'.
... Python environment
	To create an environment to use Otto, type the following commands on terminal.

   		'cd /path_to_parent_folder/otto_sample_prep'
   		'conda env create -f environment.yml'

4.2. Generate SSH keys for an SSH connection with OT-2
... Reference: https://support.opentrons.com/s/article/Setting-up-SSH-access-to-your-OT-2
... Never share your SSH private key with others
4.3. Open 'master.ipyb'
4.4. Follow the instructions to create or read a sample list.
4.5 Submit a job to OT-2


5. What master.ipynb does
# STEP 1: Creates a sample list (/preparation/stocks.csv)
# STEP 2: Create Sample objects
# STEP 3: Create a labware configuration file
# STEP 4: USER MUST EDIT THE FILE at /preparation/lab_configuration.csv
# STEP 5: Read the labware configuration file
# STEP 6: Assign a unique well position to each sample
# STEP 7: (Optional) Create a sample sheet (Sample name, well number, composition, etc.)
# STEP 8: (Optional) Simulate the protocol
# STEP 9: Upload the experimental protocol to the OT-2
# STEP 10: Run the protocol
... It calls subprocess to give following commands:

ssh -i ~/.ssh/ot2_ssh_key root@10.49.35.97
cd /var/lib/jupyter/notebooks/otto_sample_prep && opentrons_execute ./startProtocol.py &




