# 2026_Matsuzawa_Susceptibility

**Author:** Takumi Matsuzawa
**Last updated:** April 27, 2026
**Contact** tm688_at_cornell.edu or takumi.matz_at_gmail.com

## 1) Description

This repository contains scripts and a notebook used in:

> Matsuzawa et al., bioRxiv (2026), *Metabolites Shift Equilibria of Biomolecular Condensates*.

- bioRxiv link: https://www.biorxiv.org/content/10.64898/2026.01.14.699531v2
- The codebase is designed to automate preparation of phase-separated samples using the Opentrons OT-2 liquid handling robot.

## 2) High-level workflow

1. Load and calibrate OT-2 hardware (pipettes, tip racks, plates, tube racks).
2. Prepare concentrated stock solutions.
3. Create a sample list CSV (recommended via `master.ipynb`).
4. Generate protocol files from `master.ipynb`.
5. Copy generated files to OT-2 (e.g., with SSH/SCP).
6. Log in to OT-2 and run `startProtocol.py`.
7. Monitor the run on OT-2 to ensure safe operation.

## 3) System specification used in this study

### OT-2 / Opentrons
- Firmware version: `v1.1.0-25e5ceaa`
- Supported Protocol API versions: `v2.0` to `v2.17`
- Robot server version: `v7.2.2`
- Protocol API version used: `2.13`

### Control computer
- OS: Windows 11
- OT-2 connection: USB A/B cable

## 4) Detailed setup and operation

### 4.1 Hardware setup

1. Attach 2 µL or 20 µL pipettes to OT-2 (follow Opentrons setup wizard).
2. Run positional/deck calibration in Opentrons software.
   Reference: https://support.opentrons.com/s/article/Get-started-Calibrate-the-deck
3. Load labware in deck slots (tip racks, plates, tube racks).
4. Register loaded labware in `preparation/labware_configuration.csv`.

#### `preparation/labware_configuration.csv` notes

- `type` must be one of: `plate`, `pipette`, `tiprack`, `stock`
- `name` must be an Opentrons-recognized labware/instrument name (e.g., `p300_single_gen2`, `thermoscientificnunc_96_wellplate_2000ul`)
- `deck location` must be in OT-2 deck positions `1`-`11`
- For `pipette`, `deck location` must be `left` or `right`
- `starting_tip` is the first tip position to use (e.g., `A1`, `A2`, `B1`)
- For `stock`, `chemical_name` must match entries in `preparation/stocks.csv`
- Enter stock volumes for each solution
- For Falcon tubes, use a supported tuberack (e.g., `opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical`) and set:
  - `rack_location`
  - `tube_type` (`falcon15` or `falcon50`)

### 4.2 Software setup

1. Create the Conda environment from `environment.yml`:

```bash
cd /path_to_parent_folder/2026_Matsuzawa_Susceptibility
conda env create -f environment.yml
```

2. Generate SSH keys for OT-2 access.
   Reference: https://support.opentrons.com/s/article/Setting-up-SSH-access-to-your-OT-2

> ⚠️ Never share your SSH private key.

3. Open `master.ipynb`.
4. Follow notebook instructions to create/read sample lists.
5. Submit the generated protocol job to OT-2.

## 5) What `master.ipynb` does

The notebook follows this sequence:

1. Create/read sample information (including `preparation/stocks.csv`)
2. Create `Sample` objects
3. Create a labware configuration template
4. **User edits** the labware configuration file (`preparation/labware_configuration.csv`)
5. Read labware configuration
6. Assign unique well positions to samples
7. (Optional) Create sample sheet (sample name, well, composition, etc.)
8. (Optional) Simulate protocol
9. Upload generated protocol to OT-2
10. Run the protocol

Typical remote commands used are:

```bash
ssh -i ~/.ssh/ot2_ssh_key root@<OT2_IP_ADDRESS>
cd /var/lib/jupyter/notebooks/2026_Matsuzawa_Susceptibility && opentrons_execute ./startProtocol.py &
```

## 6) Repository files (quick reference)
- `preparation/`: Input CSVs for samples, stocks, and labware configuration
- 'samples/': Location for the created sample list. Each sample list is stored in a pickle. For manually created sample recipes, place the sample sheet under 'samples/csv/'. 
- `master.ipynb`: Main workflow notebook
- `startProtocol.py`: Protocol entry point executed on OT-2
- `sample_handling.py`, `rotbot_control.py`: Core helper scripts
- `environment.yml`: Reproducible Conda environment definition