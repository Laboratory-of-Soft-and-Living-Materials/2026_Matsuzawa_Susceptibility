# 2026_Matsuzawa_Susceptibility

**Author:** Takumi Matsuzawa
**Last updated:** July 15, 2026
**Contact:** `tm688_at_cornell.edu` or `takumi.matz_at_gmail.com`

## 1. Description

This repository contains scripts and a Jupyter notebook used in:

> Matsuzawa et al. (2026), *Susceptibility and Regulation of Biomolecular
Condensates by Solutes*, PNAS (accepted).

* **Preprint:** https://www.biorxiv.org/content/10.64898/2026.01.14.699531v3
* The codebase is designed to automate the preparation of phase-separated samples using the Opentrons OT-2 liquid-handling robot.

## 2. High-Level Workflow

1. Load and calibrate the OT-2 hardware, including pipettes, tip racks, plates, and tube racks.
2. Prepare concentrated stock solutions.
3. Create a sample-list CSV file, preferably using `master.ipynb`.
4. Generate protocol files from `master.ipynb`.
5. Copy the generated files to the OT-2, for example using SSH or SCP.
6. Log in to the OT-2 and run `startProtocol.py`.
7. Monitor the OT-2 during the run to ensure safe operation.

## 3. System Specifications Used in This Study

### 3.1 OT-2 / Opentrons

* **Firmware version:** `v1.1.0-25e5ceaa`
* **Supported Protocol API versions:** `v2.0` to `v2.17`
* **Robot server version:** `v7.2.2`
* **Protocol API version used:** `2.13`

### 3.2 Control Computer

* **Operating system:** Windows 11
* **OT-2 connection:** USB A-to-B cable

## 4. Detailed Setup and Operation

### 4.1 Hardware Setup

1. Attach 2 µL or 20 µL pipettes to the OT-2 by following the Opentrons setup wizard.
2. Run positional and deck calibration using the Opentrons software.
   See the [Opentrons deck calibration guide](https://support.opentrons.com/s/article/Get-started-Calibrate-the-deck).
3. Load the required labware into the OT-2 deck slots, including tip racks, plates, and tube racks.
4. Register the loaded labware in `preparation/labware_configuration.csv`.

#### Notes for `preparation/labware_configuration.csv`

* `type` must be one of the following:

  * `plate`
  * `pipette`
  * `tiprack`
  * `stock`
* `name` must be an Opentrons-recognized labware or instrument name, such as:

  * `p300_single_gen2`
  * `thermoscientificnunc_96_wellplate_2000ul`
* `deck location` must correspond to an OT-2 deck position from `1` to `11`.
* For a `pipette`, `deck location` must be either `left` or `right`.
* `starting_tip` specifies the first tip position to use, such as `A1`, `A2`, or `B1`.
* For a `stock`, `chemical_name` must match an entry in `preparation/stocks.csv`.
* Enter the stock volume for each solution.
* For Falcon tubes, use a supported tube rack, such as `opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical`, and specify:

  * `rack_location`
  * `tube_type`, which must be either `falcon15` or `falcon50`

### 4.2 Software Setup

1. Create the Conda environment from `environment.yml`:

   ```bash
   cd /path_to_parent_folder/2026_Matsuzawa_Susceptibility
   conda env create -f environment.yml
   ```

2. Generate SSH keys for OT-2 access.
   See the [Opentrons SSH setup guide](https://support.opentrons.com/s/article/Setting-up-SSH-access-to-your-OT-2).

   > [!WARNING]
   > Never share your SSH private key.

3. Open `master.ipynb`.

4. Follow the instructions in the notebook to create or load a sample list.

5. Submit the generated protocol job to the OT-2.

## 5. What `master.ipynb` Does

The notebook follows the sequence below:

1. Creates or reads sample information, including `preparation/stocks.csv`.
2. Creates `Sample` objects.
3. Creates a labware-configuration template.
4. Prompts the user to edit `preparation/labware_configuration.csv`.
5. Reads the completed labware configuration.
6. Assigns a unique well position to each sample.
7. Optionally creates a sample sheet containing sample names, well positions, compositions, and related information.
8. Optionally simulates the protocol.
9. Uploads the generated protocol to the OT-2.
10. Runs the protocol.

Typical remote commands include:

```bash
ssh -i ~/.ssh/ot2_ssh_key root@<OT2_IP_ADDRESS>
cd /var/lib/jupyter/notebooks/2026_Matsuzawa_Susceptibility
opentrons_execute ./startProtocol.py &
```

## 6. Repository Files

| File or directory    | Description                                                                              |
| -------------------- | ---------------------------------------------------------------------------------------- |
| `preparation/`       | Input CSV files for samples, stocks, and labware configuration                           |
| `samples/`           | Storage location for generated sample lists; each sample list is stored as a pickle file |
| `samples/csv/`       | Location for manually created sample-recipe CSV files                                    |
| `master.ipynb`       | Main workflow notebook                                                                   |
| `startProtocol.py`   | Protocol entry point executed on the OT-2                                                |
| `sample_handling.py` | Core sample-handling utilities                                                           |
| `robot_control.py`   | Core robot-control utilities                                                             |
| `environment.yml`    | Reproducible Conda environment definition                                                |
