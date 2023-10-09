# Creating host guest systems in batch and writing Gaussian input files

## Workflow

1. _(For the first time only)_ Create virtual environment (`venv`) from `create_host_guest.yml`:

    ```
    conda env create -f create_host_guest.yml
    ```

1. Activate `venv`:

    ```
    conda activate create_host_guest
    ```

1. To run script, enter the directory path to the hosts and either 1) filepath to a guest or 2) the SMILES string of the guest:

    ```
    python create_host_guest.py -i <dirpath_to_hosts> -g <filepath_to_guest>

     python create_host_guest.py -i <dirpath_to_hosts> -g <smiles_of_guest>
    ```

1. For help:

    ```
    python create_host_guest.py -h
    ```

## Test run

1. Run:

    ```
    python create_host_guest.py -i data/test_run/hosts/ -g '[Cl-]' --multiplicity 2
    ```

## Features in progress

1. Method to avoid atomic overlaps.
2. Calculate multiplicity automatically.