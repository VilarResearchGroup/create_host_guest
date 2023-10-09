# General
import sys
import os

# Error message for common issues.
err_msg = '''
Exited {}.

Try either of the following:
    
    - Ensure the virtual environment is activated.
    
        $ conda activate create_host_guest

    - Make sure to use python and not python3. Python 3.7.13 was used to write this script.
    
        $ python write_gaussian_input.py -i <path_to_input_csv>

    - Start a new command line session.
'''.format(os.path.basename(__file__))

try:
    import pandas as pd
except ModuleNotFoundError:
    sys.exit(err_msg)

import glob
import math
import sys
import shutil
import subprocess
import argparse
from tqdm import tqdm
from time import sleep
from pathlib import Path

# Open Babel
from openbabel import pybel

# stk
import py3Dmol
from IPython.display import Image
import matplotlib.pyplot as plt
import subprocess
import time
import stk
import stko
import spindry as spd
# %matplotlib inline

from rdkit import Chem 
from rdkit.Chem import AllChem as rdkit
from collections import defaultdict
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
IPythonConsole.ipython_3d = True

import py3Dmol
from IPython.display import Image
import matplotlib.pyplot as plt
import subprocess
import time
import stk
import stko
import spindry as spd
# %matplotlib inline

# Script from write_gaussian_input repo
import write_gaussian_input

# Handle command-line argument processing
parser = argparse.ArgumentParser(
    description = 'Description in progress'
    )
parser.add_argument('--host_dir', '-i',
                        default = 'data/',
                        help = 'Input directory containing SDF files')
parser.add_argument('--guest', '-g',
                        help = 'Provide the guest as either a 1) SMILES string, SDF filepath or SMI filepath.')
parser.add_argument('--output_dir', '-o',
                        default = 'output/stk_gjf/',
                        help = '''Name of directory for output GIF files to be saved.
                        CAUTION: If directory already exists, the contents will be deleted.''')
parser.add_argument('--job_type',
                        default = 'optimisation',
                        help = 'Gaussian job type. Current options:\n   optimisation: Geometry optimisation')
parser.add_argument('--functional',
                        default = 'b3lyp',
                        help = 'Functional. Current options:\n   b3lyp')
parser.add_argument('--basis_set',
                        default = '6-311++g(d,p)',
                        help = 'Basis set. Current options:\n   6-311++g(d,p)')
parser.add_argument('--multiplicity',
                        default = 1,
                        help = 'Multiplicity (int).')
parser.add_argument('--nprocshared',
                        default = 32,
                        help = 'nprocshared (int).')
parser.add_argument('--mem',
                        default = '100gb',
                        help = 'Memory (str).')


def main():

    print('Running {}.'.format(os.path.basename(__file__)))

    args = parser.parse_args()
    if not args.host_dir:
        parser.error('''Please provide a directory path to host SDF files as an arguement:\n
    $ python {} -i <host_dirpath>
                     '''.format(__file__))
    if not args.guest:
        parser.error('''Please provide a filepath to the guest SDF file or the SMILES string as an arguement:\n
    $ python {} -g <guest_filepath>
    $ python {} -g '<smiles_string>'
                     '''.format(__file__))

    # Dir variables
    data_dir = 'data/'
    output_dir = 'output/'
    xyz_dir = write_gaussian_input.create_dir(output_dir, 'stk_xyz/')
    mol_dir = write_gaussian_input.create_dir(output_dir, 'stk_mol/')
    gjf_dir = write_gaussian_input.create_dir(output_dir, 'stk_gjf/')

    # Load guest
    if '.sdf' or '.smi' in args.guest:
        guest = load_stk_guest(args.guest)
    elif isinstance(args.guest, str):
        guest = load_stk_guest(args.guest)
    else:
        parser.error('''Please provide a valid guest as either a .sdf, .smi or SMILES string:\n
    $ python {} -g <guest_filepath>
    $ python {} -g '<smiles_string>'
                     '''.format(__file__))

    # Loop over all hosts in host dir
    for sdf in glob.glob(args.host_dir+'*.sdf'):

        # Load host
        host = stk.BuildingBlock.init_from_file(sdf)

        # Extract host ID
        mol_id = Path(sdf).stem

        xyz_filepath = xyz_dir + mol_id + '.xyz'
        mol_filename = mol_dir + mol_id + '.mol'

        # Construct host guest complex
        complex = construct_complex(host, guest, xyz_filepath, mol_filename)

        # Define charge of overall complex
        charge = write_gaussian_input.get_mol_charge(mol_filename)
        
        # Write Gaussian input file from method in write_gaussian_input repo
        write_gaussian_input.write_gaussian_input(
            xyz_file=xyz_filepath,
            gjf_dr=gjf_dir,
            job_type=args.job_type,
            functional=args.functional,
            basis_set=args.basis_set,
            charge=charge,
            multiplicity=args.multiplicity,
            nprocshared=args.nprocshared,
            mem=args.mem,
        )

    gjf_count = len(os.listdir(gjf_dir))
    xyz_count = len(os.listdir(xyz_dir))
    print('Successfully converted {} GJF files from {} inputted XYZ files.'.format(gjf_count, xyz_count))

    print('Finished running {}.'.format(os.path.basename(__file__)))

    # Check if there are any atomic clashes
        # Translate or rotate guest until satisfied


def construct_complex(host, guest, xyz_filename, mol_filename):

    complex = stk.ConstructedMolecule(
        topology_graph=stk.host_guest.Complex(
            host=host,
            guests=guest,
        ),
    )

    writer = stk.XyzWriter()
    writer.write(molecule=complex, path=xyz_filename)

    writer = stk.MolWriter()
    writer.write(molecule=complex, path=mol_filename)

    return complex


def load_stk_bb(input):

    if '.sdf' in input:
        bb = stk.BuildingBlock.init_from_file(input)
    elif isinstance(input, str):
        try:
            bb = stk.BuildingBlock(smiles=input)
        except:
            print('Enter a valid guest SMILES string.')
            sys.exit()
    
    return bb


def load_stk_guest(input):

    if '.sdf' in input:
        guest = stk.host_guest.Guest(
            building_block=stk.BuildingBlock.init_from_file(input)
        )
    elif isinstance(input, str):
        try:
            guest = stk.host_guest.Guest(
                building_block=stk.BuildingBlock(smiles=input)
            )
        except:
            print('Enter a valid guest SMILES string.')
            sys.exit()

    return guest


if __name__ == '__main__':
    main()
