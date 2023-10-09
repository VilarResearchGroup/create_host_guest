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

# RDKit
from rdkit import Chem 
from rdkit.Chem import AllChem as rdkit
from collections import defaultdict
from rdkit.Chem import rdFMCS
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdDistGeom
IPythonConsole.ipython_3d = True

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


# Handle command-line argument processing
parser = argparse.ArgumentParser(
    description = 'Description in progress'
    )
parser.add_argument('--input_csv', '-i',
                        help = '''Input CSV containing an entry ID (column index = 0) and
                        SMILES strings (column name = 'smiles').''')
parser.add_argument('--output_dir', '-o',
                        default = 'output/schro_sdf/',
                        help = '''Name of directory for output SDF files to be saved.
                        CAUTION: If directory already exists, the contents will be deleted.''')
parser.add_argument('--sdf_method',
                        default = 'schrodinger',
                        help = '''Method used to convert SMILES strings to SDF files.
                        Options:
                            schrodinger:    Schrodinger API
                            obabel:         Open Babel
                        Schrodinger API recommended for complex macrocyclic and cage structures.''')
parser.add_argument('--index_col',
                        default = 'receptor_id',
                        help = 'Column name for entry index in input csv.')
parser.add_argument('--smi_col',
                        default = 'smiles',
                        help = 'Column name for SMILES in input csv.')
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
    if not args.input_csv:
        parser.error('''Please provide an input CSV as an arguement:\n
    $ python {} -i <input_csv_filepath>
                     '''.format(__file__))

    # # Ensure input csv string ends with a '/'
    # smi_df = args.input_csv
    # smi_df = smi_df + '/' if not smi_df.endswith('/') else smi_df

    # Dir variables
    data_dir = 'data/'
    output_dir = 'output/'
    xyz_dir = output_dir + 'obabel_xyz/'

    # Generate SDF files from SMILES
    print('Using {} to convert SMILES strings to SDF files.'.format(args.sdf_method))
    if args.sdf_method == 'schrodinger':
        sdf_dir = gen_schrodinger_sdf(args.input_csv, args.output_dir, args.sdf_method, args.index_col, args.smi_col)
        gjf_dir = create_dir(output_dir, 'schro_gjf/')
    elif args.sdf_method == 'obabel':
        input_csv = pd.read_csv(args.input_csv, index_col=args.index_col)[args.smi_col]
        sdf_dir = gen_obabel_sdf(input_csv, output_dir=output_dir)
        gjf_dir = create_dir(output_dir, 'obabel_gjf/')

    # Write Gaussian input (.gjf) files from XYZ files
    for xyz in glob.glob(output_dir+'*_xyz/*.xyz'):
        mol_id = Path(xyz).stem
        charge = get_mol_charge(sdf_dir+mol_id+'.sdf')
        write_gaussian_input(
            xyz_file=xyz_dir+mol_id+'.xyz',
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


def write_gaussian_input(xyz_file, gjf_dr, job_type='optimisation', functional='b3lyp',
                         basis_set='6-311++g(d,p)',charge=0, multiplicity=1,
                         nprocshared=32, mem='100gb'):
    """"""
    
    if job_type == 'optimisation':
        job = 'opt'
    else:
        raise ValueError('Please enter a valid option:\noptimisation\n<to_be_added>')
    
    xyz_cord = []
    with open(xyz_file) as f:
        for line in f.readlines()[2:]:
            xyz_cord.append(line)
    
    filename = Path(xyz_file).stem
    gjf_file = gjf_dr + filename + '.gjf'
    title = filename + '_' + job
    
    try:
        os.remove(gjf_file)
    except OSError:
        pass        
    
    with open(gjf_file, 'a') as f:
        f.write('%nprocshared={}\n'.format(str(nprocshared)))
        f.write('%mem={}\n'.format(mem))
        f.write('%chk={}F.chk\n'.format(xyz_file[-14:-4]))
        f.write('# {} {}/{}\n'.format(job, functional, basis_set))
        # f.write('# {} {}/{} geom=connectivity\n'.format(job, functional, basis_set))
        f.write('\n')
        f.write('{}\n'.format(title))
        f.write('\n')
        f.write('{} {}\n'.format(charge, multiplicity))
        for line in xyz_cord:
            f.write(str(line))
        f.write('\n')


def get_mol_charge(sdf_file):
    """Return charge (int) of a SDF file.
    
    Parameter:
    sdf_file (str): filepath to SDF file.
    
    Return:
    charge (int): charge of an inputted molecule.
    """
    
    mol = Chem.SDMolSupplier(sdf_file)[0]
    charge = Chem.rdmolops.GetFormalCharge(mol)
    
    return charge


def gen_schrodinger_sdf(input_csv, output_dir, sdf_method, index_col, smi_col):

    # Loop over all 
    os.system('bash run_schrodinger_api.sh -s gen_schrodinger_sdf.py -i {} -o {} -m {} -x {} -c {}'.format(input_csv, output_dir, sdf_method, index_col, smi_col))

    sdf_dir = output_dir

    return sdf_dir


def gen_obabel_sdf(input_csv, output_dir):

    sdf_dir = create_dir(output_dir, 'obabel_sdf/')
    xyz_dir = create_dir(output_dir, 'obabel_xyz/')

    err = []

    try:
        # Setup tqdm tracking bar
        pbar = tqdm(total=len(input_csv))

        for mol_id, smi in input_csv.iteritems():
            
            pbar.set_description('Converting {0}.smi to {0}.sdf'.format(mol_id))

            # Read SMILES string.
            mol = pybel.readstring('smi', smi)
            
            # Set entry name in SDF file.
            mol.title = mol_id
            
            # Generate 3D coordinates.
            mol.make3D()
            
            # Write output SDF file.
            output_file = open(sdf_dir+mol_id+'.sdf', "w")
            output_file.write(mol.write('sdf'))
            
            # Write output XYZ file.
            output_file = open(xyz_dir+mol_id+'.xyz', "w")
            output_file.write(mol.write('xyz'))

            pbar.update()
        
        pbar.close()

    except OSError:
        print('Could not convert a SMILES string of {}. Use another method.'.format(mol_id))

    # Check which SMILES could not be converted.
    output_sdf = []
    for sdf in glob.glob(sdf_dir+'*.sdf'):
        basename = os.path.basename(sdf)
        mol_id = os.path.splitext(basename)[0]
        output_sdf.append(mol_id)

    for mol_id in input_csv.index:
        if mol_id not in output_sdf:
            err.append(mol_id)

    sdf_count = len(os.listdir(sdf_dir))
    print('Successfully converted {} SDF files from {} inputted SMILES strings.'.format(sdf_count, len(input_csv)))
    print('Could not convert the following SMILES strings:')
    for smi in err: print(smi)

    return sdf_dir

def create_dir(tail_dir, head_dir):
    
    new_dir = tail_dir + head_dir
    if os.path.exists(new_dir) and os.path.isdir(new_dir):
        shutil.rmtree(new_dir)
        os.makedirs(new_dir)
        print('Replaced old', new_dir)
    else:
        os.makedirs(new_dir)
        print('Created', new_dir)
    
    return new_dir


if __name__ == '__main__':
    main()
