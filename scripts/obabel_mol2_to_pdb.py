import argparse
import pyprind
import subprocess
from mlxtend.file_io import find_files
import os

parser = argparse.ArgumentParser(
    description='Converts Ligand files from MOL2 to PDB files using OpenBabel.',

    formatter_class=argparse.RawTextHelpFormatter
    )


parser.add_argument('-i', '--input', type=str, help='Input directory with MOL2 files.')
parser.add_argument('-o', '--output', type=str, help='Output directory for PDB files.')
parser.add_argument('-e', '--extension', type=str, default='.pdb', help='Optional exension, e.g., "_lig.pdb" (default: ".pdb").')
parser.add_argument('-r', '--recursive', action='store_true', help='Search directories recursively.')


args = parser.parse_args()


if not args.input:
    print('Please provide an input directory via the -i flag. \n')
    args.help() 
    quit()

if not args.output:
    print('Please provide an input directory via the -o flag. \n')
    args.help() 
    quit()
    
if not os.path.isdir(args.output):
    os.mkdir(args.output)
    
mol2_files = find_files(substring='', 
               path=args.input, 
               recursive=args.recursive, 
               check_ext='.mol2', 
               ignore_invisible=True)
               
               
pbar = pyprind.ProgBar(len(mol2_files))

for path in mol2_files:
    
    mol2_name = os.path.basename(path)
    
    if args.recursive:
        case_dir = os.path.split(os.path.split(path)[0])[1]
        new_dir = os.path.join(args.output, case_dir)
    
    else:
        new_dir = args.output
    
    new_file = mol2_name.split('.mol2')[0] + args.extension
    new_path = os.path.join(new_dir, new_file)

    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)
    
    t = subprocess.call('obabel -i mol2 %s -o pdb > %s 2>/dev/null' % (path, new_path), shell=True)
    
    pbar.update()