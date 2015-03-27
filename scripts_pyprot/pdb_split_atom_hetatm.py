#!/usr/bin/env python

# Sebastian Raschka 2014
#
# Python PyProt script to autmatically create separate PDB files
# from ATOM and HETATM lines in a PDB file.
#
# run
# ./pdb_split_atom_hetatm.py -h
# for help
#

import argparse
import pyprot
import os

parser = argparse.ArgumentParser(
    description='Autmatically creates separate PDB files from ATOM and HETATM lines in a PDB file.',
    formatter_class=argparse.RawTextHelpFormatter
    )

parser.add_argument('-i', '--input', type=str, help='Input directory.', required=True)
parser.add_argument('-o', '--output', type=str, help='Output directory.', required=True)
parser.add_argument('-c', '--conect', action='store_true', help='Writes CONECT records to ligand file.')

args = parser.parse_args()

atom_out = os.path.join(args.output, 'pdb_atom')
hetatm_out = os.path.join(args.output, 'pdb_hetatm')

if not os.path.exists(args.output):
    os.mkdir(args.output)


for d in (atom_out, hetatm_out):
    if not os.path.exists(d):
        os.mkdir(d)

pdb_list = [os.path.join(args.input, pdb) for pdb in os.listdir(args.input)
            if pdb.endswith('.pdb')]



n = len(pdb_list)

if not n:
    print('{0}\nPDB list is empty. Please check the input directory for PDB files.\n{0}'.format(50* '-'))
    parser.print_help()
    quit()    


for pdb in pdb_list:
    pdb_obj = pyprot.Pdb(pdb)
    pdb_atom = pyprot.Pdb(pdb_obj.atom_ter)
    print(pdb_obj.conect)
    if args.conect:
        pdb_hetatm = pyprot.Pdb(pdb_obj.hetatm + pdb_obj.conect)
    else:
        pdb_hetatm = pyprot.Pdb(pdb_obj.conect)        

    pdb_atom.save_pdb(os.path.join(atom_out, os.path.basename(pdb)))
    pdb_hetatm.save_pdb(os.path.join(hetatm_out, os.path.basename(pdb)))


