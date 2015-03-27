# Sebastian Raschka 03/26/2015
#
# Tool to extract protein residues or atoms within a specified radius around a
# ligand's heavy atoms based on a template protein-ligand complex.
#

import os

class Pdb():
    """ Object that allows operations with protein files in PDB format. """

    def __init__(self):
        self.cont = []
        
    def read_pdbfile(self, file_path):
        with open(file_path, 'r') as pdb_file:
            self.cont = [row.strip() for row in pdb_file.read().split('\n') if row.strip()]
        self._process_pdb()
        return self

    def _process_pdb(self):
        if self.cont:
             self.atom = [row for row in self.cont if row.startswith('ATOM')]
             self.atom_ter = [row for row in self.cont if row.startswith(('ATOM', 'TER'))]
             self.hetatm = [row for row in self.cont if row.startswith('HETATM')]
             self.mainchain = [row for row in self.atom if  row[13:15] in ('CA', 'N ', 'C ', 'O ')]
             self.calpha = [row for row in self.mainchain if row[13:15] == 'CA']
             self.heavy_atom = [row for row in self.atom if 'H' not in row[12:16]]
             self.heavy_hetatm = [row for row in self.hetatm if 'H' not in row[12:16]]
             
    def grab_radius(self, radius, coordinates):
        """
        Grabs those atoms that are within a specified
        radius given a 3D-coordinate.

        Parameters
        ----------
        
        radius : `int` or `float`.
          Radius in Angstroms.
          
        Coordinates : `list` 
          A list of x, y, z coordinates , e.g., `[1.0, 2.4, 4.0]`

        Returns
        ----------

        atom_cont : `list`.
          List of PDB file contents that are within the specified radius.

        """
        in_radius = []
        for line in self.atom + self.hetatm:
            xyz_coords = self._get_xyz_coords(line)
            distance = (sum([(coordinates[i]-xyz_coords[i])**2 for i in range(3)]))**0.5
            if distance <= radius:
                in_radius.append(line)
        return in_radius  
        
    def _get_xyz_coords(self, pdb_line):
        return [float(pdb_line[30:38]), float(pdb_line[38:46]), float(pdb_line[46:54])]

    def _map_residues(self):
        residue_map = dict()
        for l in self.atom + self.hetatm:
            if not l[17:26] in residue_map:
                residue_map[l[17:26]] = [l]
            else:
                residue_map[l[17:26]].append(l)
        return residue_map



if __name__ == '__main__':
    
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Extracts a protein substructure based on the residues or atoms surrounding the ligand\'s heavy atoms.',
        epilog='Example:\n'\
                'bindingsite_scoop_templated.py -p my_struct.pdb -l my_struct.pdb -t target.pdb -o substructure.pdb\n',
        formatter_class=argparse.RawTextHelpFormatter
        )

    parser.add_argument('-p', '--template_protein', type=str, metavar='PDB', help='Template Protein, PDB.')
    parser.add_argument('-l', '--template_ligand', type=str, metavar='PDB', help='Template Ligand, PDB.', required=True)
    parser.add_argument('-t', '--target_structure', type=str, metavar='PDB', help='Target structure.', required=True)
    parser.add_argument('-o', '--output_pdb', type=str, metavar='PDB', help="Path to the output PDB file", required=True)
    parser.add_argument('-a', '--apply_to_dir', action='store_true', default=False, help="Treats --output_pdb and --target_structure as directories (default: False)")
    parser.add_argument('-e', '--extract', type=str, default='residues', help='"atoms" or "residues" within radius (default="residues")')
    parser.add_argument('-r', '--radius', type=float, default=9.0, help="Radius around ligand's heavy atoms (default=9.0)")

    args = parser.parse_args()  
    
    if args.apply_to_dir:
        target_structures = [os.path.join(args.target_structure, f) for f in os.listdir(args.target_structure)]
        if not os.path.isdir(args.output_pdb):
            os.mkdir(args.output_pdb)
        output_pdbs = [os.path.join(args.output_pdb, f) for f in os.listdir(args.target_structure)]

    else:
        template_proteins = [args.template_protein]
        template_ligands = [args.template_ligand]
        target_structures = [args.target_structure]
        output_pdbs = [args.output_pdb]

    if args.extract not in ('atoms', 'residues'):
        raise ValueError('--extract must be "atoms" or "residues"')

    # Load template structures
    temp_prot = Pdb().read_pdbfile(args.template_protein)
    temp_lig = Pdb().read_pdbfile(args.template_ligand)    
    
    # get ligand coordinates and template atoms
    lig_coords = [temp_lig._get_xyz_coords(l) for l in temp_lig.heavy_hetatm]
    template_atoms = set()
    for c in lig_coords:
        template_atoms = template_atoms.union({l[17:26] for l in temp_prot.grab_radius(radius=args.radius, coordinates=c)})

    # apply to target structures
    for ts, output in zip(target_structures, output_pdbs):
        
        tar_struc = Pdb().read_pdbfile(ts) 
        target_atoms = [l for l in tar_struc.atom + tar_struc.hetatm if l[17:26] in template_atoms]

        # residues within radius
        if args.extract == 'residues':
            mapped_residues = tar_struc._map_residues()
            residues = []
            for l in target_atoms:
                residues.extend(mapped_residues[l[17:26]])

            target_atoms = set(residues)
        
        with open(output, 'w') as out_file:
            for line in sorted(target_atoms):
                out_file.write(line+'\n')
