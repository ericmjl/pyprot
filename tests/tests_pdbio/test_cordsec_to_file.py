import pyprot

def test_atom_ter_hetatm_trio():
    
    in_pdb = "./tests/tests_pdbio/data/atom_ter_hetatm_in.pdb"
    out_pdb = "./tests/tests_pdbio/data/atom_ter_hetatm_out.pdb"
    expected_pdb = "./tests/tests_pdbio/data/atom_ter_hetatm_out_expected.pdb"
    
    pdb1 = pyprot.Pdb(in_pdb)
    pdb1.coordsec_to_file(out_pdb)
    
    with open(out_pdb, 'r') as in_1, open(expected_pdb, 'r') as in_2:
        for line1, line2 in zip(in_1, in_2):
            assert(line1 == line)