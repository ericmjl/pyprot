"""
Sebastian Raschka 2014

Unit tests for RMSD method in PdbStats class
from pyprot.pdbstats.rmsd

"""

import pyprot

lig1 = pyprot.Pdb("./tests/data/pdbs/lig_conf_1.pdb")
lig2 = pyprot.Pdb("./tests/data/pdbs/lig_conf_2.pdb")
pdb1 = pyprot.Pdb("./tests/data/pdbs/1T48_995.pdb")
pdb2 = pyprot.Pdb("./tests/data/pdbs/1T49_995.pdb")
pdb3 = pyprot.Pdb("./tests/data/pdbs/short_RIV_1.pdb")

def test_rmsd():
    assert(lig1.rmsd(lig1, ligand=True) == 0)
    assert(lig1.rmsd(lig2, ligand=True) == 1.9959)
    assert(lig1.rmsd(lig2, ligand=True, atoms = "all") == 2.6444)
    assert(pdb1.rmsd(pdb1) == 0.0)
    assert(pdb1.rmsd(pdb2) == 0.7377)
    assert(pdb1.rmsd(pdb2, atoms = "ca") == 0.4785)

def test_rmsd_old():
    assert(lig1.rmsd_(lig1, ligand=True) == 0)
    assert(lig1.rmsd_(lig2, ligand=True) == 1.9959)
    assert(lig1.rmsd_(lig2, ligand=True, atoms="all") == 2.6444)
    assert(pdb1.rmsd_(pdb1) == 0.0)
    assert(pdb1.rmsd_(pdb2) == 0.7377)
    assert(pdb1.rmsd_(pdb2, atoms = "ca") == 0.4785)

def test_rmsd_fail():
    assert(pdb1.rmsd(pdb3, ligand=False) == None)