![pyprot Logo](./images/logos/molecule_logo.png)

**PyProt is a Python package for working with protein structure files. It comes with a collection of ready-to-use scripts for the most common file operations and protein analyses.**



[[download pyprot.zip](https://github.com/rasbt/pyprot/archive/master.zip)] [[link to pyprot on GitHub](http://htmlpreview.github.io/?https://github.com/rasbt/pyprot/blob/master/README.html)]

<hr>
## ReadMe Contents

- [Scripts and command line tools](#scripts-and-command-line-tools)
- [Tutorials](#tutorials)
- [API documentation](#api-documentation)
- [Installation](#installation)

<hr>



<br>
<br>




## Scripts and command line tools

PyProt provides ready-to-use command line scripts that are using the underlying `pyprot` objects to work with PDB and MOL2 files.  
The scripts are located in the subdirectory `./scripts` and can be used after `pyprot` was successfully installed.   

### List of command line tools

- Working with PDB files
    - [Center of Mass](./docs/tools/pdb_center_of_mass.md)
    - [Grab atoms within a radius](./docs/tools/pdb_grab_atom_radius.md)
    - [Root-mean-square deviation (RMSD)](./docs/tools/pdb_rmsd.md)
    - [PDB to FASTA converter](./docs/tools/pdb_to_fasta.md)
    - [PDB atom and residue renumbering](./docs/tools/pdb_renumber.md)
    - [B-factor statistics](./docs/tools/pdb_bfactor_stats.md)
    - [PDB downloader](./docs/tools/pdb_download.md)
    - [Separating ATOM and HETATM sections](./docs/tools/pdb_split_atom_hetatm.md)
    - [Extracting Ligands from PDB files](./docs/tools/pdb_extract_ligands.md)

  
- Working with MOL2 files
    - [Transfer charges](./docs/tools/mol2_transfer_charge.md)
    - [Split multi-MOL2 files](./docs/tools/mol2_split.md)
    - [MOL2 functional group filter](./docs/tools/mol2_filter_funcgroups.md)
    - [MOL2 intermolecular functional group screening](./docs/tools/mol2_screening_intermol_funcgroup.md)

<br>
<br>


## API documentation

Please find the API documentation at [http://rasbt.github.io/pyprot/](http://rasbt.github.io/pyprot/).


<br>
<br>



## Tutorials

In progress ...

<br>
<br>

## Installation

PyProt was build and tested in Python 3.

The `pyprot` package can be installed  via 
	
	pip install pyprot
	
Alternatively, you can download `pyprot` directly from GitHub ([pyprot-master.zip](https://github.com/rasbt/pyprot/archive/master.zip)), unzip the archive, and install it via

	python setup.py install

After successful installation, the `pyprot` API is ready to use.

<hr>
	
If you are also interested in using the available scripts and tools you have two options:

**a)** If you downloaded `pyprot` from GitHub, you can use the scripts from the [`pyprot/scripts/`](./scripts) subdirectory.  

**b)** If  you installed `pyprot` via pip, you can download the the scripts bundle separately ([pyprot-scripts.zip](https://github.com/rasbt/pyprot/blob/master/scripts/pyprot-scripts.zip?raw=true)).


For more details, please see the separate [**Installation Documentation**](./docs/pyprot_installation.md).