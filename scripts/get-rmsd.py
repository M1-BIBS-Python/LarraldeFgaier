#!/usr/bin/env python
# coding: utf-8
"""
Usage:
    get-rmsd.py -r RECEPTOR -l LIGAND -t TEST [options]

Positional Arguments:
    -r RECEPTOR         The path to the PDB file of
                        the native receptor conformation.
    -l LIGAND           The path to the PDB file of
                        the native ligand conformation.
    -t TEST             The path to the PDB file of
                        the tested ligand conformation.

Optional Arguments:
    -h, --help              Print this message.
    -i INTERFACE_THRESHOLD  The distance under which
                            residues are considered to be
                            at the contact interface.
"""
from __future__ import print_function

# stdlib imports
import sys
import os

# update sys.path to make dockerasmus importable
# locally although it is in the parent directory
SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))
MAINDIR = os.path.dirname(SCRIPTDIR)
sys.path.insert(0, MAINDIR)

# Try importing non standard dependencies
try:
    import docopt
except ImportError:
    sys.exit('Could not import docopt - is it installed ?')

from dockerasmus.pdb import Protein

if __name__ == "__main__":

    args = docopt.docopt(__doc__)

    # Import PDB files into Protein objects
    receptor = Protein.from_pdb_file(args['-r'])
    ligand = Protein.from_pdb_file(args['-l'])
    test = Protein.from_pdb_file(args['-t'])

    # Create two Proteins containing only the residues
    # in the Receptor/Ligand interface
    interface_ligand = ligand.copy()
    interface_test = test.copy()

    print(
        "Number of atoms in Ligand: ",
        len(list(test.iteratoms()))
    )

    print(
        "Number of atoms in Receptor: ",
        len(list(receptor.iteratoms()))
    )

    # Remove the residues not in the interface from the copies
    interface_residues = {res1 for res1, _ in ligand.interface(receptor)}
    for chain in ligand.itervalues():
        for res in chain.itervalues():
            if res not in interface_residues:
                del interface_test[chain.id][res.id]
                del interface_ligand[chain.id][res.id]

    print(
        "Number of atoms in Ligand (interface only): ",
        len(list(interface_test.iteratoms()))
    )

    # Compute and display RMSD
    print(
        "            Ligand-only RMSD: ",
        ligand.rmsd(test)
    )
    print(
        "Receptor/Ligand complex RMSD: ",
        (receptor+ligand).rmsd(receptor+test)
    )
    print(
        "  Interface only Ligand RMSD: ",
        interface_ligand.rmsd(interface_test)
    )
