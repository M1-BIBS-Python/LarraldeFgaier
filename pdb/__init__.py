# coding: utf-8
"""
pdb
===
An Object model of a PDB protein

Use `Protein.from_pdb_file` to quickly import a protein from a pdb file,
or manually create a protein (see)
"""

from __future__ import absolute_import
from __future__ import unicode_literals

from .protein import Protein
from .residual import Residual
from .chain import Chain
from .atom import Atom

__author__ = "althonos"
__author_email__ = "martin.larralde@ens-cachan.fr"
__version__ = "0.1.0"
__license__ = "GPLv3"

__all__ = ["Protein", "Residual", "Chain", "Atom"]
