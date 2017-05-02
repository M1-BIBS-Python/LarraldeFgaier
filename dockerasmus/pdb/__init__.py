# coding: utf-8
"""
pdb
===

The Object model of a PDB protein.
"""

from __future__ import absolute_import
from __future__ import unicode_literals

from .protein import Protein
from .residue import Residue
from .chain import Chain
from .atom import Atom

__author__ = "althonos"
__author_email__ = "martin.larralde@ens-cachan.fr"
__version__ = "0.1.0"
__license__ = "GPLv3"

__all__ = ["Protein", "Residue", "Chain", "Atom"]
