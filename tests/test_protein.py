# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import unittest

from pdb import Protein
from .utils import DATADIR


class TestProperties(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.arginine_prot = Protein.from_pdb_file(
            os.path.join(DATADIR, 'arginine.pdb')
        )

    def test_centre_de_masse(self):
        pass
