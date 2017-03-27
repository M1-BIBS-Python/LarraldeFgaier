# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest
import math

from dockerasmus import utils


class TestMaybeImport(unittest.TestCase):

    def test_missing_import(self):
        self.assertIs(utils.maybe_import('fail'), None)
        self.assertIs(utils.maybe_import('fail2'), None)

    def test_present_import(self):
        self.assertIs(utils.maybe_import('math'), math)
        self.assertIs(utils.maybe_import('unittest'), unittest)
