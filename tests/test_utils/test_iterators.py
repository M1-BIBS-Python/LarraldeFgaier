# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import unittest

from dockerasmus.utils import iterators


class TestInfiniteWords(unittest.TestCase):

    def test_default_behaviour(self):
        # wraps after the end of the alphabet is reached
        expected = [u"ZX", u"ZY", u"ZZ", u"AAA", u"AAB"]
        it = iterators.wordrange(u"ZX", u"AAC")
        self.assertEqual(expected, list(it))

    def test_wrong_order(self):
        # no iterator if start > stop in alphabetic order
        self.assertEqual(list(iterators.wordrange("B", "A")), [])


class TestNth(unittest.TestCase):

    def test_behaviour(self):
        self.assertEqual(iterators.nth(range(10), 2), 2)
        self.assertEqual(iterators.nth(range(2, 10), 3), 5)

    def test_default(self):
        self.assertIs(iterators.nth(range(10), 11), None)
        self.assertEqual(iterators.nth(range(10), 11, 'default'), 'default')
