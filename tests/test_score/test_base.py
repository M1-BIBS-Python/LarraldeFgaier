# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import sys
import unittest
import functools
import contextlib
import importlib
import warnings

from dockerasmus.score.base import ScoringFunction
from dockerasmus.pdb import Protein

from ..utils import mock


class TestScoringFunction(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        class TestFunction(ScoringFunction):
            _backends = ['math', 'fail', 'os']
            def _setup_fail(self, fail):
                self._score = lambda: 1
            def _setup_os(self, os):
                self._score = lambda: 2
            def _setup_math(self, math):
                self._score = lambda: 3
            def __call__(self):
                return self._score()

        class TestFunction2(ScoringFunction):
            _backends = ["fail"]
            def _setup_fail(self, fail):
                self._score = lambda x: x
            def __call__(self):
                return self._score()

        class TestFunction3(ScoringFunction):
            _backends = ["fail", "math"]
            def _setup_fail(self, fail):
                self._score = lambda: 1
            def _setup_math(self, math):
                self._score = math.sqrt
            def __call__(self, x):
                return self._score(x)


        cls.TestFunction = TestFunction
        cls.TestFunction2 = TestFunction2
        cls.TestFunction3 = TestFunction3

    def test_unimportable_backend(self):
        with self.assertRaises(ValueError) as ctx:
            func = self.TestFunction('fail')
        self.assertTrue(str(ctx.exception).startswith('Unavailable'))

    def test_unknown_backend(self):
        with self.assertRaises(ValueError) as ctx:
            func = self.TestFunction('unknown')
        self.assertTrue(str(ctx.exception).startswith('Unknown'))

    def test_backend_order(self):
        func = self.TestFunction()
        self.assertEqual(func._backend, 'math')

    def test_force_backend(self):

        func = self.TestFunction(force_backend='os')
        self.assertEqual(func._backend, 'os')
        self.assertEqual(func(), 2)

    def test_no_available_backend(self):
        with self.assertRaises(RuntimeError):
            func = self.TestFunction2()

    def test_backend_fallback(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('ignore')
            func = self.TestFunction3()
        self.assertEqual(func._backend, 'math')
        self.assertEqual(func(4), 2)
