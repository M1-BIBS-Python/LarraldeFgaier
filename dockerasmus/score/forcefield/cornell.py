# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import warnings
import functools
import importlib

from ..base import ScoringFunction
from ...utils.matrices import distance


__all__ = ["score"]


class _CornellScoringFunction(ScoringFunction):
    """Score a docking configuration of the VDW term of Cornell's force field

    .. math::

       S = \sum_{i<j}{\\frac{A_{i,j}}{r_{i,j}^{12}}
           - 2*\\frac{B_{i,j}}{r_{i,j}^{6}}
           + \\frac{q_i q_j}{\epsilon * r_{i,j}}}


    Arguments:
        protein1 (Protein): the reference protein (receptor)
        protein2 (Protein): the moving protein (ligand)

    Keyword Arguments:
        diel (float): the value of the dielectric constant for the medium.
            Default value uses results obtained after simulations in
            `Cossins et al. <https://dx.doi.org/10.1371/journal.pcbi.1002066>`_.
            for the cytoplasm of E. Coli. [default: 65.0]

    Example:
        >>> cornell.score(barnase, barstar)
        -22.4...

    .. note::
        Available backends: `theano`, `numpy`

    .. seealso::
        Source concept of the forcefield described in
        `Cornell et al. <https://dx.doi.org/10.1021/ja00124a002>`_
    """
    _backends = ['theano', 'numpy']

    def _setup_theano(self, theano):
        from ..components.theano import lennard_jones, coulomb
        self._lennard_jones = lennard_jones
        self._coulomb = coulomb

    def _setup_numpy(self, numpy):
        from ..components.numpy import lennard_jones, coulomb
        self._lennard_jones = lennard_jones
        self._coulomb = coulomb

    def _score(self, v_eps1, v_eps2, v_rad1, v_rad2, v_q1, v_q2, mx_pos1, mx_pos2, diel):
        mx_distance = distance(mx_pos1, mx_pos2)
        return self._lennard_jones(v_eps1, v_eps2, v_rad1, v_rad2, mx_distance)\
             + self._coulomb(v_q1, v_q2, mx_distance, diel)

    def __call__(self, receptor, ligand, diel=65.0):
        return float(self._score(
            receptor.atom_epsilon(), ligand.atom_epsilon(),
            receptor.atom_radius(), ligand.atom_radius(),
            receptor.atom_charges(), ligand.atom_charges(),
            receptor.atom_positions(), ligand.atom_positions(),
            diel,
        ))


score = _CornellScoringFunction()
