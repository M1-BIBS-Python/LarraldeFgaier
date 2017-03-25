# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals


import os

_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)), "{}.csv")
_transtype = lambda k,v: (k, float(v))


with open(_CSV.format('atomic_epsilon'), 'rb') as f:
    ATOMIC_EPSILON = dict(
        _transtype(*line.decode('utf-8').strip().split(','))
            for line in f
    )


with open(_CSV.format('atomic_masses'), 'rb') as f:
    ATOMIC_MASSES = dict(
        _transtype(*line.decode('utf-8').strip().split(','))
            for line in f
    )


with open(_CSV.format('atomic_radius'), 'rb') as f:
    ATOMIC_RADIUS = dict(
        _transtype(*line.decode('utf-8').strip().split(','))
            for line in f
    )


with open(_CSV.format('aminoacid_charges'), 'rb') as f:
    headers = next(f).decode('utf-8').split(',')[1:]
    AMINOACID_CHARGES = {
        line.split(',')[0]: {
            atom: float(charge)
                for atom, charge in zip(headers, line.split(',')[1:])
                    if charge
        }
        for line in (l.decode('utf-8').strip() for l in f)
    }


with open(_CSV.format('aminoacid_epsilon'), 'rb') as f:
    _headers = next(f).decode('utf-8').split(',')[1:]
    AMINOACID_EPSILON = {
        line.split(',')[0]: {
            atom: ATOMIC_EPSILON[ref]
                for atom, ref in zip(headers, line.split(',')[1:])
                    if ref in ATOMIC_EPSILON
        }
        for line in (l.decode('utf-8').strip() for l in f)
    }


with open(_CSV.format('aminoacid_radius'), 'rb') as f:
    _headers = next(f).decode('utf-8').split(',')[1:]
    AMINOACID_RADIUS = {
        line.split(',')[0]: {
            atom: ATOMIC_RADIUS[ref]
                for atom, ref in zip(headers, line.split(',')[1:])
                    if ref in ATOMIC_RADIUS
        }
        for line in (l.decode('utf-8').strip() for l in f)
    }


# Remove local variables and imports
del os, absolute_import, unicode_literals
del f, headers
