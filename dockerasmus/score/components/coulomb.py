# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

from .base import BaseComponent
from .. import requirements


class Coulomb(BaseComponent):
    """A scoring component modeling electrostatic forces.

    Reference:
        `Coulomb, C. A.
        "Premier mémoire sur l’électricité et le magnétisme". Histoire de
        l’Académie Royale des Sciences, 569-577 (1785).
        <https://books.google.com/books?id=by5EAAAAcAAJ&pg=PA569>`_
    """
    backends = ["theano", "numpy"]
    requires = [
        requirements.charge,
        requirements.distance
    ]
    parameters = {
        'diel',
    }

    def _setup_theano(self, theano):
        ### Dielectric constant
        diel = theano.tensor.dscalar('diel')
        ### Charge matrix from protein vectors
        v_q1, v_q2 = theano.tensor.dvectors('v_q1', 'v_q2')
        mx_q = theano.tensor.outer(v_q1, v_q2)
        ### Atomwise distance matrix
        mx_distance = theano.tensor.dmatrix('mx_distance')
        ### Final function
        self._call = theano.function(
            [v_q1, v_q2, mx_distance, diel],
            theano.tensor.sum(theano.tensor.triu(
                mx_q/(diel*mx_distance)
            ))
        )

    def _setup_numpy(self, numpy):
        def call(v_q1, v_q2, mx_distance, diel):
            mx_q = numpy.outer(v_q1, v_q2)
            return numpy.sum(numpy.triu(mx_q/(diel*mx_distance), 1))
        self._call = call

    def __call__(self, charge, distance, diel=65.0):
        return self._call(charge[0], charge[1], distance, diel)
