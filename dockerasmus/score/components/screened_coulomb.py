# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

from .base import BaseComponent


class ScreenedCoulomb(BaseComponent):
    """A scoring component modeling screened electrostatic forces.

    Reference:
        `Mehler, E. L., and G. Eichele.
        “Electrostatic Effects in Water-Accessible Regions of Proteins.”
        Biochemistry 23, no. 17 (August 1, 1984): 3887–91.
        doi:10.1021/bi00312a015. <http://dx.doi.org/10.1021/bi00312a015>`_
    """
    backends = ["theano", "numpy"]

    def _setup_theano(self, theano):
        ### Dielectric constant
        diel = theano.tensor.dscalar('diel')
        ### Screening parameters
        A, k, l = theano.tensor.dscalars('A', 'k', 'l')
        B = diel - A
        ### Charge matrix from protein vectors
        v_q1, v_q2 = theano.tensor.dvectors('v_q1', 'v_q2')
        mx_q = theano.tensor.outer(v_q1, v_q2)
        ### Atomwise distance matrix
        mx_distance = theano.tensor.dmatrix('mx_distance')
        ### Effective dielectric permittivity
        mx_perm = A + B / (1 + k * theano.tensor.exp(-l * B * mx_distance))
        ### Final function
        self._call = theano.function(
            [v_q1, v_q2, mx_distance, diel, A, k, l],
            theano.tensor.sum(theano.tensor.triu(
                mx_q/(mx_perm*mx_distance)
            ))
        )

    def _setup_numpy(self, numpy):
        def call(v_q1, v_q2, mx_distance, diel, A, k, l):
            B = diel - A
            mx_perm = A + B / (1 + k * numpy.exp(-l * B * mx_distance))
            mx_q = numpy.outer(v_q1, v_q2)
            return numpy.sum(numpy.triu(mx_q/(mx_perm*mx_distance), 1))
        self._call = call

    def __call__(self, charge, distance, diel=65.0, A=-8.5525, k=7.7839, l=0.003627):
        return self._call(
            charge[0], charge[1], distance, diel, A, k, l,
        )
