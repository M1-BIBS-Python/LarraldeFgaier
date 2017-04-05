# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import division

from .base import BaseComponent
from .. import requirements


class LennardJones(BaseComponent):
    """A scoring component modeling Van der Waals interactions.

    Reference:
        `Jones, J. E.
        "On the Determination of Molecular Fields. II. From the Equation of
        State of a Gas". Proceedings of the Royal Society of London A:
        Mathematical, Physical and Engineering Sciences 106, 463–477 (1924).
        <https://dx.doi.org/10.1098%2Frspa.1924.0082>`_
    """

    backends = ["theano", "numpy"]

    def _setup_theano(self, theano):
        ### Epsilon matrix from protein vectors
        v_eps1, v_eps2 = theano.tensor.dvectors('v_eps1', 'v_eps2')
        mx_eps = theano.tensor.sqrt(theano.tensor.outer(v_eps1, v_eps2))
        ### Radius matrix from protein vectors
        v_rad1, v_rad2 = theano.tensor.dvectors('v_rad1', 'v_rad2')
        mx_rad = theano.tensor.transpose(theano.tensor.repeat(
            theano.tensor.reshape(
                v_rad1, (1, v_rad1.shape[0])
            ), v_rad2.shape[0], axis=0,
        )) + v_rad2
        ### Atomwise distance matrix
        mx_distance = theano.tensor.dmatrix('mx_distance')
        ### Van der Waals constants
        mx_A = mx_eps * mx_rad**12
        mx_B = 2 * mx_eps * mx_rad**6
        ### Final function
        self._call = theano.function(
            [v_eps1, v_eps2, v_rad1, v_rad2, mx_distance],
            theano.tensor.sum(theano.tensor.triu(
                mx_A/(mx_distance**12) - mx_B/(mx_distance**6), k=1,
            ))
        )

    def _setup_numpy(self, numpy):
        from ...utils.matrices import compose
        def call(v_eps1, v_eps2, v_rad1, v_rad2, mx_distance):
            ### Epsilon matrix from protein vectors
            mx_epsilon = numpy.sqrt(numpy.outer(v_eps1, v_eps2))
            ### Radius matrix from protein vectors (to the power of 6)
            mx_radius_6 = compose(lambda r1, r2: r1+r2, v_rad1, v_rad2)**6
            ### Van der Waals constants
            mx_B = 2 * mx_epsilon * mx_radius_6
            mx_A = 0.5 * mx_B * mx_radius_6
            ### Atomwise distance matrix
            mx_distance_6 = mx_distance**6
            ### Final function
            return numpy.sum(numpy.triu(
                mx_A/(mx_distance_6**2)-mx_B/(mx_distance_6), k=1,
            ))
        self._call = call

    def __call__(self, epsilon, vdw_radius, distance):
        return self._call(
            epsilon[0], epsilon[1],
            vdw_radius[0], vdw_radius[1],
            distance,
        )
