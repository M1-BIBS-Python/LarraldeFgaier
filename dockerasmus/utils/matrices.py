# coding: utf-8
from __future__ import unicode_literals
from __future__ import absolute_import

import numpy


def compose(function, x, y):
    """Create an array A: i,j --> func(x[i], y[j]) from vectors x and y

    About 5x quicker than using `numpy.fromfunction` where the function
    accesses x[i] and y[j] for each (i,j).
    """
    return function(x.reshape(1, *x.shape).repeat(len(y), 0).T, y)


def distance(u, v):
    """Create a matrix of euclidian distance between two sets of points.

    Given a matrix u of shape (n, 3) and a matrix v of dimension (m, 3),
    returns a matrix d of shape (n, m) where d[i,j] is the euclidian
    distance between the i-th point of u and the j-th point of v.

    Exemple:
        >>> x = numpy.array([ (0, 0), (1, 1) ])
        >>> y = numpy.array([ (4, 3), (0, 1) ])
        >>> distance(x,y)[0, 1]   # (0,0) / (0,1)
        1.0
        >>> distance(x,y)[1, 0]   # (1,1) / (4,3)
        3.6055512754639891
        >>> distance(x,y)
        array([[ 5.   ,  1.   ],
               [ 3.606,  1.   ]])

    """

    # Compute distance for each component (x, y for a 2D space,
    # x,y,z for a 3D space, etc.) using the compose function
    # and then the
    d_components = [
        compose(lambda x,y: (x-y)**2, u_k, v_k)
            for u_k, v_k in zip(u.T, v.T)
    ]

    return numpy.sqrt(sum(d_components))