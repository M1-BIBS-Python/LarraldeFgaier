# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import numpy


def compose(function, x, y):
    """Create an array A: i,j --> func(x[i], y[j]) from vectors x and y.

    About 5x quicker than using `numpy.fromfunction` where the function
    accesses x[i] and y[j] for each (i,j).
    """
    return function(x.reshape(1, *x.shape).repeat(len(y), 0).T, y)


def distance(u, v):
    """Create a matrix of euclidean distance between two sets of points.

    Given a matrix u of shape (n, 3) and a matrix v of dimension (m, 3),
    returns a matrix d of shape (n, m) where d[i,j] is the euclidian
    distance between the i-th point of u and the j-th point of v.

    Arguments:
        u (`numpy.ndarray`): a vector of positions (abcissa of
            generated distance matrix), shape :math:`(m, k)`.
        v (`numpy.ndarray`): a vector of positions (ordinate of
            generated distance matrix), shape :math:`(n, k)`.

    Returns:
        `numpy.ndarray`: a distance matrix :math:`d` of dimensions
        :math:`(m, n)` where :math:`d[i,j]` is the euclidean
        distance between the i-th point of ``u`` and the j-th point
        of ``v``.

    Example:
        >>> x = numpy.array([ (0, 0), (1, 1) ])
        >>> y = numpy.array([ (4, 3), (0, 1) ])
        >>> distance(x,y)[0, 1]   # (0,0) / (0,1)
        1.0
        >>> distance(x,y)[1, 0]   # (1,1) / (4,3)
        3.6055512754639891
        >>> distance(x,y)
        array([[ 5.   ,  1.   ],
               [ 3.606,  1.   ]])

    Note:
        This method is a slow equivalent of `scipy.spatial.distance.cdist`,
        but only depends on `numpy`, and only works for 2D arrays (vectors
        of spatial positions).
    """
    d_components = [
        compose(lambda x,y: (x-y)**2, u_k, v_k)
            for u_k, v_k in zip(u.T, v.T)
    ]
    return numpy.sqrt(sum(d_components))


def normalized(a, axis=-1, order=2):
    """Return an array of normalized vectors.

    Arguments:
        a (`numpy.ndarray`): an array of vectors

    Keyword Arguments:
        axis (`int`): the axis along which to norm.
        order (`int`): the order of the norm.

    Example:
        >>> a = numpy.array([ (0, 0), (3, 3), (0, 2) ])
        >>> normalized(a)
        array([[ 0.   ,  0.   ],
               [ 0.707,  0.707],
               [ 0.   ,  1.   ]])
    """
    l2 = numpy.atleast_1d(numpy.linalg.norm(a, order, axis))
    l2[l2==0] = 1
    return a / numpy.expand_dims(l2, axis)
