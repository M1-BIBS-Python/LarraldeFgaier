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
    """

    # Compute distance for each component (x, y for a 2D space,
    # x,y,z for a 3D space, etc.) using the compose function
    # and then the
    # d_components = [
    #     compose(lambda x,y: (x-y)**2, u_k, v_k)
    #         for u_k, v_k in zip(u.T, v.T)
    # ]
    #
    # return numpy.sqrt(sum(d_components))
    mx_arr_u = numpy.repeat(u, v.shape[0], axis=0)
    mx_arr_v = numpy.tile(v, (u.shape[0], 1))
    # Calculate vector of euclidian distance
    v_d = numpy.sqrt(numpy.sum((mx_arr_u-mx_arr_v)**2, axis=1))
    # Rearrange vector as a matrix where mx_d[i,j] is the distance between
    # the i-th atom of protein 1 and the j-th atom of protein 2
    return numpy.reshape(v_d, (u.shape[0], v.shape[0]))


def normalized(a, axis=-1, order=2):
    """Returns an array of normalized vectors

    Arguments:
        a (numpy.array): an array of vectors

    Keyword Arguments:
        axis (int): the axis on which to norm **[default: -1]**
        order (int): the order of the norm **[default: 2]**

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
