# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import division

import theano


def norm(a, axis=-1):
    """Return an array of norms from of an array of vectors.
    """
    return theano.tensor.sqrt(
        theano.tensor.sum(theano.tensor.sqr(a), axis=-1, keepdims=True)
    )


def normalized(a, axis=-1):
    """Return an array of normalized vectors.
    """
    return a / norm(a, axis=axis)


def distance(a,b):
    """Return a matrix of distances from a matrix of coordinates.
    """
    mx_arr_a = theano.tensor.repeat(a, b.shape[0], axis=0)
    mx_arr_b = theano.tensor.tile(b, (a.shape[0], 1))
    # Calculate vector of euclidian distance
    v_d = theano.tensor.sqrt(theano.tensor.sum((mx_arr_a-mx_arr_b)**2, axis=1))
    # Rearrange vector as a matrix where mx_d[i,j] is the distance between
    # the i-th atom of protein 1 and the j-th atom of protein 2
    return theano.tensor.reshape(v_d, (a.shape[0], b.shape[0]))
