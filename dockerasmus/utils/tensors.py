# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import theano
import numpy


def norm(a):
    return theano.tensor.sqrt(
        theano.tensor.sum(theano.tensor.sqr(a), axis=-1, keepdims=True)
    )


def normalized(a, axis=-1, order=2):
    return a / a.sum(axis=axis, keepdims=True)
    #return a / theano.tensor.nlinalg.norm(a)
    #return a / norm(a)
    #return a / a.sum(axis=-1, keepdims=True)#.reshape((a.shape[0], 1))

def distance(a,b):
    mx_arr_a = theano.tensor.repeat(a, b.shape[0], axis=0)
    mx_arr_b = theano.tensor.tile(b, (a.shape[0], 1))
    # Calculate vector of euclidian distance
    v_d = theano.tensor.sqrt(theano.tensor.sum((mx_arr_a-mx_arr_b)**2, axis=1))
    # Rearrange vector as a matrix where mx_d[i,j] is the distance between
    # the i-th atom of protein 1 and the j-th atom of protein 2
    return theano.tensor.reshape(v_d, (a.shape[0], b.shape[0]))

# def tcompose(function, x, y):
#     """Create an array A: i,j --> func(x[i], y[j]) from vectors x and y
#
#     About 5x quicker than using `numpy.fromfunction` where the function
#     accesses x[i] and y[j] for each (i,j).
#     """
#     return function(
#             x.reshape(1, *x.shape).repeat(len(y), 0).T, y
#
#     )
