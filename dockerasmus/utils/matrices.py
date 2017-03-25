# coding: utf-8
from __future__ import unicode_literals
from __future__ import absolute_import


def compose(function, x, y):
    """Create an array A: i,j --> func(x[i], y[j]) from vectors x and y

    About 5x quicker than using `numpy.fromfunction` where the function
    accesses x[i] and y[j] for each (i,j).
    """
    return function(x.reshape(1, *x.shape).repeat(len(y), 0).T, y)
