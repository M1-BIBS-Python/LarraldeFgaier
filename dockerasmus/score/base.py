# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import abc
import six

@six.add_metaclass(abc.ABCMeta)
class ScoringFunction(object):

    @abc.abstractmethod
    def __call__(self, protein1, protein2, **kwargs):
        pass
