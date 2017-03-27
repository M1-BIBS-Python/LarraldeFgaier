# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import abc
import six
import warnings

from .. import utils


@six.add_metaclass(abc.ABCMeta)
class ScoringFunction(object):

    _backends = []
    _backend = None

    def __init__(self, force_backend=None):

        # the actual scoring function
        self._score = None

        # try forcing a specific backend if required
        if force_backend is not None:
            if force_backend not in self._backends:
                raise ValueError("Unknown backend: {}".format(force_backend))
            backend_module = utils.maybe_import(force_backend)
            if backend_module is None:
                raise ValueError("Unavailable backend: {}".format(force_backend))
            else:
                backend = force_backend
                setup_function = getattr(self, "_setup_{}".format(force_backend))
                setup_function(backend_module)
                self._backend = force_backend

        # find the best available backend (first in the _backend list)
        else:
            unavailable_backends = []
            for backend in self._backends:
                backend_module = utils.maybe_import(backend)
                # the backend was imported
                if backend_module is not None:
                    setup_function = getattr(self, "_setup_{}".format(backend))
                    setup_function(backend_module)
                    self._backend = backend
                    if unavailable_backends: # if we are not using the best backend
                        warnings.warn("Unavailable backends: {}, using {}".format(
                            ', '.join(unavailable_backends), backend,
                        ), UserWarning)
                    break
                # the backend is unavailable
                else:
                    unavailable_backends.append(backend)

        # raise an error if no available backend was found.
        if self._score is None:
            raise RuntimeError(
                "Could not find any available backend for {} "
                "among: {}".format(
                    type(self).__name__, ', '.join(self._backends)
                )
            )


    @abc.abstractmethod
    def __call__(self, protein1, protein2, **kwargs):
        pass
