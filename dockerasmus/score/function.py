# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import logging

from . import requirements
from .components.base import BaseComponent


class ScoringFunction(object):
    """The generalisation of a scoring function

    The ScoringFunction will first fetch dynamically the
    requirements that must be pre-computed before calling each
    scoring component: for instance, since both `components.Coulomb`
    and `components.LennardJones` require the atom-wise distance
    matrix between the two proteins, this matrix will be
    created only once. This *magic* behaviour is obtained by
    matching the names of the arguments of each scoring component
    with an actual function in the `requirements` module.

    Attributes:
        components (`list`): the list of individual components
            that are used independently to compute the final
            score.
        requirements (`set`): a set of arguments that must be
            preprocessed to compute the score, based on the
            requirements of the individual scoring components.

    Examples:

        Non-bound terms of Cornell's scoring function:

        >>> f = ScoringFunction(LennardJones, Coulomb)
        >>> f(barnase, barstar)
        -84.9...

        Using weights on individual components:

        >>> g = ScoringFunction(LennardJones, Fabiola, weights=[1, 3])
        >>> g(barnase, barstar)
        -118.54...
    """

    def __init__(self, *components, **kwargs):
        self.components = []
        self.weights = kwargs.get('weights') or [1 for _ in range(len(components))]
        for component in components:
            if issubclass(component, BaseComponent):
                logging.debug("Creating new {} instance...".format(component.__name__))
                self.components.append(component())
            elif isinstance(component, BaseComponent):
                logging.debug("Registering {} instance...".format(component.__class__.__name__))
                self.components.append(component)
            else:
                raise TypeError("Invalid component: {}".format(component))
        self.requirements = {
            req:getattr(requirements, req)
                for c in components for req in c.args()
        }

    def __call__(self, protein1, protein2, **parameters):
        requirements = self._compute_requirements(protein1, protein2)
        score = 0
        for weight, component in zip(self.weights, self.components):
            args = self._filter_requirements(component, requirements)
            kwargs = self._filter_parameters(component, parameters)
            score += weight*component(*args, **kwargs)
        return score

    def _compute_requirements(self, protein1, protein2):
        return {req: func(protein1, protein2) for req,func in self.requirements.items()}

    @staticmethod
    def _filter_requirements(component, requirements):
        return [requirements[arg] for arg in component.args()]

    @staticmethod
    def _filter_parameters(component, parameters):
        return {k:v for k,v in parameters.items() if k in component.kwargs()}
