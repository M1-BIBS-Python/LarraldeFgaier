# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import logging

from . import components
from . import requirements


class ScoringFunction(object):
    """The generalisation of a scoring function

    Attributes:
        components (list): the list of individual components
            that are used independently to compute the final
            score.
        requirements (set): a set of arguments that must be
            preprocessed to compute the score, based on the
            requirements of the individual scoring components.

    Examples:

        Non-bound terms of Cornell's scoring function:

        >>> f = ScoringFunction(LennardJones, Coulomb)
        >>> f(barnase, barstar)
        -22.4...

        Using weights on individual components:

        >>> g = ScoringFunction(LennardJones, Fabiola, weights=[1, 3])
        >>> g(barnase, barstar)
        -49.1...
    """

    def __init__(self, *components, **kwargs):
        self.components = []
        self.weights = kwargs.get('weights') or [1 for _ in range(len(components))]
        for component in components:
            logging.debug("Initializing {}".format(component.__name__))
            self.components.append(component())
        self.requirements = {req for c in components for req in c.requires}
        self.parameters = {param for c in components for param in c.parameters}

    def __call__(self, protein1, protein2, **parameters):
        requirements = self._compute_requirements(protein1, protein2)
        score = 0
        for weight, component in zip(self.weights, self.components):
            args = self._filter_requirements(component, requirements)
            kwargs = self._filter_parameters(component, parameters)
            score += weight*component(*args, **kwargs)
        return score

    def _compute_requirements(self, protein1, protein2):
        return {req: req(protein1, protein2) for req in self.requirements}

    @staticmethod
    def _filter_requirements(component, requirements):
        return [requirements[req] for req in component.requires]

    @staticmethod
    def _filter_parameters(component, parameters):
        return {k:v for k,v in parameters.items() if k in component.parameters}
