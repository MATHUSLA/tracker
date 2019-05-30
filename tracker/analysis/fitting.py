# -*- coding: utf-8 -*- #
#
# tracker/analysis/fitting.py
#
#
# MIT License
#
# Copyright (c) 2018-2019 Brandon Gomes
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

"""
MATHUSLA Tracker Analysis Library.

"""

# -------------- Standard Library -------------- #

from contextlib import contextmanager

# -------------- External Library -------------- #

from iminuit import Minuit
from uptrack.fitting import ParameterEnumBase, Fitter

# -------------- Tracker  Library -------------- #

from ..util.math import nominal_value, std_dev


__all__ = (
    'CartesianParameterType',
    'MinuitFitter',
)


class CartesianParameterType(ParameterEnumBase):
    """"""

    @classmethod
    def from_cartesian(cls, coordinate):
        """Get Parameter from Cartesian Coordinate."""
        return cls(coordinate.value)


class MinuitFitter(Fitter):
    """
    MINUIT Fitter.

    """

    @classmethod
    def parameters_to_dict(cls, parameters):
        """Convert ParameterSet into Dictionary for Insertion into MINUIT."""
        parameter_dict = {}
        for name, value in parameters:
            name = name.lower()
            parameter_dict[name] = nominal_value(value)
            parameter_dict['error_' + name] = std_dev(value)
            if isinstance(value, Bounded):
                parameter_dict['limit_' + name] = (nominal_value(value.minimum),
                                                   nominal_value(value.maximum))
        return parameter_dict

    @classmethod
    def from_minuit_parameters(cls, parameters):
        """Convert from MINUIT Parameters to Parameter Set."""
        return NotImplemented

    def __init__(self, function, *, errordef=0.5, print_level=0, max_iterations=10000, strategy=2):
        """Initialize Minuit Fitter."""
        self.function = function
        self.errordef = errordef
        self.print_level = print_level
        self.max_iterations = max_iterations
        self.strategy = strategy

    def fit(self, parameters, data):
        """Perform MINUIT Fit."""
        parameter_dict = self.parameters_to_dict(parameters)
        self.minuit = Minuit(partial(self.function, data),
                             errordef=self.errordef,
                             print_level=self.print_level,
                             **parameter_dict)
        self.minuit.set_strategy(self.strategy)
        minimum, out_parameters = self.minuit.migrad(ncall=self.max_iterations)
        return (minimum.is_valid,
                self.from_minuit_parameters(out_parameters),
                self.minuit.covariance,
                {'minimum': minimum, 'minuit': self.minuit})
