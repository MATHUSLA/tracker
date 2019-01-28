# -*- coding: utf-8 -*- #
#
# tracker/analysis/track.py
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

# -------------- External Library -------------- #

from uptrack.track import Track as TrackBase
from uptrack.fitting import Fitter, ParameterEnumBase, ParameterSetBase

# -------------- Tracker  Library -------------- #

from ..util import value_or
from ..util.math import umath, ThreeVector, FourVector, CartesianCoordinate


class TrackFitter(Fitter):
    """
    Base TrackFitter Object.

    """

    def track_squared_residual(self, parameters, point):
        """"""
        dt = (point.z.n - parameters.z0.n) / parameters.vz.n
        t_residual = (dt + parameters.t0.n - point.t.n) / point.t.s
        x_residual = (dt * parameters.vx.n + parameters.x0.n - point.x.n) / point.x.s
        y_residual = (dt * parameters.vy.n + parameters.y0.n - point.y.n) / point.y.s
        return t_residual ** 2.0 + 12.0 * x_residual ** 2.0 + 12.0 * y_residual ** 2.0

    def gaussian_nll(self, parameters, points):
        """"""
        return 0.5 * sum(self.track_squared_residual(parameters, point) for point in points)

    def fit(self, parameters, data):
        """Perform Fit."""


default_track_fitter = TrackFitter()


class Track(TrackBase):
    """
    MATHUSLA Track Object.

    """

    class Parameter(ParameterEnum):
        """"""
        T0 = 0
        X0, Y0, Z0, VX, VY, VZ

        @classmethod
        def from_cartesian(cls, coordinate):
            """Get Parameter from Cartesian Coordinates."""
            return cls(coordinate.value)

    class ParameterSet(ParameterSetBase, friend=Parameter):
        """"""

    def __init__(self, points, direction, fitter=None, *args, geometry=None, **kwargs):
        """Initialize Fitable Object."""
        super().__init__(self,
                         points,
                         value_or(fitter, default_track_fitter),
                         *args,
                         geometry=geometry,
                         **kwargs)
        self._fixed_parameter = set(Parameter.from_cartesian(direction))

    @property
    def default_parameters(self):
        """Default Fit Parameters."""
        first = self.front
        dr = self.back - first
        return first.t,
               first.x,
               first.y,
               first.z,
               dr.x / dr.t,
               dr.y / dr.t,
               dr.z / dr.t

    @property
    def parameters(self):
        """Get Parameters of Fit."""
        return Parameter.total_set

    @parameters.setter
    def parameters(self, parameters):
        """Set Parameters of Fit."""
        return NotImplemented

    @property
    def free_parameters(self):
        """Get Free Parameters of Fit."""
        return self.parameters - self._fixed_parameter

    @free_parameters.setter
    def free_parameters(self, free):
        """Set Free Parameters of Fit."""
        return NotImplemented

    @property
    def fixed_parameter(self):
        """Get Fixed Parameter of Fit."""
        return self._fixed_parameter

    @fixed_parameter.setter
    def fixed_parameter(self, fixed):
        """Set Fixed Parameter of Fit."""
        self._fixed_parameter = set(fixed)

    @property
    def fixed_parameters(self):
        return self.fixed_parameters

    @fixed_parameters.setter
    def fixed_parameters(self, fixed):
        self.fixed_parameter = fixed

    def divergence_response(self):
        """Respond to Divergences in Fit."""
        return NotImplemented

    @property
    def t0(self):
        """"""
        return self.fit_parameters.t0

    @property
    def x0(self):
        """"""
        return self.fit_parameters.x0

    @property
    def y0(self):
        """"""
        return self.fit_parameters.y0

    @property
    def z0(self):
        """"""
        return self.fit_parameters.z0

    @property
    def vx(self):
        """"""
        return self.fit_parameters.vx

    @property
    def vy(self):
        """"""
        return self.fit_parameters.vy

    @property
    def vz(self):
        """"""
        return self.fit_parameters.vz

    @property
    def origin(self):
        """Get Origin Point of Track."""
        return FourVector(self.t0, self.x0, self.y0, self.z0)

    @property
    def ray(self):
        """Get Vector Along Track."""
        return ThreeVector(self.vx, self.vy, self.vz)

    @property
    def unit(self):
        """"""
        return self.ray.unit

    @property
    def angle(self):
        """"""
        if Parameter.T0 in self._fixed_parameter:
            return 0  #FIXME: how to handle T direction parameterization?
        elif Parameter.X0 in self._fixed_parameter:
            return umath.acos(self.unit.x)
        elif Parameter.Y0 in self._fixed_parameter:
            return umath.acos(self.unit.y)
        elif Parameter.Z0 in self._fixed_parameter:
            return umath.acos(self.unit.z)
        else:
            raise TypeError('BAD!')

    def position(self, parameter):
        """Get Position of Track at Parameter."""
        if Parameter.T0 in self._fixed_parameter:
            dt = t - self.t0
        elif Parameter.X0 in self._fixed_parameter:
            dt = (parameter - self.x0) / self.vx
        elif Parameter.Y0 in self._fixed_parameter:
            dt = (parameter - self.y0) / self.vy
        elif Parameter.Z0 in self._fixed_parameter:
            dt = (parameter - self.z0) / self.vz
        else:
            raise TypeError('BAD!')
        return self.origin + dt * FourVector(1, self.vx, self.vy, self.vz)
