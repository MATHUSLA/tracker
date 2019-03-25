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
from uptrack.fitting import ParameterSetBase

# -------------- Tracker  Library -------------- #

from .fitting import CartesianParameterType, MinuitFitter
from ..util import value_or, partial
from ..util.math import umath, ThreeVector, FourVector, CartesianCoordinate


__all__ = (
    'track_squared_residual',
    'gaussian_nll',
    'DefaultFitter',
    'DEFAULT_FITTER',
    'TrackParameter',
    'TrackParameterSet',
    'Track'
)


def track_squared_residual(self, t0, x0, y0, z0, vx, vy, vz, point):
    """Track Squared Residual Definition. (ONLY for Z0 Fixed)."""
    dt = (point.z.n - z0.n) / vz.n
    t_residual = (dt + t0.n - point.t.n) / point.t.s
    x_residual = (dt * vx.n + x0.n - point.x.n) / point.x.s
    y_residual = (dt * vy.n + y0.n - point.y.n) / point.y.s
    return t_residual ** 2.0 + 12.0 * x_residual ** 2.0 + 12.0 * y_residual ** 2.0


def gaussian_nll(self, points, t0, x0, y0, z0, vx, vy, vz):
    """Gaussian Negative Log Likelihood Track Fit."""
    residual = partial(track_squared_residual, t0, x0, y0, z0, vx, vy, vz)
    return 0.5 * sum(map(residual, points))


DefaultFitter = partial(MinuitFitter, gaussian_nll)
DEFAULT_FITTER = DefaultFitter()


class TrackParameter(CartesianParameterType):
    """"""
    T0, X0, Y0, Z0, VX, VY, VZ


class TrackParameterSet(ParameterSetBase, parameter_type=TrackParameter):
    """"""


class Track(TrackBase, parameter_set=TrackParameterSet, parameter_properties=True):
    """
    MATHUSLA Track Object.

    """

    def __init__(self, points, direction, fitter=None, *args, geometry=None, **kwargs):
        """Initialize Fitable Object."""
        super().__init__(self,
                         points,
                         value_or(fitter, DEFAULT_FITTER)
                         *args,
                         geometry=geometry,
                         **kwargs)
        self._fixed_parameter = Parameter.from_cartesian(direction)

    @property
    def default_parameters(self):
        """Default Fit Parameters."""
        first = self.front
        dr = self.back - first
        return TrackParameterSet(first.t,
                                 first.x,
                                 first.y,
                                 first.z,
                                 dr.x / dr.t,
                                 dr.y / dr.t,
                                 dr.z / dr.t)

    @property
    def default_parameter_freedoms(self):
        """"""
        freedoms = {parameter: True for parameter in self.parameters}
        freedoms[TrackParameter.Z0] = False
        return freedoms

    def divergence_response(self):
        """Respond to Divergences in Fit."""
        return NotImplemented

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
        if TrackParameter.T0 == self._fixed_parameter:
            return 0  #FIXME: how to handle T direction parameterization?
        elif TrackParameter.X0 == self._fixed_parameter:
            return umath.acos(self.unit.x)
        elif TrackParameter.Y0 == self._fixed_parameter:
            return umath.acos(self.unit.y)
        elif TrackParameter.Z0 == self._fixed_parameter:
            return umath.acos(self.unit.z)
        else:
            raise TypeError('BAD!')

    def position(self, parameter):
        """Get Position of Track at Parameter."""
        if TrackParameter.T0 == self._fixed_parameter:
            dt = t - self.t0
        elif TrackParameter.X0 == self._fixed_parameter:
            dt = (parameter - self.x0) / self.vx
        elif TrackParameter.Y0 == self._fixed_parameter:
            dt = (parameter - self.y0) / self.vy
        elif TrackParameter.Z0 == self._fixed_parameter:
            dt = (parameter - self.z0) / self.vz
        else:
            raise TypeError('BAD!')
        return self.origin + dt * FourVector(1, self.vx, self.vy, self.vz)
