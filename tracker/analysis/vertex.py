# -*- coding: utf-8 -*- #
#
# tracker/analysis/vertex.py
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

from uptrack.vertex import Vertex as VertexBase
from uptrack.fitting import ParameterSetBase

# -------------- Tracker  Library -------------- #

from .fitting import CartesianParameterType, MinuitFitter
from ..util import value_or, partial


__all__ = (
    'vertex_track_r3_distance',
    'vertex_track_distances',
    'gaussian_nll',
    'DefaultFitter',
    'DEFAULT_FITTER',
    'VertexParameterSet',
    'Vertex'
)


def vertex_track_r3_distance(self, parameters, track):
    """"""
    track_point = track.at_t(parameters.t)
    # TODO: hypot function
    return hypot(track_point.x - parameters.x,
                 track_point.y - parameters.y,
                 track_point.z - parameters.z)

def vertex_track_distances(self, parameters, tracks):
    """"""
    for track in tracks:
        yield vertex_track_r3_distance(parameters, track)


def gaussian_nll(self, parameters, tracks):
    """"""
    # TODO: log function
    return sum(0.5 * distance.n ** 2.0 + log(distance.s)
               for distance in vertex_track_distances(parameters, tracks))


DefaultFitter = partial(MinuitFitter, gaussian_nll)
DEFAULT_FITTER = DefaultFitter()


class VertexParameter(CartesianParameterType):
    """"""
    T, X, Y, Z


class VertexParameterSet(ParameterSetBase, parameter_type=VertexParameter):
    """"""


class Vertex(VertexBase, parameter_set=VertexParameterSet, parameter_properties=True):
    """"""

    def __init__(self, tracks, fitter=None, *args, geometry=None, **kwargs):
        """"""
        super().__init__(self,
                         tracks,
                         value_or(fitter, DEFAULT_FITTER),
                         *args,
                         geometry=geometry,
                         **kwargs)

    @property
    def default_parameters(self):
        """Default Fit Parameters."""
        average_point = sum(track.at_t(track.t0) for track in self.tracks) / len(self.tracks)
        return VertexParameterSet(average_point.t,
                                  average_point.x,
                                  average_point.y,
                                  average_point.z)

    def divergence_response(self):
        """Respond to Divergences in Fit."""
        return NotImplemented

    @property
    def position(self):
        """Get Position of Vertex."""
        return FourVector(self.t, self.x, self.y, self.z)
