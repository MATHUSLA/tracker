# -*- coding: utf-8 -*- #
#
# tracker/geometry.py
#
#
# MIT License
#
# Copyright (c) 2018 Brandon Gomes
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
MATHUSLA Tracking Geometry.

"""

# -------------- Standard Library -------------- #

# -------------- External Library -------------- #

from uptrack.geometry import Geometry as GeometryBase
from uptrack.geometry import Volume as VolumeBase

# -------------- Tracker  Library -------------- #


class GEANT4Geometry(GeometryBase):
    """
    GEANT4 Base Geometry Object.

    """

    @classmethod
    def from_file(cls, path):
        """Load Geometry from File."""

    def save(self, path, *args, **kwargs):
        """Save Geometry to File."""

    def __getitem__(self, key):
        """Get Volume by Key."""

    def __setitem__(self, key, value):
        """Set Volume by Key."""

    def __delitem__(self, key):
        """Delete Volume by Key."""

    def __iter__(self):
        """Get Volume Iterator."""

    def __len__(self):
        """Get Number of Volumes."""

    def volume_around(self, point):
        """Get Volume Around a Point."""

    def is_inside_volume(self, point, key):
        """Check If Point is Inside Volume."""


class GEANT4Volume(VolumeBase):
    """
    GEANT4 Base Volume Object.

    """

    @property
    def center(self):
        """Get Geometric Center of Volume."""

    @property
    def time_resolution(self):
        """Get Time Resolution of Volume."""

    @property
    def bounding_box(self):
        """Get Bounding Box around Volume."""
