# -*- coding: utf-8 -*- #
#
# tracker/ui/canvas.py
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
MATHUSLA Tracking Canvas.

"""

# -------------- External Library -------------- #

from uptrack.ui import CanvasBase
from ROOT import TCanvas

# -------------- Tracker  Library -------------- #


class ROOTCanvas(CanvasBase):
    """"""

    @classproperty
    def default_width(cls):
        """Default Width of Canvas."""
        return 700

    @classproperty
    def default_height(cls):
        """Default Height of Canvas."""
        return 500

    def __post_init__(self):
        """Post-Initialize ROOT Canvas."""
        self._canvas = TCanvas(self.name, self.title, self.width, self.height)
        self._point_set = set()
        self._line_set = set()

    def add_point(self, point, size=1, color=NamedColor.Black):
        """Add Point To Canvas."""
        return NotImplemented

    def add_line(self, start, end, width=1, color=NamedColor.Black):
        """Add Line To Canvas."""
        return NotImplemented

    @property
    @abstractmethod
    def is_empty(self):
        """Check if Canvas is Empty."""
        return not self._point_set and not self._line_set

    def draw(self):
        """"""
        return NotImplemented

    def clear(self):
        """"""
        return NotImplemented

    def save(self, path):
        """"""
        return NotImplemented

    def load(self, path, name):
        """"""
        return NotImplemented


class MPLCanvas(CanvasBase):
    """"""
