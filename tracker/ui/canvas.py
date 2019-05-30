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

# -------------- Standard Library -------------- #

from collections import defaultdict

# -------------- External Library -------------- #

from uptrack.ui import CanvasBase, NamedColor
from rootpy import ROOT
from rootpy.ROOT import TCanvas, TPolyLine3D, TPolyMarker3D, TView, TView3D, TColor
from rootpy.io import root_open

# -------------- Tracker  Library -------------- #

from ..util import classproperty


__all__ = (
    'ROOTCanvas',
    'MPLCanvas'
)


def _to_TColor_id(color):
    """Conver Color to ROOT Color."""
    return TColor.GetColor(*color.rgb)


class ROOTCanvas(CanvasBase):
    """
    ROOT Canvas.

    """

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
        super().__post_init__()
        self._canvas = TCanvas(self.name, self.title, self.width, self.height)
        self._view = TView.CreateView()
        self._point_map = defaultdict(list)
        self._line_set = set()
        self._has_updated = False

    def add_point(self, point, size=1, color=NamedColor.Black):
        """Add Point To Canvas."""
        self._point_set[(color, size)].append(point)

    def add_line(self, start, end, width=1, color=NamedColor.Black):
        """Add Line To Canvas."""
        poly_line = TPolyLine3D()
        poly_line.SetNextPoint(start.x, start.y, start.z)
        poly_line.SetNextPoint(end.x, end.y, end.z)
        poly_line.SetLineWidth(width)
        poly_line.SetLineColor(_to_TColor_id(color))

    @property
    def is_empty(self):
        """Check if Canvas is Empty."""
        return not self._point_set and not self._line_set

    def draw(self):
        """Draw Canvas."""
        self._canvas.cd()

        for (color, size), points in self._point_map.items():
            polymarker = TPolyMarker3D(len(points), 20)
            for point in points:
                polymarker.SetNextPoint(point.x, point.y, point.z)
            polymarker.SetMarkerColor(_to_TColor_id(color))
            polymarker.SetMarkerSize(size)
            polymarker.Draw()

        for poly_line in self._line_set:
            poly_line.Draw()

        if not self._has_updated:
            self._view.ShowAxis()
            self._canvas.cd()
            axis = TAxis3D.GetPadAxis()
            if axis:
                axis.SetLabelColor(ROOT.kBlack)
                axis.SetAxisColor(ROOT.kBlack)
                axis.SetTitleOffset(2)
                axis.SetXTitle('X ()')
                axis.SetYTitle('Y ()')
                axis.SetZTitle('Z ()')

        self._canvas.Modified()
        self._canvas.Update()
        self._has_updated = True

    def clear(self):
        """Clear the Canvas."""
        if self._has_updated:
            self._canvas.cd()
            self._canvas.Clear()
            self._canvas.Modified()
            self._canvas.Update()
            self._reset_view()
            self._point_set.clear()
            self._line_set.clear()
            self._has_updated = False

    def save(self, path):
        """Save Canvas to ROOT File."""
        with root_open(path, 'UPDATE') as file:
            file.cd()
            file.WriteTObject(self._canvas.Clone())

    def load(self, path, name):
        """Load Canvas from ROOT File."""
        return NotImplemented


class MPLCanvas(CanvasBase):
    """
    Matplotlib Canvas.

    """
