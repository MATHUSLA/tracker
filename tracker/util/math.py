# -*- coding: utf-8 -*- #
#
# tracker/util/math.py
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
MATHUSLA Tracker Math Utilities.

"""

# -------------- External Library -------------- #

from aenum import Enum
# import numpy as np
# import sympy as sym
# import glorentz as glo
# import uncertainties.umath as umath
# from uncertainties import ufloat, unumpy

# -------------- Tracker  Library -------------- #

# ThreeVector, FourVector, GL3_BASIS, GL4_BASIS = glo.GL3GL4


class CartesianCoordinate(Enum):
    """"""

    T = 0
    X = 1
    Y = 2
    Z = 3
