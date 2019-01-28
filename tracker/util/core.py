# -*- coding: utf-8 -*- #
#
# tracker/util/core.py
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
Core Utilities.

"""

# -------------- Standard Library -------------- #

import operator

# -------------- UpTrack  Library -------------- #


def between(min_value, value, max_value, partial_order=operator.le):
    """Check if Value is Inbetween Min and Max."""
    return partial_order(min_value, value) and partial_order(value, max_value)


def value_or(value, default):
    """Return Value or Default if Value is None."""
    return value if value is not None else default


def instance_of(types):
    """Returns Function which Checks Type."""
    def inner(obj):
        return isinstance(obj, types)
    return inner


def subclass_of(types):
    """Returns Function which Checks Subclass."""
    def inner(obj):
        return issubclass(obj, types)
    return inner


class classproperty(property):
    """Class Property."""

    def __get__(self, obj, objtype=None):
        """Wrap Getter Function."""
        return super().__get__(objtype)

    def __set__(self, obj, value):
        """Wrap Setter Function."""
        return super().__set__(type(obj), value)

    def __delete__(self, obj):
        """Wrap Deleter Function."""
        super().__delete__(type(obj))
