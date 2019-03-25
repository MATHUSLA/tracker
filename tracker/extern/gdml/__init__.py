# -*- coding: utf-8 -*- #
#
# tracker/extern/gdml/__init__.py
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
GDML Parsing Library.

"""

# -------------- External Library -------------- #

from lxml import etree
from path import Path

# -------------- UpTrack  Library -------------- #

from .geometry import *


__all__ = (
    'GDML_SCHEMA',
    'GDMLStore',
    'parse_gdml'
)


GDML_SCHEMA = etree.XMLSchema(file=Path(__file__).parent / 'schema/gdml.xsd')


def _clear_return(element, *return_values):
    """Return Results after Clearing Element."""
    element.clear()
    if len(return_values) == 1:
        return return_values[0]
    return return_values


def parse_definitions(elem):
    """Parse Definitions."""
    return _clear_return(elem, list(elem))


def parse_materials(elem):
    """Parse Materials Definitions."""
    return _clear_return(elem, list(elem))


def parse_solids(elem):
    """Parse Solid Structure."""
    return _clear_return(elem, list(elem))


def parse_structure(elem):
    """Parse Geometry Structure."""
    return _clear_return(elem, list(elem))


def parse_setup(elem):
    """Parse World Volume Setup."""
    return _clear_return(elem, {'name': elem.attrib['name'],
                                'version': elem.attrib['version'],
                                'world_ref': elem[0].attrib['ref']})


class _GDMLFirstPass(etree.TreeBuilder):
    """"""

    def start(self, tag, attrib):
        """"""
        return super().start(tag, attrib)


def _default_gdml_parser():
    """Default GDML Parser with First Pass Target."""
    return etree.XMLParser(target=_GDMLFirstPass())


class GDMLStore:
    """GDML Storage Object."""

    @classmethod
    def from_file(cls, path):
        """Read GDML Store from File Path."""
        return cls(etree.parse(path, parser=_default_gdml_parser()))

    @classmethod
    def _parse_validated_tree(cls, tree):
        """Parse Validated Tree."""
        dictionary = {}
        for element in tree:
            if element.tag == 'define':
                dictionary['definitions'] = parse_definitions(element)
            elif element.tag == 'materials':
                dictionary['materials'] = parse_materials(element)
            elif element.tag == 'solids':
                dictionary['solids'] = parse_solids(element)
            elif element.tag == 'structure':
                dictionary['structure'] = parse_structure(element)
            elif element.tag == 'setup':
                dictionary['setup'] = parse_setup(element)
        return dictionary

    def __init__(self, tree, *, skip_validation=False):
        """Initalize GDML Store from Tree."""
        self._tree = tree
        if not skip_validation:
            GDML_SCHEMA.assertValid(tree)
        self.__dict__.update(type(self)._parse_validated_tree(tree))

    @property
    def tree(self):
        """Get XML Tree associated to GDML Store."""
        return self._tree

    @property
    def geometry(self):
        """Get Assoicated Geometry from GDML File."""
        return GDMLGeometry()

    def clear(self):
        """Clear GDML Store."""
        self.definitions.clear()
        self.materials.clear()
        self.solids.clear()
        self.structure.clear()
        self.setup.clear()


def parse_gdml(path):
    """Parse GDML File from Path."""
    return GDMLStore.from_file(path)
