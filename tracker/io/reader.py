# -*- coding: utf-8 -*- #
#
# tracker/io/reader.py
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
MATHUSLA Tracker Reader IO.

"""

# -------------- Standard Library -------------- #

from abc import ABC, abstractmethod
from collections.abc import Iterator
from contextlib import contextmanager

# -------------- External Library -------------- #

import uproot

# -------------- Tracker  Library -------------- #


__all__ = (
    'Traversal',
    'EventTraversal',
    'Reader',
    'EventReader',
    'ROOTEventReader'
)


class Traversal(Iterator):
    """"""

    def __init__(self, obj, *, collect=False, storage_factory=list):
        """"""
        self.obj = obj
        self.collect = collect
        self.store = storage_factory()
        self.started = False

    def clear(self):
        """"""
        self.store.clear()

    def reset(self):
        """"""
        self.clear()
        self.started = False

    @abstractmethod
    def start(self):
        """"""

    @abstractmethod
    def next_item(self):
        """"""

    def start_loop(self):
        """"""
        while not self.started:
            self.started = self.start()

    def __next__(self):
        """"""
        self.start_loop()
        item = self.next_item()
        if self.collect:
            self.store.append(item)
        return item

    def collect_all(self, *, reset_first=False):
        """"""
        if reset_first:
            self.reset()
        self.start_loop()
        while True:
            try:
                self.store.append(self.next_item())
            except StopIteration:
                break
        return self.store


class EventTraversal(Traversal):
    """"""


class Reader(ABC):
    """"""

    @abstractmethod
    def open(self, path, *args, **kwargs):
        """"""

    @abstractmethod
    def close(self):
        """"""

    @contextmanager
    def load(self, path, *args, **kwargs):
        """"""
        try:
            yield self.open(path, *args, **kwargs)
        finally:
            return self.close()


class EventReader(Reader):
    """"""


class ROOTEventReader(EventReader):
    """"""

    def open(self, path, name=None):
        """"""
        obj = uproot.open(path)
        if name:
            obj = obj[name]
        self.stored_object = obj
        return EventTraversal(obj)

    def close(self):
        """"""
        del self.stored_object
