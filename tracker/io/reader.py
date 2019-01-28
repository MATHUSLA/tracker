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

# -------------- External Library -------------- #

from abc import ABC, abstractmethod
from collections.abc import AsyncIterator
from contextlib import contextmanager

# -------------- Tracker  Library -------------- #


class AsyncTraversal(AsyncIterator):
    """"""

    def __init__(self, path, *, collect=False, default_storage=list):
        """"""
        self.path = path
        self.collect = collect
        self.store = default_storage()
        self.started = False

    def clear(self):
        """"""
        self.store.clear()

    def reset(self):
        """"""
        self.clear()
        self.started = False

    @abstractmethod
    async def start(self):
        """"""

    @abstractmethod
    async def next_item(self):
        """"""

    async def start_loop(self):
        """"""
        while not self.started:
            self.started = await self.start()

    async def __anext__(self):
        """"""
        await self.start_loop()
        item = await self.next_item()
        if self.collect:
            self.store.append(item)
        return item

    async def collect_all(self, *, reset_first=False):
        """"""
        if reset_first:
            self.reset()
        await self.start_loop()
        while True:
            try:
                self.store.append(await self.next_item())
            except StopAsyncIteration:
                break
        return self.store


class EventTraversal(AsyncTraversal):
    """"""


class Reader(ABC):
    """"""

    @abstractmethod
    await def open(self, path, *args, **kwargs):
        """"""

    @abstractmethod
    await def close(self):
        """"""

    @asynccontextmanager
    await def load(self, path, *args, **kwargs):
        """"""
        try:
            await self.open(path, *args, **kwargs)
        finally:
            return await self.close()


class EventReader(Reader):
    """"""


class ROOTEventReader(EventReader):
    """"""
