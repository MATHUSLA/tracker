# -*- coding: utf-8 -*- #
#
# tracker/io/script.py
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
MATHUSLA Tracking Script.

"""

# -------------- Standard Library -------------- #

from dataclasses import dataclass, field
from typing import List, Tuple

# -------------- Tracker  Library -------------- #

from ..units import U
from ..util import classproperty


@dataclass
class ScriptOptions:
    """"""

    geometry_file: str = ''
    geometry_map_file: str = ''
    geometry_time_file: str = ''
    default_time_error: float = 0

    data_directories: List[str] = field(default_factory=list)
    data_timing_offsets: List[float] = field(default_factory=list)
    data_file_extension: str = ''
    data_t_key: str = ''
    data_x_key: str = ''
    data_y_key: str = ''
    data_z_key: str = ''
    data_dt_key: str = ''
    data_dx_key: str = ''
    data_dy_key: str = ''
    data_dz_key: str = ''
    data_detector_key: str = 'Detector'
    data_track_id_key: str = 'Track'
    data_parent_id_key: str = 'Parent'
    data_e_key: str = ''
    data_px_key: str = ''
    data_py_key: str = ''
    data_pz_key: str = ''

    statistics_directory: str = ''
    statistics_file_prefix: str = 'statistics'
    statistics_file_extension: str = 'root'
    merge_input: bool = False

    time_smearing: bool = True
    simulated_efficiency: float = 1
    simulated_noise_rate: float = 0
    event_time_window: Tuple[float, float] = field(default=(0, 0))
    layer_axis: Coordinate = Coordinate.Z
    layer_depth: float = 0
    line_width: float = 1
    seed_size: int = 3
    event_density_limit: float = 1
    event_overload_limit: float = 2
    track_density_limit: float = 1

    verbose_output: bool = False
    draw_events: bool = False

    @classproperty
    def comment_character(cls):
        """"""
        return '#'

    @classproperty
    def space_character(cls):
        """"""
        return ' '

    @classproperty
    def key_value_separator(cls):
        """"""
        return ':'

    @classproperty
    def continuation_string(cls):
        """"""
        return '...'

    @classproperty
    def continuation_line_character(cls):
        """"""
        return '-'
