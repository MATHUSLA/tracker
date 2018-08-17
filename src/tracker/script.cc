/*
 * src/tracker/script.cc
 *
 * Copyright 2018 Brandon Gomes
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <tracker/script.hh>

namespace MATHUSLA { namespace TRACKER {

namespace script { /////////////////////////////////////////////////////////////////////////////

//__Parse Key Value Pair into File Path_________________________________________________________
void parse_file_path(const std::string& key,
                     const std::string& value,
                     path_type& out,
                     bool exit_on_error) {
  util::error::exit_when(exit_on_error && !util::io::path_exists(value),
    "[FATAL ERROR] Missing File Path \"", value, "\" for Key \"", key, "\" in Tracking Script.\n");
  out = value;
}
//----------------------------------------------------------------------------------------------

//__Parse Data Key Value Pair___________________________________________________________________
void parse_data_keys(const std::string& key,
                     const std::string& value,
                     std::string& t_key,
                     std::string& x_key,
                     std::string& y_key,
                     std::string& z_key,
                     bool exit_on_error) {
  std::vector<std::string> data_keys;
  util::string::split(value, data_keys, ",;");
  const auto size = data_keys.size();
  const auto& end = data_keys.cend();
  util::error::exit_when(exit_on_error && std::find(data_keys.cbegin(), end, "") != end,
    "[FATAL ERROR] Invalid Data Argument for \"", key, "\" in Tracking Script.\n"
    "              Expected 1 argument \"T\" or 4 arguments \"T, X, Y, Z\".\n");
  if (size == 1) {
    t_key = data_keys[0];
  } else if (size == 4) {
    t_key = data_keys[0];
    x_key = data_keys[1];
    y_key = data_keys[2];
    z_key = data_keys[3];
  } else if (exit_on_error) {
    util::error::exit(
      "[FATAL ERROR] Invalid Data Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected 1 argument \"T\" or 4 arguments \"T, X, Y, Z\" "
                    "but received ", size, ".\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Boolean Key Value_____________________________________________________________________
void parse_boolean(const std::string& key,
                   const std::string& value,
                   bool& out,
                   bool exit_on_error) {
  const auto value_lowercase = util::string::tolower(value);
  if (value_lowercase == "true" || value_lowercase == "t" || value_lowercase == "1") {
    out = true;
  } else if (value_lowercase == "false" || value_lowercase == "f" || value_lowercase == "0") {
    out = false;
  } else if (exit_on_error) {
    util::error::exit(
      "[FATAL ERROR] Invalid Boolean Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible boolean value.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Real Key Value________________________________________________________________________
void parse_real(const std::string& key,
                const std::string& value,
                real& out,
                bool exit_on_error) {
  try {
    out = static_cast<real>(std::stold(value));
  } catch (...) {
    util::error::exit_when(exit_on_error,
      "[FATAL ERROR] Invalid Real Number Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to floating point value.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Positive Real Key Value_______________________________________________________________
void parse_positive_real(const std::string& key,
                         const std::string& value,
                         real& out,
                         bool exit_on_error) {
  try {
    out = static_cast<real>(std::stold(value));
    if (out < 0) {
      out = 0;
      throw std::invalid_argument{"Real Number Must be Greater than or Equal to Zero."};
    }
  } catch (...) {
    util::error::exit_when(exit_on_error,
      "[FATAL ERROR] Invalid Real Number Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to positive floating point value.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Real Range Key Value__________________________________________________________________
void parse_real_range(const std::string& key,
                      const std::string& value,
                      real_range& out,
                      bool exit_on_error) {
  try {
    std::vector<std::string> values;
    util::string::split(value, values, ",;");
    if (values.size() != 2)
      throw std::invalid_argument{"Size Mismatch."};

    out = real_range{static_cast<real>(std::stold(values[0])),
                     static_cast<real>(std::stold(values[1]))};

  } catch (...) {
    util::error::exit_when(exit_on_error,
      "[FATAL ERROR] Invalid Real Number Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to Real range (min, max).\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Integer Key Value_____________________________________________________________________
void parse_integer(const std::string& key,
                   const std::string& value,
                   integer& out,
                   bool exit_on_error) {
  try {
    out = static_cast<integer>(std::stoll(value));
  } catch (...) {
    util::error::exit_when(exit_on_error,
      "[FATAL ERROR] Invalid Integer Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to integral value.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Size Type Key Value___________________________________________________________________
void parse_size_type(const std::string& key,
                     const std::string& value,
                     std::size_t& out,
                     bool exit_on_error) {
  try {
    out = static_cast<std::size_t>(std::stoull(value));
  } catch (...) {
    util::error::exit_when(exit_on_error,
      "[FATAL ERROR] Invalid Size Type Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected argument convertible to positive integral value.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse R4 Key Value__________________________________________________________________________
void parse_r4_point(const std::string& key,
                    const std::string& value,
                    r4_point& out,
                    bool use_units,
                    bool exit_on_error) {
  std::vector<std::string> point;
  util::string::split(value, point, ",;");
  const auto size = point.size();
  const auto& end = point.cend();
  util::error::exit_when(exit_on_error && (size != 4 || (std::find(point.cbegin(), end, "") != end)),
    "[FATAL ERROR] Invalid R4 Point Argument for \"", key, "\" in Tracking Script.\n"
    "              Expected 4 arguments \"T, X, Y, Z\" but received ", size, ".\n");
  try {
    out = r4_point{
      (use_units ? units::time   : 1.0L) * static_cast<real>(std::stold(point[0])),
      (use_units ? units::length : 1.0L) * static_cast<real>(std::stold(point[1])),
      (use_units ? units::length : 1.0L) * static_cast<real>(std::stold(point[2])),
      (use_units ? units::length : 1.0L) * static_cast<real>(std::stold(point[3])) };
  } catch (...) {
    util::error::exit_when(exit_on_error,
      "[FATAL ERROR] Invalid R4 Point Argument for \"", key, "\" in Tracking Script.\n"
      "              Expected 4 arguments \"T, X, Y, Z\" convertible to floating point values.\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse R3 Coordinate Key Value_______________________________________________________________
void parse_r3_coordinate(const std::string& key,
                         const std::string& value,
                         Coordinate& coordinate,
                         bool exit_on_error) {
  if (value == "x" || value == "X") {
    coordinate = Coordinate::X;
  } else if (value == "y" || value == "Y") {
    coordinate = Coordinate::Y;
  } else if (value == "z" || value == "Z") {
    coordinate = Coordinate::Z;
  } else {
    util::error::exit_when(exit_on_error,
      "[FATAL ERROR] Invalid R3 Coordinate Argument for \"", key, "\" in Tracking Script.\n",
      "              Expected X/x, Y/y, or Z/z but received \"", value, "\".\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse R4 Coordinate Key Value_______________________________________________________________
void parse_r4_coordinate(const std::string& key,
                         const std::string& value,
                         Coordinate& coordinate,
                         bool exit_on_error) {
  if (value == "t" || value == "T") {
    coordinate = Coordinate::T;
  } else if (value == "x" || value == "X") {
    coordinate = Coordinate::X;
  } else if (value == "y" || value == "Y") {
    coordinate = Coordinate::Y;
  } else if (value == "z" || value == "Z") {
    coordinate = Coordinate::Z;
  } else {
    util::error::exit_when(exit_on_error,
      "[FATAL ERROR] Invalid R4 Coordinate Argument for \"", key, "\" in Tracking Script.\n",
      "              Expected T/t, X/x, Y/y, or Z/z but received \"", value, "\".\n");
  }
}
//----------------------------------------------------------------------------------------------

//__Parse Line from Tracking Script_____________________________________________________________
void parse_line(const std::string& line,
                std::string& key,
                std::string& value) {
  bool on_key = true;
  for (const auto& ch : line) {
    if (ch == reserved::comment_character) break;
    if (ch == reserved::space_character) continue;
    if (ch == reserved::key_value_separator) {
      on_key = false;
      continue;
    }
    if (on_key) key.push_back(ch);
    else value.push_back(ch);
  }
}
//----------------------------------------------------------------------------------------------

//__Check that a Line has a Continuation String_________________________________________________
bool is_continuation_header(const std::string& key,
                            const std::string& value) {
  return value == reserved::continuation_string;
}
//----------------------------------------------------------------------------------------------

//__Check that a Line has a Continuation String_________________________________________________
bool is_continuation_body(const std::string& key,
                          const std::string& value) {
  return key.front() == reserved::continuation_line_character;
}
//----------------------------------------------------------------------------------------------

//__Parse Line from a Continuation______________________________________________________________
void parse_continuation_line(const std::string& line,
                             std::string& entry) {
  bool on_entry = false;
  for (const auto& ch : line) {
    if (ch == reserved::comment_character) break;
    if (ch == reserved::space_character) continue;
    if (!on_entry && ch == reserved::continuation_line_character) {
      on_entry = true;
      continue;
    }
    if (on_entry) entry.push_back(ch);
  }
}
//----------------------------------------------------------------------------------------------

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */
