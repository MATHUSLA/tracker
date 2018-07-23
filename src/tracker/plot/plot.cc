/*
 * src/tracker/plot/plot.cc
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

#include <tracker/plot.hh>

#include <ROOT/TApplication.h>

#include "../helper/root.hh"

namespace MATHUSLA { namespace TRACKER {

namespace plot { ///////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Plotting Application________________________________________________________________________
TApplication* _app = nullptr;
thread_local bool _app_on = false;
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Start Plotting Environment__________________________________________________________________
void init(bool on) {
  root::helper::init(false);
  _app_on = on;
  if (!_app && _app_on) {
    static int argc = 1;
    static std::array<char*, 1> argv{strdup("app")};
    _app = new TApplication("app", &argc, argv.data());
  }
}
//----------------------------------------------------------------------------------------------

//__End Plotting Environment____________________________________________________________________
void end() {
  if (_app) {
    try {
      root::helper::set_batch_mode(false);
      _app->Run(true);
      delete _app;
      _app = nullptr;
    } catch (...) {}
  }
}
//----------------------------------------------------------------------------------------------

//__Check if Plotting is On_____________________________________________________________________
bool is_on() {
  return _app_on;
}
//----------------------------------------------------------------------------------------------

//__Default Colors______________________________________________________________________________
const color color::WHITE   = {255, 255, 255};
const color color::BLACK   = {  0,   0,   0};
const color color::RED     = {255,   0,   0};
const color color::GREEN   = {  0, 255,   0};
const color color::BLUE    = {  0,   0, 255};
const color color::CYAN    = {  0, 255, 255};
const color color::MAGENTA = {255,   0, 255};
const color color::YELLOW  = {255, 255,   0};
//----------------------------------------------------------------------------------------------

} /* namespace plot */  ////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */
