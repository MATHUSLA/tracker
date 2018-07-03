/*
 * src/tracker/plot/value_tag.cc
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

#include <ROOT/TFile.h>
#include <ROOT/TNamed.h>

namespace MATHUSLA { namespace TRACKER {

namespace plot { ///////////////////////////////////////////////////////////////////////////////

//__Value Tag Key Constructor___________________________________________________________________
value_tag::value_tag(const std::string& key) : _key(key), _value("") {}
//----------------------------------------------------------------------------------------------

//__Value Tag Key-Value Constructor_____________________________________________________________
value_tag::value_tag(const std::string& key,
                     const std::string& value) : _key(key), _value(value) {}
//----------------------------------------------------------------------------------------------

//__Set Key from Value Tag______________________________________________________________________
void value_tag::key(const std::string& new_key) {
  _key = new_key;
}
//----------------------------------------------------------------------------------------------

//__Get Key from Value Tag______________________________________________________________________
const std::string& value_tag::key() const {
  return _key;
}
//----------------------------------------------------------------------------------------------

//__Set Value from Value Tag____________________________________________________________________
void value_tag::value(const std::string& new_value) {
  _value = new_value;
}
//----------------------------------------------------------------------------------------------

//__Get Value from Value Tag____________________________________________________________________
const std::string& value_tag::value() const {
  return _value;
}
//----------------------------------------------------------------------------------------------

//__Save Value Tag to File______________________________________________________________________
bool value_tag::save(const std::string& path) const {
  TFile file(path.c_str(), "UPDATE");
  if (!file.IsZombie()) {
    file.cd();
    TNamed named(_key, _value);
    named.Write();
    file.Close();
    return true;
  }
  return false;
}
//----------------------------------------------------------------------------------------------

//__Load Value Tag From File____________________________________________________________________
value_tag value_tag::load(const std::string& path,
                          const std::string& key) {
  value_tag out;
  TFile file(path.c_str(), "READ");
  if (!file.IsZombie()) {
    TNamed* test = nullptr;
    file.GetObject(key.c_str(), test);
    if (test) {
      out._key = test->GetName();
      out._value = test->GetTitle();
    }
    file.Close();
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Equality of Value Tags______________________________________________________________________
bool value_tag::operator==(const value_tag& other) {
  return _key == other._key && _value == other._value;
}
//----------------------------------------------------------------------------------------------

} /* namespace plot */  ////////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */