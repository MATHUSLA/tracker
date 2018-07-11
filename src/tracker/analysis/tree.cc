/*
 * src/tracker/analysis/tree.cc
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

#include <tracker/analysis/tree.hh>

#include "../helper/analysis.hh"
#include "../helper/root.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Analysis Data Tree Implementation___________________________________________________________
struct tree::impl {
  TTree* _tree;

  impl() = default;

  impl(const std::string& name,
       const std::string& title) : _tree(new TTree(name.c_str(), title.c_str())) {}

  impl(const impl& other) : _tree(other._tree->CloneTree()) {}

  impl(impl&& other) noexcept = default;

  impl& operator=(const impl& other) {
    if (this != &other) {
      _tree = other._tree->CloneTree();
    }
    return *this;
  }

  impl& operator=(impl&& other) noexcept = default;

  ~impl() = default;
};
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Default Constructor______________________________________________________
tree::tree() : _impl() {}
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Named Constructor________________________________________________________
tree::tree(const std::string& name) : tree(name, name) {}
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Named Constructor________________________________________________________
tree::tree(const std::string& name,
           const std::string& title) : _impl(std::make_unique<impl>(name, title)) {}
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Copy Constructor_________________________________________________________
tree::tree(const tree& other) : _impl(new impl(*other._impl)) {}
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Move Constructor_________________________________________________________
tree::tree(tree&& other) noexcept = default;
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Copy Assignment Operator_________________________________________________
tree& tree::operator=(const tree& other) {
  if (this != &other)
    _impl.reset(new impl(*other._impl));
  return *this;
}
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Move Assignment Operator_________________________________________________
tree& tree::operator=(tree&& other) noexcept = default;
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Destructor_______________________________________________________________
tree::~tree() = default;
//----------------------------------------------------------------------------------------------

//__Get Tree Name_______________________________________________________________________________
const std::string tree::name() const {
  return _impl->_tree->GetName();
}
//----------------------------------------------------------------------------------------------

//__Get Tree Title______________________________________________________________________________
const std::string tree::title() const {
  return _impl->_tree->GetTitle();
}
//----------------------------------------------------------------------------------------------

//__Set Tree Name_______________________________________________________________________________
void tree::name(const std::string& name) {
  _impl->_tree->SetName(name.c_str());
}
//----------------------------------------------------------------------------------------------

//__Set Tree Title______________________________________________________________________________
void tree::title(const std::string& title) {
  _impl->_tree->SetTitle(title.c_str());
}
//----------------------------------------------------------------------------------------------

//__Save Tree to File___________________________________________________________________________
bool tree::save(const std::string& path) const {
  TFile file(path.c_str(), "UPDATE");
  if (!file.IsZombie()) {
    file.cd();
    file.WriteTObject(_impl->_tree->Clone());
    file.Close();
    return true;
  }
  return false;
}
//----------------------------------------------------------------------------------------------

//__Load Tree from File_________________________________________________________________________
tree tree::load(const std::string& path,
                const std::string& name) {
  tree out;
  TFile file(path.c_str(), "READ");
  if (!file.IsZombie()) {
    TTree* test = nullptr;
    file.GetObject(name.c_str(), test);
    if (test) {
      out._impl->_tree = test;
    }
    file.Close();
  }
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */
