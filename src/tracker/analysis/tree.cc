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

#include <ROOT/TFile.h>
#include <ROOT/TTree.h>

// TODO: improve this dependency
#include "../helper/root.hh"

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Analysis Data Tree Implementation___________________________________________________________
struct tree::impl {
  struct TTreeType : public TTree {
    TTreeType() : TTree() { root::helper::init(); }
    TTreeType(const char* name,
              const char* title) : TTree(name, title) { root::helper::init(); }
    using TTree::BranchImp;
    using TTree::BranchImpRef;
  };

  TTreeType* _tree;

  impl() = default;
  impl(const std::string& name,
       const std::string& title) : _tree(new TTreeType(name.c_str(), title.c_str())) {}
  impl(const impl& other) : _tree(dynamic_cast<TTreeType*>(other._tree->CloneTree())) {}
  impl(impl&& other) noexcept = default;
  ~impl() = default;

  impl& operator=(const impl& other) {
    if (this != &other)
      _tree = dynamic_cast<TTreeType*>(other._tree->CloneTree());
    return *this;
  }
  impl& operator=(impl&& other) noexcept = default;
};
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Default Constructor______________________________________________________
tree::tree() : _impl(std::make_unique<impl>()) {}
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

//__Draw Tree Branches__________________________________________________________________________
void tree::draw(const key_type& key) const {
  _impl->_tree->Draw(key.c_str());
}
//----------------------------------------------------------------------------------------------

//__Draw Tree Branches__________________________________________________________________________
void tree::draw(std::initializer_list<key_type> keys) const {
  std::string combined_key;
  auto begin = keys.begin();
  const auto end = keys.end();
  if (end == begin)
    return;

  combined_key.reserve(64UL * (end - begin));
  while (begin != end - 1)
    combined_key.append(*begin++ + ":");
  combined_key.append(*begin);

  draw(combined_key);
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
    auto test = static_cast<TTree*>(file.Get(name.c_str()));
    if (test) {
      test->SetDirectory(nullptr);
      out._impl->_tree = static_cast<impl::TTreeType*>(test);
    }
    file.Close();
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Add Branch to Tree__________________________________________________________________________
void* tree::create_memory_slot(const key_type& key,
                               void* memory,
                               const std::type_info& info) {
  const auto classtype = TBuffer::GetClass(info);
  const auto datatype = TDataType::GetType(info);
  _impl->_tree->BranchImpRef(key.c_str(), classtype, datatype, memory, 32000, 99);
  return memory;
}
//----------------------------------------------------------------------------------------------

//__Get Branch Memory from Tree_________________________________________________________________
void* tree::load_memory_slot(const key_type& key,
                             void* memory) {
  _impl->_tree->SetBranchAddress(key.c_str(), memory);
  return memory;
}
//----------------------------------------------------------------------------------------------

//__Get Entry from Tree_________________________________________________________________________
void tree::operator[](const std::size_t index) const {
  _impl->_tree->GetEntry(index);
}
//----------------------------------------------------------------------------------------------

//__Get Size of Tree____________________________________________________________________________
std::size_t tree::size() const {
  const auto out = _impl->_tree->GetEntries();
  return out <= 0 ? 0UL : static_cast<std::size_t>(out);
}
//----------------------------------------------------------------------------------------------

//__Get Count of Key____________________________________________________________________________
std::size_t tree::count(const key_type& key) const {
  // TODO: implement
  return 0;
}
//----------------------------------------------------------------------------------------------

//__Fill Tree___________________________________________________________________________________
void tree::fill() {
  _impl->_tree->Fill();
}
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */
