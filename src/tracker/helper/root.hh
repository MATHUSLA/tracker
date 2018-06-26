/*
 * src/tracker/helper/root.hh
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

#ifndef TRACKER__HELPER__ROOT_HH
#define TRACKER__HELPER__ROOT_HH
#pragma once

#include <ROOT/TFile.h>
#include <ROOT/TSystemDirectory.h>
#include <ROOT/TTree.h>
#include <ROOT/TError.h>
#include <ROOT/TKey.h>

namespace MATHUSLA { namespace TRACKER {

namespace root { namespace helper { ////////////////////////////////////////////////////////////

namespace error { //////////////////////////////////////////////////////////////////////////////

//__Empty Error Handler_________________________________________________________________________
inline void empty_handler(int, Bool_t, const char*, const char*) {}
//----------------------------------------------------------------------------------------------

} /* namespace error */ ////////////////////////////////////////////////////////////////////////

//__Initialize ROOT Environment_________________________________________________________________
inline void init() {
  static bool initialized = false;
  if (!initialized) {
    gErrorIgnoreLevel = kFatal;
    SetErrorHandler(error::empty_handler);
    initialized = true;
  }
}
//----------------------------------------------------------------------------------------------

//__Recursive TSystemFile Traversal_____________________________________________________________
inline void collect_paths(TSystemDirectory* dir,
                          std::vector<std::string>& paths,
                          const std::string& ext) {
  if (!dir || !dir->GetListOfFiles()) return;
  for (const auto&& obj : *dir->GetListOfFiles()) {
    const auto file = dynamic_cast<TSystemFile*>(obj);
    const auto name = std::string(file->GetName());
    if (!file->IsDirectory()) {
      if (ext.empty() || name.substr(1 + name.find_last_of('.')) == ext)
        paths.push_back(std::string(file->GetTitle()) + "/" + name);
    } else if (!(name == "." || name == "..")) {
      collect_paths(dynamic_cast<TSystemDirectory*>(file), paths, ext);
    }
  }
}
//----------------------------------------------------------------------------------------------

//__Search and Collect ROOT File Paths__________________________________________________________
inline const std::vector<std::string> search_directory(const std::string& path,
                                                       const std::string& ext) {
  std::vector<std::string> paths{};
  collect_paths(new TSystemDirectory("search", path.c_str()), paths, ext);
  return paths;
}
//----------------------------------------------------------------------------------------------

//__ROOT File Open______________________________________________________________________________
template<class UnaryFunction>
UnaryFunction while_open(const std::string& path,
                         const std::string& mode,
                         UnaryFunction f) {
  auto file = TFile::Open(path.c_str(), mode.c_str());
  if (file && !file->IsZombie()) {
    file->cd();
    f(file);
    file->Close();
  }
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__ROOT File Traverse Keys_____________________________________________________________________
template<class BinaryFunction>
BinaryFunction traverse_keys(const std::string& path,
                             const std::string& mode,
                             BinaryFunction f) {
  while_open(path, mode, [&](const auto file) {
    TIter next(file->GetListOfKeys());
    TKey* key = nullptr;
    while ((key = dynamic_cast<TKey*>(next())))
      f(file, key);
  });
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

//__Get TKey Classname__________________________________________________________________________
inline const std::string get_key_classname(const TKey* key) {
  return key->GetClassName();
}
//----------------------------------------------------------------------------------------------

namespace tree { ///////////////////////////////////////////////////////////////////////////////

//__Safe Get TTree Object From File_____________________________________________________________
inline TTree* get_tree(TFile* file,
                       const TKey* key) {
  return (file && key && get_key_classname(key) == "TTree")
    ? dynamic_cast<TTree*>(file->Get(key->GetName()))
    : nullptr;
}
//----------------------------------------------------------------------------------------------

//__Get TTree Object From File Without Checking Type____________________________________________
inline TTree* unchecked_get_tree(TFile* file,
                                 const TKey* key) {
  return dynamic_cast<TTree*>(file->Get(key->GetName()));
}
//----------------------------------------------------------------------------------------------

//__Set TTree Branch Base Function______________________________________________________________
inline void set_branches(TTree* tree,
                         const std::string& name,
                         Double_t* value) {
  tree->SetBranchAddress(name.c_str(), value);
}
//----------------------------------------------------------------------------------------------

//__Recursively Set TTree Branches______________________________________________________________
template<class... Args>
inline void set_branches(TTree* tree,
                         const std::string& name,
                         Double_t* value,
                         Args ...args) {
  set_branches(tree, name, value);
  set_branches(tree, args...);
}
//----------------------------------------------------------------------------------------------

//__Linear Traverse Entries in TTree____________________________________________________________
template<class NullaryFunction>
NullaryFunction traverse_entries(TTree* tree,
                                 NullaryFunction f) {
  const auto entries = tree->GetEntries();
  if (entries < 0)
    return std::move(f);

  const auto size = static_cast<std::size_t>(entries);
  for (std::size_t i = 0; i < size; ++i) {
    tree->GetEntry(i);
    f();
  }
  return std::move(f);
}
//----------------------------------------------------------------------------------------------

} /* namespace tree */ /////////////////////////////////////////////////////////////////////////

} } /* namespace root::helper */ ///////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__HELPER__ROOT_HH */
