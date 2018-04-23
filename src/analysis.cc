#include "analysis.hh"

#include <ROOT/TFile.h>
#include <ROOT/TSystemDirectory.h>
#include <ROOT/TTree.h>
#include <ROOT/TKey.h>

namespace MATHUSLA { namespace TRACKER {

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Recursive TSystemFile Traversal_____________________________________________________________
void _collect_paths(TSystemDirectory* dir,
                    std::vector<std::string>& paths,
                    const std::string& ext) {
  if (!dir || !dir->GetListOfFiles()) return;
  for (const auto&& obj : *dir->GetListOfFiles()) {
    const auto file = static_cast<TSystemFile*>(obj);
    const auto&& name = std::string(file->GetName());
    if (!file->IsDirectory()) {
      if (ext == "" || name.substr(1 + name.find_last_of(".")) == ext)
        paths.push_back(std::string(file->GetTitle()) + "/" + name);
    } else if (!(name == "." || name == "..")) {
      _collect_paths(static_cast<TSystemDirectory*>(file), paths, ext);
    }
  }
}
//----------------------------------------------------------------------------------------------
} /* annonymous namespace */ ///////////////////////////////////////////////////////////////////

//__Search and Collect ROOT File Paths__________________________________________________________
std::vector<std::string> Analysis::search_directory(const std::string& path) {
  std::vector<std::string> paths{};
  _collect_paths(new TSystemDirectory("data", path.c_str()), paths, "root");
  return paths;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
Analysis::event_vector Analysis::import_events(const std::string& path,
                                               const std::array<std::string, 4>& keys) {
  auto out = event_vector{};
  auto file = TFile::Open(path.c_str());
  TIter next(file->GetListOfKeys());
  TKey* key = nullptr;
  TTree* tree = nullptr;
  while ((key = static_cast<TKey*>(next()))) {
    if (std::string(key->GetClassName()) == "TTree") {
      real t, x, y, z;
      tree = static_cast<TTree*>(file->Get(key->GetName()));
      tree->SetBranchAddress(keys[0].c_str(), &t);
      tree->SetBranchAddress(keys[1].c_str(), &x);
      tree->SetBranchAddress(keys[2].c_str(), &y);
      tree->SetBranchAddress(keys[3].c_str(), &z);
      auto points = r4_point_vector{};
      const auto&& size = tree->GetEntries();
      for (auto i = 0; i < size; ++i) {
        tree->GetEntry(i);
        points.push_back({t, x, y, z});
      }
      out.push_back(points);
    }
  }
  delete key;
  delete tree;
  return out;
}
//----------------------------------------------------------------------------------------------

} } /* namespace MATHUSLA::TRACKER */
