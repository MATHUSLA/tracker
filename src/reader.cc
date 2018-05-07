#include "reader.hh"

#include <array>
#include <fstream>

#include "ROOT/TFile.h"
#include "ROOT/TKey.h"
#include "ROOT/TSystemDirectory.h"
#include "ROOT/TTree.h"

#include "geometry.hh"
#include "units.hh"

namespace MATHUSLA { namespace TRACKER {

namespace reader { /////////////////////////////////////////////////////////////////////////////

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
//__ROOT File Key Traversal_____________________________________________________________________
template<class BinaryFunction>
BinaryFunction _traverse_file(const std::string& path, BinaryFunction f) {
  auto file = TFile::Open(path.c_str());
  TIter next(file->GetListOfKeys());
  TKey* key = nullptr;
  while ((key = static_cast<TKey*>(next())))
    f(file, key);
  delete key;
  return std::move(f);
}
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

namespace root { ///////////////////////////////////////////////////////////////////////////////

//__Search and Collect ROOT File Paths__________________________________________________________
std::vector<std::string> search_directory(const std::string& path) {
  std::vector<std::string> paths{};
  _collect_paths(new TSystemDirectory("data", path.c_str()), paths, "root");
  return paths;
}
//----------------------------------------------------------------------------------------------

//__Import Detector Map from File_______________________________________________________________
detector_map import_detector_map(const std::string& path) {
  std::ifstream file(path);
  detector_map out{};
  std::string line;
  while (std::getline(file, line)) {
    // TODO: how about std::isspace?
    const auto& space = line.find(" ");
    if (space == std::string::npos || line[0] == '#') {
      continue;
    } else {
      try {
        const auto& name = line.substr(0, space);
        const auto& id = std::stoll(line.substr(1 + space));
        out.insert({id, name}); // overwrites duplicate names
      } catch (...) {}
    }
  }
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
analysis::event_vector import_events(const std::string& path,
                                     const point_keys& keys) {
  analysis::event_vector out{};
  TTree* tree = nullptr;
  _traverse_file(path, [&](const auto& file, const auto& key) {
    if (std::string(key->GetClassName()) == "TTree") {
      Double_t t, x, y, z;
      tree = static_cast<TTree*>(file->Get(key->GetName()));
      tree->SetBranchAddress(keys[0].c_str(), &t);
      tree->SetBranchAddress(keys[1].c_str(), &x);
      tree->SetBranchAddress(keys[2].c_str(), &y);
      tree->SetBranchAddress(keys[3].c_str(), &z);
      analysis::event_points points{};
      const auto&& size = tree->GetEntries();
      for (auto i = 0; i < size; ++i) {
        tree->GetEntry(i);
        points.push_back({
          t * units::time,
          x * units::length,
          y * units::length,
          z * units::length});
      }
      out.push_back(points);
    }
  });
  delete tree;
  return out;
}
//----------------------------------------------------------------------------------------------

//__Import Events from ROOT File________________________________________________________________
analysis::event_vector import_events(const std::string& path,
                                     const detector_keys& keys,
                                     const detector_map& map) {
  analysis::event_vector out{};
  TTree* tree = nullptr;
  _traverse_file(path, [&](const auto& file, const auto& key) {
    if (std::string(key->GetClassName()) == "TTree") {
      Double_t t, detector;
      tree = static_cast<TTree*>(file->Get(key->GetName()));
      tree->SetBranchAddress(keys[0].c_str(), &t);
      tree->SetBranchAddress(keys[1].c_str(), &detector);
      analysis::event_points points{};
      const auto&& size = tree->GetEntries();
      for (auto i = 0; i < size; ++i) {
        tree->GetEntry(i);
        const auto& search = map.find(detector);
        const auto& name = search != map.end() ? search->second : std::to_string(detector);
        const auto& center = geometry::limits_of(name).center;
        points.push_back({t, center.x, center.y, center.z});
      }
      out.push_back(points);
    }
  });
  delete tree;
  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace root */ /////////////////////////////////////////////////////////////////////////

namespace script { /////////////////////////////////////////////////////////////////////////////

namespace { ////////////////////////////////////////////////////////////////////////////////////
//__Tracking Script Allowed Keys Iterators______________________________________________________
const auto& _allowed_keys_begin = allowed_keys.cbegin();
const auto& _allowed_keys_end = allowed_keys.cend();
//----------------------------------------------------------------------------------------------
} /* anonymous namespace */ ////////////////////////////////////////////////////////////////////

//__Tracking Script Key Check___________________________________________________________________
bool is_key_allowed(const std::string& key) {
  return _allowed_keys_end != std::find(_allowed_keys_begin, _allowed_keys_end, key);
}
//----------------------------------------------------------------------------------------------

//__Tracking Script Options Parser______________________________________________________________
const tracking_options read(const std::string& path) {
  std::ifstream file(path);
  tracking_options out{};
  std::string line;
  while (std::getline(file, line)) {
    if (line.empty()) continue;

    std::string key;
    std::string value;
    bool on_key = true;
    for (const auto& ch : line) {
      if (ch == '#') break;
      if (ch == ' ') continue;
      if (ch == ':') {
        on_key = false;
        continue;
      }
      if (on_key) key.push_back(ch);
      else value.push_back(ch);
    }

    if (!value.empty() && is_key_allowed(key))
      out[key] = value;
  }

  return out;
}
//----------------------------------------------------------------------------------------------

} /* namespace script */ ///////////////////////////////////////////////////////////////////////

} /* namespace reader */ ///////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */
