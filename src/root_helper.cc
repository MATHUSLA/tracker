#include "root_helper.hh"

#include <fstream>

#include "ROOT/TSystemDirectory.h"
#include "ROOT/TTree.h"

#include "geometry.hh"
#include "units.hh"

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
  std::ifstream _file(path);
  detector_map out{};
  std::string line;
  while (std::getline(_file, line)) {
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
  traverse_file(path, [&](const auto& file, const auto& key) {
    if (std::string(key->GetClassName()) == "TTree") {
      type::real t, x, y, z;
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
          t * Units::Time,
          x * Units::Length,
          y * Units::Length,
          z * Units::Length});
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
  traverse_file(path, [&](const auto& file, const auto& key) {
    if (std::string(key->GetClassName()) == "TTree") {
      type::real t, detector;
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

} } /* namespace MATHUSLA::TRACKER */
