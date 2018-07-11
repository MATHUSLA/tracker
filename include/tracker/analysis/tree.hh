/*
 * include/tracker/analysis/tree.hh
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

#ifndef TRACKER__ANALYSIS__TREE_HH
#define TRACKER__ANALYSIS__TREE_HH
#pragma once

#include <memory>

#include <tracker/analysis/type.hh>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Analysis Data Tree__________________________________________________________________________
class tree {
public:
  tree();
  tree(const std::string& name);
  tree(const std::string& name,
       const std::string& title);
  tree(const tree& other);
  tree(tree&& other) noexcept;
  ~tree();

  tree& operator=(const tree& other);
  tree& operator=(tree&& other) noexcept;

  const std::string name() const;
  const std::string title() const;
  void name(const std::string& name);
  void title(const std::string& title);

  bool save(const std::string& path) const;

  static tree load(const std::string& path,
                   const std::string& name);

private:
  struct impl;
  std::unique_ptr<impl> _impl;
};
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__TREE_HH */
