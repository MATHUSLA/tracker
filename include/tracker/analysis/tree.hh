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
#include <typeinfo>

#include <tracker/analysis/type.hh>

#include <iostream>

namespace MATHUSLA { namespace TRACKER {

namespace analysis { ///////////////////////////////////////////////////////////////////////////

//__Analysis Data Tree__________________________________________________________________________
class tree {
public:
  using key_type = std::string;

  template<class T>
  class branch;

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

  std::size_t size() const;
  std::size_t count(const key_type& key) const;

  void operator[](const std::size_t index) const;
  void fill();

  template<class T, class ...Args>
  branch<T> add_branch(const key_type& key,
                       Args&& ...args) {
    return branch<T>(*this, key, std::forward<Args>(args)...);
  }

  template<class T, class ...Args>
  branch<T> get_branch(const key_type& key,
                       Args&& ...args) {
    return branch<T>::load(*this, key, std::forward<Args>(args)...);
  }

  template<class NullaryFunction>
  NullaryFunction traverse(NullaryFunction f) {
    const auto s = size();
    for (std::size_t i{}; i < s; ++i) {
      operator[](i);
      f();
    }
    return std::move(f);
  }
  template<class NullaryFunction>
  NullaryFunction traverse(NullaryFunction f) const { return traverse(f); }

  template<class UnaryFunction>
  UnaryFunction traverse_with_index(UnaryFunction f) {
    const auto s = size();
    for (std::size_t i{}; i < s; ++i) {
      operator[](i);
      f(i);
    }
    return std::move(f);
  }
  template<class UnaryFunction>
  UnaryFunction traverse_with_index(UnaryFunction f) const { return traverse_with_index(f); }

  void draw(const key_type& key) const;
  void draw(std::initializer_list<key_type> keys) const;

  bool save(const std::string& path) const;

  static tree load(const std::string& path,
                   const std::string& name);

protected:
  void* create_memory_slot(const key_type& key,
                           void* memory,
                           const std::type_info& info);
  void* load_memory_slot(const key_type& key,
                         void* memory);

  struct impl;
  std::unique_ptr<impl> _impl;
};
//----------------------------------------------------------------------------------------------

//__Analysis Data Tree Branch___________________________________________________________________
template<class T>
class tree::branch {
public:
  using value_type = T;

  template<class ...Args>
  branch(tree& base,
         const tree::key_type& key,
         Args&& ...args)
      : branch(key, base) {
    create_memory_slot(new T(std::forward<Args>(args)...));
  }

  template<class ...Args>
  static branch load(tree& base,
                     const tree::key_type& key,
                     Args&& ...args) {
    branch out(key, base);
    out.load_memory_slot(new T(std::forward<Args>(args)...));
    return out;
  }

  branch(branch&& other) noexcept = default;
  branch& operator=(branch&& other) noexcept = default;
  ~branch() = default;

  T& get() { return *_ptr; }
  const T& get() const { return *_ptr; }
  void set(const T& t) { *_ptr = t; }

  operator const T&() const { return *_ptr; }

  tree& base() { return _base; }
  const tree& base() const { return _base; }

protected:
  branch() = default;
  branch(const tree::key_type& key,
         tree& base)
      : _base(&base), _key(key) {}

  T* create_memory_slot(T* memory) {
    return _ptr = static_cast<T*>(_base->create_memory_slot(_key, memory, typeid(T)));
  }

  T* load_memory_slot(T* memory) {
    return _ptr = static_cast<T*>(_base->load_memory_slot(_key, memory));
  }

  T* _ptr;
  tree* _base;
  tree::key_type _key;
};
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__TREE_HH */
