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

  template<class T>
  branch<T> insert_branch(const key_type& key,
                          T* value) {
    return branch<T>(*this, key, value);
  }

  template<class T, class ...Args>
  branch<T> emplace_branch(const key_type& key,
                           Args&& ...args) {
    return branch<T>::dynamic(*this, key, std::forward<Args>(args)...);
  }

  template<class T>
  branch<T> get_branch(const key_type& key,
                       T* value) {
    return branch<T>::load(*this, key, value);
  }

  template<class T, class ...Args>
  branch<T> get_dynamic_branch(const key_type& key,
                               Args&& ...args) {
    return branch<T>::load_dynamic(*this, key, std::forward<Args>(args)...);
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

  branch(tree& base,
         const tree::key_type& key,
         T* value)
      : branch(key, base, value) {
    create_memory_slot(value);
  }

  branch(const branch&) = delete;
  branch(branch&&) noexcept = default;
  branch& operator=(const branch&) = delete;
  branch& operator=(branch&&) noexcept = default;
  ~branch() = default;

  template<class ...Args>
  static branch dynamic(tree& base,
                        const tree::key_type& key,
                        Args&& ...args) {
    branch out(key, base, new T(std::forward<Args>(args)...));
    out.create_memory_slot(out._ptr.get());
    out._ptr.get_deleter().delete_switch = true;
    return out;
  }

  static branch load(tree& base,
                     const tree::key_type& key,
                     T* value) {
    branch out(key, base, value);
    out.load_memory_slot(value);
    return out;
  }

  template<class ...Args>
  static branch load_dynamic(tree& base,
                             const tree::key_type& key,
                             Args&& ...args) {
    branch out(key, base, new T(std::forward<Args>(args)...));
    out.load_memory_slot(out._ptr.get());
    out._ptr.get_deleter().delete_switch = true;
    return out;
  }

  T& get() { return *_ptr; }
  const T& get() const { return *_ptr; }
  void set(const T& t) { *_ptr = t; }
  void move(T&& t) { *_ptr = std::move(t); }

  operator T() const { return *_ptr; }
  operator const T&() const { return *_ptr; }
  operator T&() { return *_ptr; }
  branch& operator=(const T& value) {
    *_ptr = value;
    return *this;
  }
  branch& operator=(T&& value) {
    *_ptr = value;
    return *this;
  }

  tree& base() { return *_base; }
  const tree& base() const { return *_base; }

protected:
  branch() = default;

  template<class U>
  struct switch_delete {
    switch_delete() = default;
    switch_delete(bool b) : delete_switch(b) {}
    void operator()(U* u) { if (delete_switch) delete u; }
    bool delete_switch = false;
  };

  branch(const tree::key_type& key,
         tree& base,
         T* ptr)
      : _ptr(ptr, switch_delete<T>{}), _base(&base), _key(key) {}

  T* create_memory_slot(T* memory) {
    return static_cast<T*>(_base->create_memory_slot(_key, memory, typeid(T)));
  }

  T* load_memory_slot(T* memory) {
    return static_cast<T*>(_base->load_memory_slot(_key, memory));
  }

  std::unique_ptr<T, switch_delete<T>> _ptr;
  tree* _base;
  tree::key_type _key;
};
//----------------------------------------------------------------------------------------------

} /* namespace analysis */ /////////////////////////////////////////////////////////////////////

} } /* namespace MATHUSLA::TRACKER */

#endif /* TRACKER__ANALYSIS__TREE_HH */
