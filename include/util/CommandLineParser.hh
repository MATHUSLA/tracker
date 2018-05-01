#ifndef UTIL_COMMANDLINEPARSER_HH
#define UTIL_COMMANDLINEPARSER_HH
#pragma once

#include <initializer_list>
#include <string>
#include <type_traits>
#include <vector>

#include "util.hh"

namespace MATHUSLA {

namespace CLI {
/*
class option {
public:
  enum class argument_type { none, optional, required };
  enum class name_type { one_char, one_word, multi_word };

  option() {}

  option& short_name(const char name) {
    _short = name;
    return *this;
  }

  option& long_name(const std::string& name) {
    _long = name;
    return *this;
  }

  option& description(const std::string& desc) {
    _desc = desc;
    return *this;
  }

  std::size_t count() const { return _count; }
  const std::string& argument() const { return _argument; }
private:
  char _short;
  std::string _long;
  std::string _desc;
  std::size_t _count;
  std::string _argument;
};

template<typename ...Options>
size_t parse(char* const argv[], Options&... options) {

}
*/
} /* namespace CLI */

// TODO: needs work, unstable
struct CommandLineOption {
  enum : size_t {
    Empty             = 0x000u,
    OptionalArguments = 0x001u,
    NoArguments       = 0x002u,
    RequiredArguments = 0x004u,
    NoHyphenArguments = 0x008u,

    Repeatable             = 0x010u,
    IsShortWithoutArgument = 0x020u,
    IsShortWithArgument    = 0x040u,
    IsLongWithoutArgument  = 0x080u,
    IsLongWithArgument     = 0x100u,

    Error = 0x200u
  };

  CommandLineOption();
  CommandLineOption(char short_name,
                    const std::string& long_name="",
                    const std::string& description="",
                    size_t flags=Empty);

  inline void reset() {
    count = 0;
    argument = nullptr;
  }

  inline void set_argument(char* arg, size_t flag) {
    argument = arg;
    flags |= flag;
  }

  inline bool operator==(const CommandLineOption& rhs) {
    return short_name  == rhs.short_name
        && long_name   == rhs.long_name
        && description == rhs.description
        && flags       == rhs.flags
        && count       == rhs.count
        && argument    == rhs.argument;
  }

  char          short_name;
  std::string long_name;
  const std::string description;
  size_t        flags;
  size_t        count;
  char*         argument;
};

using CommandLineOptionList = std::vector<CommandLineOption*>;

namespace CommandLineParser {

size_t parse(char* argv[], std::initializer_list<CommandLineOption*> options);

} /* namespace CommandLineParser */

} /* namespace MATHUSLA */

#endif /* UTIL_COMMANDLINEPARSER_HH */
