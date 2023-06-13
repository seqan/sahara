#pragma once

#include "Argument.h"

#include <fmt/format.h>

namespace clice {

inline auto generatePartialSynopsis(ArgumentBase const& arg) -> std::string {
    auto ret = fmt::format("[{}", arg.arg);
    for (auto child : arg.arguments) {
        ret += " " + generatePartialSynopsis(*child);
    }
    if (arg.type_index == std::type_index(typeid(nullptr_t))) {
        //!Nothing to do, this is a flag and doesn't take any parameters
    } else if (arg.mapping) {
        ret += fmt::format(" [{}]", fmt::join(*arg.mapping, "|"));
    } else if (arg.type_index == std::type_index(typeid(char))) {
        ret += " _char_";
    } else if (arg.type_index == std::type_index(typeid(bool))) {
        ret += " [true|false]";
    } else if (arg.type_index == std::type_index(typeid(int8_t))) {
        ret += " _int8_";
    } else if (arg.type_index == std::type_index(typeid(uint8_t))) {
        ret += " _uint8_";
    } else if (arg.type_index == std::type_index(typeid(int16_t))) {
        ret += " _int16_";
    } else if (arg.type_index == std::type_index(typeid(uint16_t))) {
        ret += " _uint16_";
    } else if (arg.type_index == std::type_index(typeid(int32_t))) {
        ret += " _int32_";
    } else if (arg.type_index == std::type_index(typeid(uint32_t))) {
        ret += " _uint32_";
    } else if (arg.type_index == std::type_index(typeid(int64_t))) {
        ret += " _int64_";
    } else if (arg.type_index == std::type_index(typeid(uint64_t))) {
        ret += " _uint64_";
    } else if (arg.type_index == std::type_index(typeid(float))) {
        ret += " _float_";
    } else if (arg.type_index == std::type_index(typeid(double))) {
        ret += " _double_";
    } else if (arg.type_index == std::type_index(typeid(std::string))) {
        ret += " _string_";
    } else if (arg.type_index == std::type_index(typeid(std::filesystem::path))) {
        ret += " _path_";
    } else if (arg.type_index == std::type_index(typeid(std::vector<bool>))) {
        ret += " [true|false] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<int8_t>))) {
        ret += " [_int8_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<uint8_t>))) {
        ret += " [_uint8_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<int16_t>))) {
        ret += " [_int16_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<uint16_t>))) {
        ret += " [_uint16_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<int32_t>))) {
        ret += " [_int32_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<uint32_t>))) {
        ret += " [_uint32_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<int64_t>))) {
        ret += " [_int64_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<uint64_t>))) {
        ret += " [_uint64_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<float>))) {
        ret += " [_float_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<double>))) {
        ret += " [_double_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<std::string>))) {
        ret += " [_string_] ...";
    } else if (arg.type_index == std::type_index(typeid(std::vector<std::filesystem::path>))) {
        ret += " [_path_] ...";

    } else {
        ret += " _unknown_";
    }
    ret += "]";
    return  ret;
}
inline auto generateSynopsis() -> std::string {
    auto ret = argv0;
    for (auto const& arg : Register::getInstance().arguments) {
        ret += " " + generatePartialSynopsis(*arg);
    }
    return ret;
}

inline auto generateSplitSynopsis() -> std::string {
    auto ret = std::string{};
    auto bases = std::vector<ArgumentBase const*>{};
    for (auto const& arg : Register::getInstance().arguments) {
        if (!arg->arg.empty() and arg->arg[0] != '-') {
            bases.emplace_back(arg);
        } else {
            ret += " " + generatePartialSynopsis(*arg);
        }
    }
    if (ret.size() > 0) {
        ret = fmt::format("{} {}\n", argv0, ret);
    }
    for (auto arg : bases) {
        auto s = generatePartialSynopsis(*arg);
        assert(s.size() >= 2);
        s.pop_back();
        s.erase(s.begin());
        ret = fmt::format("{}{} {}\n", ret, argv0, s);
    }
    return ret;
}

inline auto generateHelp() -> std::string {
    auto ret = std::string{};

    ret = generateSynopsis() + "\n\n";
    ret = generateSplitSynopsis() + "\n\n";

    auto& args = clice::Register::getInstance().arguments;
    auto f = std::function<void(std::vector<clice::ArgumentBase*>, std::string)>{};
    f = [&](auto const& args, std::string ind) {
        size_t longestWord = {};
        for (auto arg : args) {
            longestWord = std::max(longestWord, arg->arg.size());
        }

        for (auto arg : args) {
            if (arg->arg.empty() or arg->arg[0] == '-') continue;
            ret = ret + fmt::format("{}{:<{}} - {}\n", ind, arg->arg, longestWord, arg->desc);
            f(arg->arguments, ind + "    ");
        }
        for (auto arg : args) {
            if (arg->arg.empty() or arg->arg[0] != '-') continue;
            ret = ret + fmt::format("{}{:<{}}    - {}\n", ind, arg->arg, longestWord, arg->desc);
            f(arg->arguments, ind + "    ");
        }
    };
    f(args, "");
    return ret;
}
}
