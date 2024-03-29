#pragma once

#include "parseString.h"

#include <algorithm>
#include <any>
#include <filesystem>
#include <functional>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <typeinfo>
#include <typeindex>
#include <unordered_set>
#include <vector>

namespace clice {

inline std::string argv0; // Parser will fill this

struct ArgumentBase {
    ArgumentBase*            parent{};
    std::vector<std::string> args;
    std::string              desc;
    std::optional<std::vector<std::string>> mapping;
    std::vector<std::string>   tags;
    std::optional<std::string> completion;
    std::vector<ArgumentBase*> arguments; // child parameters
    bool                       symlink;   // a symlink for example to "slix-env" should actually call "slix env"
    std::type_index            type_index;

    std::function<void()> init;
    std::function<void(std::string_view)> fromString;
    std::function<void()> cb;
    size_t                cb_priority;

    ArgumentBase() = delete;
    ArgumentBase(ArgumentBase* parent, std::type_index idx);
    virtual ~ArgumentBase();
    ArgumentBase(ArgumentBase const&) = delete;
    ArgumentBase(ArgumentBase&&) = delete;
    auto operator=(ArgumentBase const&) -> ArgumentBase& = delete;
    auto operator=(ArgumentBase&&) -> ArgumentBase& = delete;
};

struct Register {
    std::vector<ArgumentBase*> arguments;

    static auto getInstance() -> Register& {
        static Register instance;
        return instance;
    }
};

inline ArgumentBase::ArgumentBase(ArgumentBase* parent, std::type_index idx)
    : parent{parent}
    , type_index{idx}
{
    if (parent) {
        parent->arguments.push_back(this);
    } else {
        Register::getInstance().arguments.push_back(this);
    }
}

inline ArgumentBase::~ArgumentBase() {
    if (parent) {
        auto& arguments = parent->arguments;
        arguments.erase(std::remove(arguments.begin(), arguments.end(), this), arguments.end());
    } else {
        auto& arguments = Register::getInstance().arguments;
        arguments.erase(std::remove(arguments.begin(), arguments.end(), this), arguments.end());
    }
}

struct ListOfStrings : std::vector<std::string> {
    ListOfStrings() {}
    ListOfStrings(char const* str) {
        emplace_back(str);
    }
    ListOfStrings(std::initializer_list<char const*> list) {
        for (auto l : list) {
            emplace_back(l);
        }
    }
};

template <typename T = nullptr_t, typename T2 = nullptr_t>
struct Argument {
    Argument<T2>* parent{};
    ListOfStrings args;
    bool          symlink;
    std::string desc;
    bool isSet{};
    T value{};
    mutable std::any anyType; // used if T is a callback
    std::function<void()> cb;
    size_t                cb_priority{100}; // lower priorities will be triggered before larger ones
    std::optional<std::unordered_map<std::string, T>> mapping;
    std::vector<std::string> tags;

    operator bool() const {
        return isSet;
    }

    auto operator*() const -> auto const& {
        if constexpr (std::is_invocable_v<T>) {
            using R = std::decay_t<decltype(value())>;
            if (!anyType.has_value()) {
                anyType = value();
            }
            return *std::any_cast<R>(&anyType);
        } else if constexpr (std::same_as<T, nullptr_t>) {
            []<bool type_available = false> {
                static_assert(type_available, "Can't dereference a flag");
            }();
        } else {
            return value;
        }
    }

    auto operator->() const -> auto const* {
        auto const& r = **this;
        return &r;
    }

    template <typename CB>
    auto run(CB _cb) -> nullptr_t {
        cb = _cb;
        return nullptr;
    }

    struct CTor {
        ArgumentBase arg;
        static auto detectType() -> std::type_index {
            if constexpr (std::is_invocable_v<T>) {
                // Get the type_index of a lambda
                return std::type_index(typeid(std::invoke_result_t<T>));
            }
            return std::type_index(typeid(T));
        };
        CTor(Argument& desc)
            : arg{desc.parent?&desc.parent->storage.arg:nullptr, detectType()}
        {
            arg.args    = desc.args;
            arg.symlink = desc.symlink;
            arg.desc    = desc.desc;
            if constexpr (std::same_as<std::filesystem::path, T>) {
                arg.completion = " -f ";
            } else if (desc.mapping) {
                std::string str;
                for (auto [key, value] : *desc.mapping) {
                    str += key + "\n";
                }
                arg.completion = str;
            }
            if (desc.mapping) {
                auto v = std::vector<std::string>{};
                for (auto const& [key, value] : *desc.mapping) {
                    v.push_back(key);
                }
                arg.mapping = v;
            }
            arg.tags = desc.tags;
            arg.init = [&]() {
                desc.isSet = true;
                arg.cb = desc.cb;
                arg.cb_priority = desc.cb_priority;
                if constexpr (std::same_as<nullptr_t, T>) {
                } else if constexpr (std::is_arithmetic_v<T>
                                     || std::same_as<std::string, T>
                                     || std::same_as<std::filesystem::path, T>
                                     || std::is_enum_v<T>) {
                    arg.fromString = [&](std::string_view s) {
                        if (desc.mapping) {
                            desc.value = desc.mapping->at(std::string{s});
                        } else {
                            desc.value = parseFromString<T>(s);
                        }
                        arg.fromString = nullptr;
                    };
                } else if constexpr (requires {{ std::declval<T>().push_back(std::declval<typename T::value_type>()) }; }) {
                    arg.fromString = [&](std::string_view s) {
                        if (desc.mapping) {
                            throw std::runtime_error("Type can't use mapping");
                        } else {
                            desc.value.push_back(parseFromString<typename T::value_type>(s));
                        }

                    };
                } else if constexpr (std::is_invocable_v<T>) {
                    arg.fromString = [&](std::string_view s) {
                        using RT = std::invoke_result_t<T>;
                        desc.anyType = parseFromString<RT>(s);
                    };

                } else {
                    []<bool type_available = false> {
                        static_assert(type_available, "Type can't be used as a value type in clice::Argument");
                    }();
                }
            };
        }
    } storage{*this};
};

}
