#pragma once

#include "Argument.h"
#include "printCompletion.h"

#include <cassert>
#include <fmt/format.h>
#include <list>
#include <map>
#include <span>

namespace clice {

namespace {
inline void makeCompletionSuggestion(std::vector<ArgumentBase*> const& activeBases, std::string_view arg) {
    // single completion
    if (activeBases.size() and activeBases.back()->fromString and activeBases.back()->completion and (arg.empty() || arg[0] != '-')) {
        fmt::print("{}", *activeBases.back()->completion);
        return;
    }
    auto options = std::vector<std::tuple<std::string, std::string>>{};
    for (auto bases : activeBases) {
        for (auto arg : bases->arguments) {
            if (!arg->args.empty()) {
                options.emplace_back(arg->args[0], arg->desc);
            }
        }
    }
    for (auto arg : Register::getInstance().arguments) {
        if (!arg->args.empty()) {
            options.emplace_back(arg->args[0], arg->desc);
        }
    }
    size_t longestWord = {};
    for (auto const& [arg, desc] : options) {
        longestWord = std::max(longestWord, arg.size());
    }
    //!TODO maybe we also want to show descriptions?

    std::ranges::sort(options);
    for (auto const& [arg, desc] : options) {
        fmt::print("{}\n", arg);
    }
}
}

auto parse(int argc, char const* const* argv, bool allowSingleDash = false) -> std::optional<std::string>;

inline auto parseSingleDash(int _argc, char const* const* _argv) -> std::optional<std::string> {
    auto args   = std::list<std::string>{};
    auto argptr = std::vector<char const*>{};
    bool allTrailing{false};
    for (int i{0}; i < _argc; ++i) {
        auto view = std::string_view{_argv[i]};
        if (allTrailing
            || view.size() <= 2
            || view[0] != '-'
            || view[1] == '-'
        ) {
            argptr.push_back(_argv[i]);
            if (view == "--") {
                allTrailing = true;
            }
        } else {
            for (size_t j{1}; j < view.size(); ++j) {
                args.push_back(std::string{"-"} + view[j]);
                argptr.push_back(args.back().c_str());
            }
        }
    }
    return parse(argptr.size(), argptr.data(), false);
}


/**
 * allowSingleDash: allows flags like "-a -b" be combined to "-ab"
 */
inline auto parse(int argc, char const* const* argv, bool allowSingleDash) -> std::optional<std::string> {
    if (allowSingleDash) {
        return parseSingleDash(argc, argv);
    }

    assert(argc > 0);
    clice::argv0 = argv[0];

    // check for symlink (only root arguments are considered)
    auto redirectedArguments = std::vector<char const*>{};
    for (auto arg : Register::getInstance().arguments) {
        if (arg->symlink and std::filesystem::path{argv0}.filename().string().ends_with("-" + arg->args[0])) {
            redirectedArguments.push_back(argv[0]);
            redirectedArguments.push_back(arg->args[0].c_str());
            for (auto i{1}; i < argc; ++i) {
                redirectedArguments.push_back(argv[i]);
            }
            redirectedArguments.push_back(nullptr);
            argc = redirectedArguments.size()-1;
            argv = redirectedArguments.data();
        }
    }

    if (auto gen = std::getenv("CLICE_GENERATE_COMPLETION"); gen != nullptr) {
        printCompletion(gen);
        exit(0);
    }

    auto activeBases = std::vector<ArgumentBase*>{}; // current commands whos sub arguments are of interest;

    auto completion = std::getenv("CLICE_COMPLETION") != nullptr;

    auto findRootArg = [&](std::string_view str) -> ArgumentBase* {
        for (auto arg : Register::getInstance().arguments) {
            auto iter = std::find(arg->args.begin(), arg->args.end(), str);
            if (iter != arg->args.end() and arg->init) {
                return arg;
            }
        }
        return nullptr;
    };

    auto findActiveArg = [&](std::string_view str, ArgumentBase* base) -> ArgumentBase* {
        for (auto arg : base->arguments) {
            auto iter = std::find(arg->args.begin(), arg->args.end(), str);
            if (iter != arg->args.end() and arg->init) {
                return arg;
            }
        }
        return nullptr;
    };


    bool allTrailing = false;
    for (int i{1}; i < argc; ++i) {
        // make suggestion about next possible tokens
        if (completion and i+1 == argc) {
            makeCompletionSuggestion(activeBases, std::string_view{argv[i]});
            continue;
        }
        if (std::string_view{argv[i]} == "--" and !allTrailing) {
            allTrailing = true;
            continue;
        }

        [&]() {
            // walk up the arguments, until one active argument has a child with fitting parameter
            for (size_t j{0}; j < activeBases.size(); ++j) {
                auto const& base = activeBases[activeBases.size()-j-1];
                if ((argv[i][0] != '-' or allTrailing) and base->fromString) {
                    base->fromString(argv[i]);
                    return;
                }
                if (auto arg = findActiveArg(argv[i], base); arg) {
                    arg->init();
                    activeBases.push_back(arg);
                    return;
                }
            }
            auto arg = findRootArg(argv[i]);
            if (arg) {
                arg->init();
                activeBases.push_back(arg);
                return;
            }
            throw std::runtime_error{std::string{"unexpected cli argument \""} + argv[i] + "\""};
        }();
    }
    if (completion) {
        exit(0);
    }

    // create list of all triggers according to priority
    auto triggers = std::map<size_t, std::vector<std::function<void()>>>{};
    auto f = std::function<void(std::vector<clice::ArgumentBase*>)>{};
    f = [&](auto const& args) {
        for (auto arg : args) {
            if (arg->cb) {
                triggers[arg->cb_priority].emplace_back([=]() {
                    arg->cb();
                });
            }
            f(arg->arguments);
        }
    };
    f(clice::Register::getInstance().arguments);
    // call triggers in priority level order
    for (auto const& [level, cbs] : triggers) {
        for (auto const& cb : cbs) {
            cb();
        }
    }

    return std::nullopt;
}

}
