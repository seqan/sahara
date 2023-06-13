#pragma once

#include <set>
#include <stdexcept>
#include <string>
#include <vector>

struct Config {
    std::string generator = "h2-k2";
    bool generator_dyn = false;
    std::filesystem::path saveOutput;
    size_t minK{0}, maxK{6}, k_stepSize{1};
    bool reverse{true};

    std::vector<std::string> algorithms;

    std::string queryPath{};
    std::string indexPath{};

    enum class Mode {
        All,
        BestHits,
    };
    Mode mode {Mode::All};
    size_t maxHitsPerQuery{0};
};
