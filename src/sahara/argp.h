#pragma once

#include <set>
#include <stdexcept>
#include <string>
#include <vector>

struct Config {
    std::string generator = "h2-k2";
    bool generator_dyn = false;
    std::filesystem::path saveOutput;
    size_t k;
    bool reverse{true};

    std::string queryPath{};
    std::string indexPath{};

    enum class Mode {
        All,
        BestHits,
    } mode;
    size_t maxHitsPerQuery{0};
};
