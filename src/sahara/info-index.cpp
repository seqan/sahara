// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "utils/error_fmt.h"
#include "VarIndex.h"

#include <channel/channel.h>
#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <clice/clice.h>
#include <fmindex-collection/suffixarray/DenseCSA.h>
#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/search/SearchNg27.h>
#include <fmindex-collection/search/SearchNg28.h>
#include <fmindex-collection/search/SearchNg28KStep.h>
#include <fmindex-collection/search/SearchNg28Options.h>
#include <fstream>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <mmser/mmser.h>
#include <string>
#include <unordered_set>


namespace {
void app();
auto cli = clice::Argument {
    .args   = "info-index",
    .desc   = "print basic information about the index",
    .cb     = app,
};

auto cliIndex = clice::Argument {
    .parent = &cli,
    .args   = {"-i", "--index"},
    .desc   = "path to the index file",
    .value  = std::filesystem::path{},
};

template <typename Alphabet>
void runInfo() {
    constexpr size_t Sigma = Alphabet::size();

    // load Index
    if (!std::filesystem::exists(*cliIndex)) {
        throw error_fmt{"no valid index path at {}", *cliIndex};
    }

    using Index = VarIndex<Alphabet>;
    auto [varIndex, storageManager] = [&]() -> std::tuple<Index, std::unique_ptr<std::any>> {
        if (cliIndex->string().ends_with(".mmser")) {
            return mmser::loadFileStream<Index>(*cliIndex);
        } else {
            auto varIndex = VarIndex<Alphabet>{};
            auto ifs     = std::ifstream{*cliIndex, std::ios::binary};
            auto archive = cereal::BinaryInputArchive{ifs};
            archive(varIndex);
            return {std::move(varIndex), std::unique_ptr<std::any>{}};
        }
    }();
    fmt::print("  samplingRate: {}\n", varIndex.samplingRate);



    auto indexSize = std::visit([](auto const& index) {
        return index.size();
    }, varIndex.vs);

    std::visit([&]<typename Index>(Index const& index) {
        auto& bwt = index.bwt;

        std::vector<uint64_t> monotonicBlocks{};
        bool currentValue = {};
        size_t number{};

        monotonicBlocks.resize(10);
        for (size_t i{0}; i < indexSize; ++i) {
            auto s = bwt.symbol(i);
            if (currentValue == s) {
                number += 1;
            } else {
                number = 1;
                currentValue = s;
            }


            for (size_t bs{1}; bs < monotonicBlocks.size(); ++bs) {
                if ((i+1) % bs == 0) {
                    if (number > bs) {
                        monotonicBlocks[bs] += 1;
                    }
                }
            }
        }

        for (size_t i{1}; i < monotonicBlocks.size(); ++i) {
            fmt::print("blocks of size {} ar in {} ({}%) of {} blocks monotonic\n", i, monotonicBlocks[i], monotonicBlocks[i] *100 / (indexSize / i), indexSize / i);
        }

    }, varIndex.vs);
}

void app() {
    auto mmser_loading = cliIndex->string().ends_with(".mmser");

    // load sigma value
    size_t sigma;
    std::string indexType;
    auto path = cliIndex->string();
    if (mmser_loading) {
        auto archive = mmser::ArchiveLoadStream{path};
        mmser::handle(archive, sigma);
        size_t samplingRate;
        mmser::handle(archive, samplingRate);
        mmser::handle(archive, indexType);
    } else {
        auto ifs     = std::ifstream{path, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        archive(sigma);
        size_t samplingRate;
        archive(samplingRate);
        archive(indexType);
    }
    auto nd = [](std::string str) {
        return str.ends_with("-nd") || str.ends_with("-nd-rev");
    };
    if (sigma == 2 && nd(indexType)) runInfo<ivs::dna2>();
    else if (sigma == 3 && !nd(indexType)) runInfo<ivs::d_dna2>();
    else if (sigma == 4 &&  nd(indexType)) runInfo<ivs::dna4>();
    else if (sigma == 5 && !nd(indexType)) runInfo<ivs::d_dna4>();
    else if (sigma == 5 &&  nd(indexType)) runInfo<ivs::dna5>();
    else if (sigma == 6 && !nd(indexType)) runInfo<ivs::d_dna5>();
    else throw error_fmt{"unknown index with {} letters, index type {}", sigma, indexType};
}
}
