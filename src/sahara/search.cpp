// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "utils/StopWatch.h"
#include "utils/error_fmt.h"

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <clice/clice.h>
#include <fmindex-collection/suffixarray/DenseCSA.h>
#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/search/SearchNg26.h>
#include <fstream>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <string>
#include <unordered_set>


namespace {
void app();
auto cli = clice::Argument {
    .args   = "search",
    .desc   = "search for a given pattern",
    .cb     = app,
};

auto cliQuery = clice::Argument {
    .parent = &cli,
    .args   = {"-q", "--query"},
    .desc   = "path to a query file",
    .value  = std::filesystem::path{},
};

auto cliIndex = clice::Argument {
    .parent = &cli,
    .args   = {"-i", "--index"},
    .desc   = "path to the index file",
    .value  = std::filesystem::path{},
};

auto cliOutput = clice::Argument {
    .parent = &cli,
    .args   = {"-o", "--output"},
    .desc   = "output path",
    .value  = std::filesystem::path{"sahara-output.txt"},
};

auto cliGenerator  = clice::Argument {
    .parent = &cli,
    .args   = {"-g", "--generator"},
    .desc   = "picking optimum search scheme generator",
    .value  = std::string{"h2-k2"},
};

auto cliDynGenerator = clice::Argument {
    .parent = &cli,
    .args   = "--dynamic_generator",
    .desc   = "should generator run expand search scheme with dynamic extension",
};

auto cliNumErrors = clice::Argument {
    .parent = &cli,
    .args   = {"-e", "--errors"},
    .desc   = "number of allowed errors (number of allowed differences insert/substitute and deletions)",
    .value  = size_t{},
};
auto cliNoReverse = clice::Argument {
    .parent = &cli,
    .args   = "--no-reverse",
    .desc   = "do not search for reversed complements",
};

enum class SearchMode { All, BestHits };
auto cliSearchMode = clice::Argument {
    .parent  = &cli,
    .args    = {"-m", "--search_mode"},
    .desc    = "search mode, all (default) or besthits",
    .value   = SearchMode::All,
    .mapping = {{{"all", SearchMode::All}, {"besthits", SearchMode::BestHits}}},
};
enum class DistanceMetric { Hamming, Levenshtein };
auto cliDistanceMetric = clice::Argument {
    .parent  = &cli,
    .args    = {"-d", "--distance-metric"},
    .desc    = "which distance metric to use. ham: hamming or lev: levenshtein(edit) distance",
    .value   = DistanceMetric::Levenshtein,
    .mapping = {{{"ham", DistanceMetric::Hamming}, {"lev", DistanceMetric::Levenshtein}}}
};
auto cliMaxHits = clice::Argument {
    .parent = &cli,
    .args   = "--max_hits",
    .desc   = "maximum number of hits per query",
    .value  = 0,
};
auto cliLimitQueries = clice::Argument {
    .parent = &cli,
    .args   = {"--limit_queries"},
    .desc   = "only run the given number of queries",
    .value  = size_t{},
};

template <typename Alphabet>
void runSearch() {
    constexpr size_t Sigma = Alphabet::size();

    auto timing = std::vector<std::tuple<std::string, double>>{};

    auto stopWatch = StopWatch();

    // load fasta file
    size_t totalSize{};
    auto queries = std::vector<std::vector<uint8_t>>{};
    for (auto record : ivio::fasta::reader {{*cliQuery}}) {
        totalSize += record.seq.size();
        queries.emplace_back(ivs::convert_char_to_rank<Alphabet>(record.seq));
        if (auto pos = ivs::verify_rank(queries.back()); pos) {
            throw error_fmt{"query '{}' ({}) has invalid character at position {} '{}'({:x})", record.id, queries.size(), *pos, record.seq[*pos], record.seq[*pos]};
        }
        if (!cliNoReverse) {
            queries.emplace_back(ivs::reverse_complement_rank<Alphabet>(queries.back()));
        }
    }
    if (cliLimitQueries) {
        queries.resize(std::min(*cliLimitQueries, queries.size()));
    }
    if (queries.empty()) {
        throw error_fmt{"query file {} was empty - abort\n", *cliQuery};
    }
    timing.emplace_back("ld queries", stopWatch.reset());


    fmt::print(
        "config:\n"
        "  query:               {}\n"
        "  index:               {}\n"
        "  generator:           {}\n"
        "  dynamic expansion:   {}\n"
        "  allowed errors:      {}\n"
        "  reverse complements: {}\n"
        "  search mode:         {}\n"
        "  max hits:            {}\n"
        "  output path:         {}\n",
        *cliQuery, *cliIndex, *cliGenerator, (bool)cliDynGenerator, *cliNumErrors, !cliNoReverse,
        (*cliSearchMode == SearchMode::BestHits?"besthits":"all"), *cliMaxHits,
        *cliOutput);


    {
        auto fwdQueries = queries.size() / (cliNoReverse?1:2);
        auto bwdQueries = queries.size() - fwdQueries;
        fmt::print("fwd queries: {}\n"
                   "bwd queries: {}\n",
                   fwdQueries, bwdQueries);
    }

    if (!std::filesystem::exists(*cliIndex)) {
        throw error_fmt{"no valid index path at {}", *cliIndex};
    }

//    auto index = fmc::BiFMIndex<Sigma/*, fmc::string::InterleavedBitvectorPrefix16*/>{};
    auto index = fmc::BiFMIndex<Sigma, fmc::string::FlattenedBitvectors_64_64k>{};

    {
        auto ifs     = std::ifstream{*cliIndex, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        size_t sigma;
        archive(sigma);
        archive(index);
    }
    timing.emplace_back("ld index", stopWatch.reset());

    auto k = *cliNumErrors;

    auto generator = [&]() {
        auto iter = fmc::search_scheme::generator::all.find(*cliGenerator);
        if (iter == fmc::search_scheme::generator::all.end()) {
            auto names = std::vector<std::string>{};
            for (auto const& [key, gen] : fmc::search_scheme::generator::all) {
                names.push_back(key);
            }
            throw error_fmt{"unknown search scheme generetaror \"{}\", valid generators are: {}", *cliGenerator, fmt::join(names, ", ")};
        }
        return iter->second.generator;
    }();

    auto loadSearchScheme = [&](int minK, int maxK, bool edit) {
        auto len = queries[0].size();
        auto oss = generator(minK, maxK, /*unused*/0, /*unused*/0);
        if (edit) {
            if (!cliDynGenerator) {
                oss = fmc::search_scheme::expand(oss, len);
            } else {
                auto partition = optimizeByWNCTopDown</*Edit=*/true>(oss, len, Sigma, index.size(), 1);
                fmt::print("partition: {}\n", partition);
                oss = fmc::search_scheme::expandByWNCTopDown</*Edit=*/true>(oss, len, Sigma, index.size(), 1);
            }
            fmt::print("node count: {}\n", fmc::search_scheme::nodeCount</*Edit=*/true>(oss, Sigma));
            fmt::print("weighted node count: {}\n", fmc::search_scheme::weightedNodeCount</*Edit=*/true>(oss, Sigma, index.size()));
        } else {
            if (!cliDynGenerator) {
                oss = fmc::search_scheme::expand(oss, len);
            } else {
                auto partition = optimizeByWNCTopDown</*Edit=*/false>(oss, len, Sigma, index.size(), 1);
                fmt::print("partition: {}\n", partition);
                oss = fmc::search_scheme::expandByWNCTopDown</*Edit=*/false>(oss, len, Sigma, index.size(), 1);
            }
            fmt::print("node count: {}\n", fmc::search_scheme::nodeCount</*Edit=*/false>(oss, Sigma));
            fmt::print("weighted node count: {}\n", fmc::search_scheme::weightedNodeCount</*Edit=*/false>(oss, Sigma, index.size()));

        }
        return oss;
    };

    auto loadSearchSchemeUsePartition = [&](int minK, int maxK, bool edit) {
        auto len       = queries[0].size();
        auto oss       = generator(minK, maxK, /*unused*/0, /*unused*/0);
        auto partition = fmc::search_scheme::createUniformPartition(oss, queries[0].size());
        if (edit) {
            if (!cliDynGenerator) {
            } else {
                partition = optimizeByWNCTopDown</*Edit=*/true>(oss, len, Sigma, index.size(), 1);
                fmt::print("partition: {}\n", partition);
            }
            fmt::print("node count: {}\n", fmc::search_scheme::nodeCount</*Edit=*/true>(oss, Sigma));
            fmt::print("weighted node count: {}\n", fmc::search_scheme::weightedNodeCount</*Edit=*/true>(oss, Sigma, index.size()));
        } else {
            if (!cliDynGenerator) {
            } else {
                partition = optimizeByWNCTopDown</*Edit=*/false>(oss, len, Sigma, index.size(), 1);
                fmt::print("partition: {}\n", partition);
            }
            fmt::print("node count: {}\n", fmc::search_scheme::nodeCount</*Edit=*/false>(oss, Sigma));
            fmt::print("weighted node count: {}\n", fmc::search_scheme::weightedNodeCount</*Edit=*/false>(oss, Sigma, index.size()));

        }
        return std::make_tuple(oss, partition);
    };


    auto resultCursors = std::vector<std::tuple<size_t, fmc::LeftBiFMIndexCursor<decltype(index)>, size_t>>{};


    bool Edit = *cliDistanceMetric == DistanceMetric::Levenshtein;
    auto res_cb = [&](size_t queryId, auto const& cursor, size_t errors) {
        resultCursors.emplace_back(queryId, cursor, errors);
    };
    if (*cliSearchMode == SearchMode::All) {
        auto [search_scheme, partition]  = loadSearchSchemeUsePartition(0, k, Edit);
        timing.emplace_back("searchScheme", stopWatch.reset());

        if (!Edit) {
            if (*cliMaxHits == 0) fmc::search_ng26::search<false>(index, queries, search_scheme, partition, res_cb);
            else                  fmc::search_ng26::search<false>(index, queries, search_scheme, partition, res_cb, *cliMaxHits);
        } else {
            if (*cliMaxHits == 0) fmc::search_ng26::search<true>(index, queries, search_scheme, partition, res_cb);
            else                  fmc::search_ng26::search<true>(index, queries, search_scheme, partition, res_cb, *cliMaxHits);
        }
    } else {
        auto search_schemes = std::vector<decltype(loadSearchSchemeUsePartition(0, k, Edit))>{};
        for (size_t j{0}; j<=k; ++j) {
            search_schemes.emplace_back(loadSearchSchemeUsePartition(j, j, Edit));
        }
        timing.emplace_back("searchScheme", stopWatch.reset());
        if (*cliMaxHits == 0) fmc::search_ng26::search_best(index, queries, search_schemes, res_cb);
        else                  fmc::search_ng26::search_best(index, queries, search_schemes, res_cb, *cliMaxHits);
    }
    timing.emplace_back("search", stopWatch.reset());

    auto results = std::vector<std::tuple<size_t, size_t, size_t, size_t>>{};
    for (auto const& [queryId, cursor, e] : resultCursors) {
        for (auto [sae, offset] : fmc::LocateLinear{index, cursor}) {
            auto [seqId, seqPos] = sae;
            results.emplace_back(queryId, seqId, seqPos+offset, e);
        }
    }

    timing.emplace_back("locate", stopWatch.reset());

    auto finishTime = std::chrono::steady_clock::now();
    {
        auto ofs = fopen(cliOutput->c_str(), "w");
        for (auto const& [queryId, seqId, pos, e] : results) {
            fmt::print(ofs, "{} {} {}\n", queryId, seqId, pos);
        }
        fclose(ofs);
    }

    timing.emplace_back("result", stopWatch.reset());

    fmt::print("stats:\n");
    double totalTime{};
    for (auto const& [key, time] : timing) {
        fmt::print("  {:<20} {:> 10.2f}s\n", key + " time:", time);
        totalTime += time;
    }
    fmt::print("  total time:          {:> 10.2f}s\n", totalTime);
    fmt::print("  queries per second:  {:> 10.0f}q/s\n", queries.size() / totalTime);
    fmt::print("  number of hits:      {:>10}\n", results.size());
}

void app() {
    // load sigma value
    size_t sigma;
    {
        auto ifs     = std::ifstream{*cliIndex, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        archive(sigma);
    }
    if (sigma == 5) {
        runSearch<ivs::d_dna4>();
    } else if (sigma == 6) {
        runSearch<ivs::d_dna5>();
    } else {
        throw error_fmt{"unknown index with {} letters", sigma};
    }
}
}
