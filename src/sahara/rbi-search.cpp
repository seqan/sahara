// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "dr_dna.h"

#include "utils/StopWatch.h"
#include "utils/error_fmt.h"

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <clice/clice.h>
#include <fmindex-collection/suffixarray/DenseCSA.h>
#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/search/all.h>
#include <fstream>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <search_schemes/expand.h>
#include <search_schemes/generator/all.h>
#include <search_schemes/nodeCount.h>
#include <search_schemes/search_schemes.h>
#include <string>
#include <unordered_set>

using namespace fmindex_collection;

namespace {
void app();
auto cli = clice::Argument{ .args   = "rbi-search",
                            .desc   = "search for a given pattern",
                            .cb     = app,
};

auto cliQuery = clice::Argument{ .parent = &cli,
                                 .args   = {"-q", "--query"},
                                 .desc   = "path to a query file",
                                 .value  = std::filesystem::path{},
};

auto cliIndex = clice::Argument{ .parent = &cli,
                                 .args   = {"-i", "--index"},
                                 .desc   = "path to the index file",
                                 .value  = std::filesystem::path{},
};

auto cliOutput = clice::Argument{ .parent = &cli,
                                  .args   = {"-o", "--output"},
                                  .desc   = "output path",
                                  .value  = std::filesystem::path{"sahara-output.txt"},
};


auto cliGenerator  = clice::Argument{ .parent = &cli,
                                      .args   = {"-g", "--generator"},
                                      .desc   = "picking optimum search scheme generator",
                                      .value  = std::string{"h2-k2"},
};

auto cliDynGenerator = clice::Argument{ .parent = &cli,
                                        .args   = "--dynamic_generator",
                                        .desc   = "should generator run expand search scheme with dynamic extension",
};

auto cliNumErrors = clice::Argument{ .parent = &cli,
                                     .args   = {"-e", "--errors"},
                                     .desc   = "number of allowed errors (number of allowed differences insert/substitute and deletions)",
                                     .value  = size_t{},
};

enum class SearchMode { All, BestHits };
auto cliSearchMode = clice::Argument{ .parent = &cli,
                                      .args   = {"-m", "--search_mode"},
                                      .desc   = "search mode, all (default) or besthits",
                                      .value  = SearchMode::All,
                                      .mapping = {{{"all", SearchMode::All}, {"besthits", SearchMode::BestHits}}},
};
auto cliMaxHits   = clice::Argument{ .parent = &cli,
                                     .args   = "--max_hits",
                                     .desc   = "maximum number of hits per query",
                                     .value  = 0,
};

void app() {
    using Alphabet = dr_dna5;
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
        "  search mode:         {}\n"
        "  max hits:            {}\n"
        "  output path:         {}\n",
        *cliQuery, *cliIndex, *cliGenerator, (bool)cliDynGenerator, *cliNumErrors,
        *cliSearchMode == SearchMode::BestHits?"besthits":"all", *cliMaxHits,
        *cliOutput);


    {
        auto fwdQueries = queries.size();
        fmt::print("fwd queries: {}\n",
                   fwdQueries);
    }

    using Table = fmindex_collection::occtable::Interleaved_32<Sigma>;
//    using Table = fmindex_collection::occtable::EprV7<Sigma>;

    if (!std::filesystem::exists(*cliIndex)) {
        throw error_fmt{"no valid index path at {}", *cliIndex};
    }

    auto index = fmindex_collection::RBiFMIndex<Table, fmindex_collection::DenseCSA>{};
    {
        auto ifs     = std::ifstream{*cliIndex, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        archive(index);
    }
    timing.emplace_back("ld index", stopWatch.reset());

    auto k = *cliNumErrors;

    auto generator = [&]() {
        auto iter = search_schemes::generator::all.find(*cliGenerator);
        if (iter == search_schemes::generator::all.end()) {
            throw error_fmt{"unknown search scheme generetaror \"{}\"", *cliGenerator};
        }
        return iter->second.generator;
    }();

    auto loadSearchScheme = [&](int minK, int maxK) {
        auto len = queries[0].size();
        auto oss = generator(minK, maxK, /*unused*/0, /*unused*/0);
        if (!cliDynGenerator) {
            oss = search_schemes::expand(oss, len);
        } else {
            oss = search_schemes::expandByWNC(oss, len, Sigma, index.size());
        }
        fmt::print("node count: {}\n", search_schemes::nodeCount</*Edit=*/true>(oss, Sigma));
        fmt::print("expected node count: {}\n", search_schemes::weightedNodeCount</*Edit=*/true>(oss, Sigma, index.size()));
        return oss;
    };

    auto resultCursors = std::vector<std::tuple<size_t, LeftRBiFMIndexCursor<decltype(index)>, size_t>>{};
    auto res_cb = [&](size_t queryId, auto cursor, size_t errors) {
        resultCursors.emplace_back(queryId, cursor, errors);
    };
    if (*cliSearchMode == SearchMode::All) {
        auto search_scheme  = loadSearchScheme(0, k);
        timing.emplace_back("searchScheme", stopWatch.reset());

        if (*cliMaxHits == 0) search_ng21::search(index, queries, search_scheme, res_cb);
        else                  search_ng21::search_n(index, queries, search_scheme, *cliMaxHits, res_cb);
    } else {
        auto search_schemes = std::vector<decltype(loadSearchScheme(0, k))>{};
        for (size_t j{0}; j<=k; ++j) {
            search_schemes.emplace_back(loadSearchScheme(j, j));
        }
        timing.emplace_back("searchScheme", stopWatch.reset());

        if (*cliMaxHits == 0) search_ng21::search_best(index, queries, search_schemes, res_cb);
        else                  search_ng21::search_best_n(index, queries, search_schemes, *cliMaxHits, res_cb);
    }
    timing.emplace_back("search", stopWatch.reset());

    auto results = std::vector<std::tuple<size_t, size_t, size_t, size_t>>{};
    for (auto const& [queryId, cursor, e] : resultCursors) {
        for (auto [seqId, pos] : LocateLinear{index, cursor}) {
            results.emplace_back(queryId, seqId, pos, e);
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
}
