#include "error_fmt.h"
#include "utils.h"

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <clice/clice.h>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/search/all.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <fstream>
#include <search_schemes/expand.h>
#include <search_schemes/generator/all.h>
#include <search_schemes/nodeCount.h>
#include <unordered_set>

using namespace fmindex_collection;

namespace {
void app();
auto cli = clice::Argument{ .arg    = "search",
                            .desc   = "search for a given pattern",
                            .cb     = app,
};

auto cliQuery = clice::Argument{ .parent = &cli,
                                 .arg    = "--query",
                                 .desc   = "path to a query file",
                                 .value  = std::filesystem::path{},
};

auto cliIndex = clice::Argument{ .parent = &cli,
                                 .arg    = "--index",
                                 .desc   = "path to the index file",
                                 .value  = std::filesystem::path{},
};

auto cliOutput = clice::Argument{ .parent = &cli,
                                  .arg    = "--output",
                                  .desc   = "output path",
                                  .value  = std::filesystem::path{"sahara-output.txt"},
};


auto cliGenerator  = clice::Argument{ .parent = &cli,
                                      .arg    = "--generator",
                                      .desc   = "picking optimum search scheme generator",
                                      .value  = std::string{"h2-k2"},
};

auto cliDynGenerator = clice::Argument{ .parent = &cli,
                                        .arg    = "--dynamic_generator",
                                        .desc   = "should generator run expand search scheme with dynamic extension",
};

auto cliNumErrors = clice::Argument{ .parent = &cli,
                                     .arg    = "--errors",
                                     .desc   = "number of allowed errors (number of allowed differences insert/substitute and deletions)",
                                     .value  = size_t{},
};
auto cliNoReverse = clice::Argument{ .parent = &cli,
                                     .arg    = "--no-reverse",
                                     .desc   = "do not search for reversed complements",
};

enum class SearchMode { All, BestHits };
auto cliSearchMode = clice::Argument{ .parent = &cli,
                                      .arg    = "--search_mode",
                                      .desc   = "do not search for reversed complements",
                                      .value  = SearchMode::All,
                                      .mapping = {{{"all", SearchMode::All}, {"besthits", SearchMode::BestHits}}},
};
auto cliMaxHits   = clice::Argument{ .parent = &cli,
                                     .arg    = "--max_hits",
                                     .desc   = "maximum number of hits per query",
                                     .value  = 0,
};

void app() {
    constexpr size_t Sigma = 5;

    auto timing = std::vector<std::tuple<std::string, double>>{};

    auto stopWatch = StopWatch();

    auto const [queries, queryInfos] = loadQueries<Sigma>(*cliQuery, !cliNoReverse);
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
        *cliSearchMode == SearchMode::BestHits?"besthits":"all", *cliMaxHits,
        *cliOutput);


    {
        auto fwdQueries = queries.size() / (cliNoReverse?1:2);
        auto bwdQueries = fwdQueries * (cliNoReverse?0:1);
        fmt::print("fwd queries: {}\n"
                   "bwd queries: {}\n",
                   fwdQueries, bwdQueries);
    }

    using Table = fmindex_collection::occtable::interleaved32::OccTable<Sigma>;
//    using Table = fmindex_collection::occtable::interleavedEPRV7::OccTable<Sigma>;

    if (!std::filesystem::exists(*cliIndex)) {
        throw error_fmt{"no valid index path at {}", *cliIndex};
    }

    auto index = fmindex_collection::BiFMIndex<Table, fmindex_collection::DenseCSA>{fmindex_collection::cereal_tag{}};
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
        return iter->second;
    }();

    auto loadSearchScheme = [&](int minK, int maxK) {
        auto len = queries[0].size();
        auto oss = generator(minK, maxK, /*unused*/0, /*unused*/0);
        if (!cliDynGenerator) {
            oss = search_schemes::expand(oss, len);
        } else {
            oss = search_schemes::expandDynamic(oss, len, Sigma, index.size());
        }
        fmt::print("node count: {}\n", search_schemes::nodeCount(oss, Sigma));
        fmt::print("expected node count: {}\n", search_schemes::expectedNodeCount(oss, Sigma, index.size()));
        return oss;
    };

    auto resultCursors = std::vector<std::tuple<size_t, LeftBiFMIndexCursor<decltype(index)>, size_t>>{};
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
