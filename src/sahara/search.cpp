#include "utils.h"

#include <clice/clice.h>
#include <cstdio>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/search/all.h>
#include <fmt/format.h>
#include <search_schemes/expand.h>
#include <search_schemes/generator/all.h>
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

    auto const [queries, queryInfos] = loadQueries<Sigma>(*cliQuery, !cliNoReverse);
    if (queries.empty()) {
        fmt::print("query file was empty - abort\n");
        return;
    }
    fmt::print("loaded {} queries (incl reverse complements)\n", queries.size());
    fmt::print("{:15}: {:>10}  ({:>10} +{:>10} ) {:>10}    - results: {:>10}/{:>10}/{:>10}/{:>10} - mem: {:>13}\n", "name", "time_search + time_locate", "time_search", "time_locate", "(time_search+time_locate)/queries.size()", "resultCt", "results.size()", "uniqueResults.size()", "readIds.size()", "memory");

    using Table = fmindex_collection::occtable::interleaved32::OccTable<Sigma>;
//    using Table = fmindex_collection::occtable::interleavedEPRV7::OccTable<Sigma>;
    auto name = Table::extension();
    fmt::print("start loading {} ...", name);
    fflush(stdout);
    auto index = loadDenseIndex<CSA, Table>(*cliIndex, /*samplingRate*/16, /*threadNbr*/1);
    fmt::print("done\n");

    auto memory = [&] () -> size_t {
        if constexpr (OccTableMemoryUsage<Table>) {
            return index.memoryUsage();
        } else {
            return 0ull;
        }
    }();
    auto k = *cliNumErrors;
    auto generator = *cliGenerator;
    auto search_scheme = [&]() {
        auto iter = search_schemes::generator::all.find(generator);
        if (iter == search_schemes::generator::all.end()) {
            throw std::runtime_error("unknown search scheme generetaror \"" + generator + "\"");
        }
        auto len = queries[0].size();
        auto oss = iter->second(0, k, 0, 0); //!TODO last two parameters of second are not being used
        auto ess = search_schemes::expand(oss, len);
        auto dss = search_schemes::expandDynamic(oss, len, 4, 3'000'000'000); //!TODO use correct Sigma and text size
        fmt::print("ss diff: {} to {}, using dyn: {}\n", search_schemes::expectedNodeCount(ess, 4, 3'000'000'000), search_schemes::expectedNodeCount(dss, 4, 3'000'000'000), cliDynGenerator);
        if (!cliDynGenerator) {
            return ess;
        } else {
            return dss;
        }
    }();
    auto search_schemes = [&]() {
        auto r = std::vector<decltype(search_scheme)>{};
        for (size_t j{0}; j<=k; ++j) {
            r.emplace_back([&]() {
                auto iter = search_schemes::generator::all.find(generator);
                if (iter == search_schemes::generator::all.end()) {
                    throw std::runtime_error("unknown search scheme generetaror \"" + generator + "\"");
                }
                auto len = queries[0].size();
                auto oss = iter->second(j, j, 0, 0); //!TODO last two parameters of second are not being used
                auto ess = search_schemes::expand(oss, len);
                auto dss = search_schemes::expandDynamic(oss, len, 4, 3'000'000'000); //!TODO use correct Sigma and text size
                if (!cliDynGenerator) {
                    return ess;
                } else {
                    return dss;
                }
            }());
        }
        return r;
    }();

    size_t resultCt{};
    StopWatch sw;
    auto results       = std::vector<std::tuple<size_t, size_t, size_t, size_t>>{};
    auto resultCursors = std::vector<std::tuple<size_t, LeftBiFMIndexCursor<decltype(index)>, size_t>>{};
    auto resultCursorsEditTranscript = std::vector<std::string>{};

    auto res_cb = [&](size_t queryId, auto cursor, size_t errors) {
        resultCursors.emplace_back(queryId, cursor, errors);
    };
    auto res_cb2 = [&](size_t queryId, auto cursor, size_t errors, auto const& actions) {
        std::string s;
        for (auto a : actions) {
            s += a;
        }
        resultCursors.emplace_back(queryId, cursor, errors);
        resultCursorsEditTranscript.emplace_back(std::move(s));
    };

    if (*cliSearchMode == SearchMode::All) {
        if (*cliMaxHits == 0) search_ng21::search(index, queries, search_scheme, res_cb);
        else                  search_ng21::search_n(index, queries, search_scheme, *cliMaxHits, res_cb);
    } else if (*cliSearchMode == SearchMode::BestHits) {
        if (*cliMaxHits == 0) search_ng21::search_best(index, queries, search_schemes, res_cb);
        else                  search_ng21::search_best_n(index, queries, search_schemes, *cliMaxHits, res_cb);
    }

    auto time_search = sw.reset();

    for (auto const& [queryId, cursor, e] : resultCursors) {
        for (auto [seqId, pos] : LocateLinear{index, cursor}) {
            results.emplace_back(queryId, seqId, pos, e);
        }
        resultCt += cursor.len;
    }
    auto time_locate = sw.reset();

    auto uniqueResults = [](auto list) {
        std::sort(begin(list), end(list));
        list.erase(std::unique(begin(list), end(list)), list.end());
        return list;
    }(results);
    auto readIds = std::unordered_set<size_t>{};
    for (auto const& [queryId, cursor, e] : resultCursors) {
        if (queryId > queries.size()/2) {
            readIds.insert(queryId - queries.size() / 2);
        } else {
            readIds.insert(queryId);
        }
    }

    fmt::print("{:15} {:3}: {:>10.3}s ({:>10.3}s+{:>10.3}s) {:>10.3}q/s - results: {:>10}/{:>10}/{:>10}/{:>10} - mem: {:>13}\n", name, k, time_search + time_locate, time_search, time_locate, queries.size() / (time_search+time_locate), resultCt, results.size(), uniqueResults.size(), readIds.size(), memory);
    {
//        if (!config.saveOutput.empty()) {
//            auto ofs = fopen(config.saveOutput.string().c_str(), "w");
//            for (auto const& [queryId, seqId, pos, e] : results) {
//                fmt::print(ofs, "{} {} {}\n", queryId, seqId, pos);
//            }
//            fclose(ofs);
//        }
    }
}
}
