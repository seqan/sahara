// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "utils/StopWatch.h"
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
#include <fstream>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <mmser/mmser.h>
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
auto cliThreads = clice::Argument {
    .parent = &cli,
    .args   = {"-t", "--threads"},
    .desc   = "number of threads running search in parallel",
    .value  = size_t{1},
};
auto cliCountOnly = clice::Argument {
    .parent = &cli,
    .args   = {"--count-only"},
    .desc   = "only count the number of results without locating them",
};
auto cliPreloadIndex = clice::Argument {
    .parent = &cli,
    .args   = {"--preload-index"},
    .desc   = "load index via copy not mmap",
};
auto cliBatchSize = clice::Argument {
    .parent = &cli,
    .args   = {"--batch_size"},
    .desc   = "numbers of queries processed in each thread",
    .value  = size_t{64},
};
auto cliNoOpt = clice::Argument {
    .parent = &cli,
    .args   = {"--no-opt-zero-error"},
    .desc   = "do not use zero error optimized code (advanced)",
    .tags   = {"advanced"},
};
auto cliNoKStep = clice::Argument {
    .parent = &cli,
    .args   = {"--no-kstep"},
    .desc   = "do not use kstep, even if available",
    .tags   = {"advanced"},
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

    using Index = VarIndex<Alphabet>;
    auto [varIndex, storageManager] = [&]() -> std::tuple<Index, std::unique_ptr<std::any>> {
        if (cliIndex->string().ends_with(".mmser")) {
            if (cliPreloadIndex) {
                return mmser::loadFileStream<Index>(*cliIndex);
            } else {
                return mmser::loadFile<Index>(*cliIndex);
            }
        } else {
            auto varIndex = VarIndex<Alphabet>{};
            auto ifs     = std::ifstream{*cliIndex, std::ios::binary};
            auto archive = cereal::BinaryInputArchive{ifs};
            archive(varIndex);
            return {std::move(varIndex), std::unique_ptr<std::any>{}};
        }
    }();
    fmt::print("  samplingRate: {}\n", varIndex.samplingRate);


    auto revTextIncluded = varIndex.type.ends_with("-rev");
    if (revTextIncluded && !cliNoReverse) {
        queries.resize(queries.size()/2);
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

    auto indexSize = std::visit([](auto const& index) {
        return index.size();
    }, varIndex.vs);

    auto loadSearchScheme = [&](int minK, int maxK, bool edit) {
        auto len = queries[0].size();
        auto oss = generator(minK, maxK, /*unused*/0, /*unused*/0);
        if (edit) {
            if (!cliDynGenerator) {
                oss = fmc::search_scheme::expand(oss, len);
            } else {
                auto partition = optimizeByWNCTopDown</*Edit=*/true>(oss, len, Sigma, indexSize, 1);
                fmt::print("partition: {}\n", partition);
                oss = fmc::search_scheme::expandByWNCTopDown</*Edit=*/true>(oss, len, Sigma, indexSize, 1);
            }
            fmt::print("node count: {}\n", fmc::search_scheme::nodeCount</*Edit=*/true>(oss, Sigma));
            fmt::print("weighted node count: {}\n", fmc::search_scheme::weightedNodeCount</*Edit=*/true>(oss, Sigma, indexSize));
        } else {
            if (!cliDynGenerator) {
                oss = fmc::search_scheme::expand(oss, len);
            } else {
                auto partition = optimizeByWNCTopDown</*Edit=*/false>(oss, len, Sigma, indexSize, 1);
                fmt::print("partition: {}\n", partition);
                oss = fmc::search_scheme::expandByWNCTopDown</*Edit=*/false>(oss, len, Sigma, indexSize, 1);
            }
            fmt::print("node count: {}\n", fmc::search_scheme::nodeCount</*Edit=*/false>(oss, Sigma));
            fmt::print("weighted node count: {}\n", fmc::search_scheme::weightedNodeCount</*Edit=*/false>(oss, Sigma, indexSize));

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
                partition = optimizeByWNCTopDown</*Edit=*/true>(oss, len, Sigma, indexSize, 1);
                fmt::print("partition: {}\n", partition);
            }
            fmt::print("node count: {}\n", fmc::search_scheme::nodeCount</*Edit=*/true>(oss, Sigma));
            fmt::print("weighted node count: {}\n", fmc::search_scheme::weightedNodeCount</*Edit=*/true>(oss, Sigma, indexSize));
        } else {
            if (!cliDynGenerator) {
            } else {
                partition = optimizeByWNCTopDown</*Edit=*/false>(oss, len, Sigma, indexSize, 1);
                fmt::print("partition: {}\n", partition);
            }
            fmt::print("node count: {}\n", fmc::search_scheme::nodeCount</*Edit=*/false>(oss, Sigma));
            fmt::print("weighted node count: {}\n", fmc::search_scheme::weightedNodeCount</*Edit=*/false>(oss, Sigma, indexSize));

        }
        return std::make_tuple(oss, partition);
    };

    std::visit([&]<typename Index>(Index const& index) {
        auto resultCursors = channel::value_mutex<std::vector<std::tuple<size_t, fmc::LeftBiFMIndexCursor<Index>, size_t>>>{};

        bool Edit = *cliDistanceMetric == DistanceMetric::Levenshtein;
        auto totalHits = size_t{};
        auto res_cb = [&](size_t queryId, auto const& cursor, size_t errors) {
            resultCursors->emplace_back(queryId, cursor, errors);
            totalHits += cursor.count();
        };
        if (*cliSearchMode == SearchMode::All) {
            if (k == 0 && *cliMaxHits == 0 && !cliNoOpt) {
                auto gqidx = channel::value_mutex<size_t>{};
                auto workers = channel::workers{*cliThreads, [&]() {
                    while (true) {
                        auto [qidx, qidx2] = [&]() -> std::tuple<size_t, size_t> {
                            auto [_, ptr] = *gqidx;
                            auto v = *ptr;
                            *ptr = std::min((*ptr)+1024, queries.size());
                            return {v, *ptr};
                        }();
                        if (qidx == qidx2) return;
                        auto report = [&](size_t queryId, auto const& cursor) {
                            res_cb(queryId+qidx, cursor, 0);
                        };
                        auto sub_queries = std::span{queries.begin()+qidx, queries.begin()+qidx2};
                        fmc::search_no_errors::search(index, sub_queries, report, *cliBatchSize);
                    }
                }};
            } else {
                auto [search_scheme, partition]  = loadSearchSchemeUsePartition(0, k, Edit);
                timing.emplace_back("searchScheme", stopWatch.reset());
                size_t maxHits = (*cliMaxHits > 0)?*cliMaxHits:std::numeric_limits<size_t>::max();

                auto gqidx = channel::value_mutex<size_t>{};
                auto workers = channel::workers{*cliThreads, [&]() {
                    while (true) {
                        auto [qidx, qidx2] = [&]() -> std::tuple<size_t, size_t> {
                            auto [_, ptr] = *gqidx;
                            auto v = *ptr;
                            *ptr = std::min((*ptr)+1024, queries.size());
                            return {v, *ptr};
                        }();
                        if (qidx == qidx2) return;
                        auto report = [&](size_t queryId, auto const& cursor, size_t e) {
                            res_cb(queryId+qidx, cursor, e);
                        };
                        auto sub_queries = std::span{queries.begin()+qidx, queries.begin()+qidx2};
                        if (cliNoKStep) {
                            if (!Edit) fmc::search_ng28/*_kstep*/::search<false>(index, sub_queries, search_scheme, partition, res_cb, maxHits);
                            else       fmc::search_ng28/*_kstep*/::search<true >(index, sub_queries, search_scheme, partition, report, maxHits);
                        } else {
                            if (!Edit) fmc::search_ng28_kstep::search<false>(index, sub_queries, search_scheme, partition, res_cb, maxHits);
                            else       fmc::search_ng28_kstep::search<true >(index, sub_queries, search_scheme, partition, report, maxHits);
                        }

                    }
                }};
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
        results.reserve(totalHits);
        auto [_, _resultCursors] = *resultCursors;
        size_t totalNumberOfHits{};

        for (auto const& [queryId, cursor, e] : *_resultCursors) {
            if (cliCountOnly) {
                totalNumberOfHits += cursor.count();
                continue;
            }

            if constexpr (std::same_as<typename Index::ADEntry, std::tuple<uint32_t, uint32_t, bool>>) {
                for (auto [seqId, seqPos, rev, offset] : fmc::LocateLinear{index, cursor}) {
                    if (cliNoReverse && rev) continue;
                    if (!rev) {
                        results.emplace_back(queryId, seqId, seqPos+offset, e);
                    } else {
                        // Found the position inside the reversed SA
                        results.emplace_back(queryId, seqId, seqPos-offset-cursor.steps+1, e);
                    }
//                    fmt::print("match: qid: {}, seqid: {} pos: {}+{}+{}, rev: {}\n", queryId, seqId, seqPos, offset, cursor.steps, rev);
                }
            } else {
                for (auto [seqId, seqPos, offset] : fmc::LocateLinear{index, cursor}) {
                    results.emplace_back(queryId, seqId, seqPos+offset, e);
                }
            }
        }

        if (!cliCountOnly) {
            timing.emplace_back("locate", stopWatch.reset());
        } else {
            timing.emplace_back("count", stopWatch.reset());
        }

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
        if (!cliCountOnly) {
            fmt::print("  number of hits:      {:>10}\n", results.size());
        } else {
            fmt::print("  number of hits:      {:>10}\n", totalNumberOfHits);
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
    if (sigma == 2 && nd(indexType)) runSearch<ivs::dna2>();
    else if (sigma == 3 && !nd(indexType)) runSearch<ivs::d_dna2>();
    else if (sigma == 4 &&  nd(indexType)) runSearch<ivs::dna4>();
    else if (sigma == 5 && !nd(indexType)) runSearch<ivs::d_dna4>();
    else if (sigma == 5 &&  nd(indexType)) runSearch<ivs::dna5>();
    else if (sigma == 6 && !nd(indexType)) runSearch<ivs::d_dna5>();
    else throw error_fmt{"unknown index with {} letters, index type {}", sigma, indexType};
}
}
