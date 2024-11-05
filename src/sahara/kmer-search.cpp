// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "AdaptiveKmerIndex.h"
#include "hash.h"
#include "utils/StopWatch.h"
#include "utils/error_fmt.h"

#include <cereal/types/unordered_map.hpp>
#include <clice/clice.h>
#include <fstream>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <search_schemes/expand.h>
#include <search_schemes/generator/all.h>
#include <search_schemes/nodeCount.h>
#include <unordered_set>

namespace {
void app();
auto cli = clice::Argument{ .args   = "kmer-search",
                            .desc   = "search for a given pattern",
                            .cb     = app,
};

auto cliQuery = clice::Argument{ .parent = &cli,
                                 .args   = "--query",
                                 .desc   = "path to a query file",
                                 .value  = std::filesystem::path{},
};

auto cliIndex = clice::Argument{ .parent = &cli,
                                 .args   = "--index",
                                 .desc   = "path to the index file",
                                 .value  = std::filesystem::path{},
};

auto cliOutput = clice::Argument{ .parent = &cli,
                                  .args   = "--output",
                                  .desc   = "output path",
                                  .value  = std::filesystem::path{"sahara-output.txt"},
};


auto cliGenerator  = clice::Argument{ .parent = &cli,
                                      .args   = "--generator",
                                      .desc   = "picking optimum search scheme generator",
                                      .value  = std::string{"h2-k2"},
};

auto cliDynGenerator = clice::Argument{ .parent = &cli,
                                        .args   = "--dynamic_generator",
                                        .desc   = "should generator run expand search scheme with dynamic extension",
};

auto cliNoReverse = clice::Argument{ .parent = &cli,
                                     .args   = "--no-reverse",
                                     .desc   = "do not search for reversed complements",
};

enum class SearchMode { All, BestHits };
auto cliSearchMode = clice::Argument{ .parent = &cli,
                                      .args   = "--search_mode",
                                      .desc   = "do not search for reversed complements",
                                      .value  = SearchMode::All,
                                      .mapping = {{{"all", SearchMode::All}, {"besthits", SearchMode::BestHits}}},
};
auto cliMaxHits   = clice::Argument{ .parent = &cli,
                                     .args   = "--max_hits",
                                     .desc   = "maximum number of hits per query",
                                     .value  = 0,
};

void app() {
    using Alphabet = ivs::d_dna5;
    constexpr size_t Sigma = Alphabet::size();
    constexpr size_t KmerSigma = 128;

    auto timing = std::vector<std::tuple<std::string, double>>{};

    auto stopWatch = StopWatch();

    fmt::print(
        "config:\n"
        "  query:               {}\n"
        "  index:               {}\n"
        "  generator:           {}\n"
        "  dynamic expansion:   {}\n"
        "  reverse complements: {}\n"
        "  search mode:         {}\n"
        "  max hits:            {}\n"
        "  output path:         {}\n",
        *cliQuery, *cliIndex, *cliGenerator, (bool)cliDynGenerator, !cliNoReverse,
        *cliSearchMode == SearchMode::BestHits?"besthits":"all", *cliMaxHits,
        *cliOutput);


    // load index file
    if (!std::filesystem::exists(*cliIndex)) {
        throw error_fmt{"no valid index path at {}", *cliIndex};
    }

    auto index = AdaptiveKmerIndex{};

    auto uniq = std::unordered_map<size_t, uint8_t>{};
    {
        auto ifs     = std::ifstream{*cliIndex, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        uint32_t fileFormatVersion;
        archive(fileFormatVersion);
        if (fileFormatVersion == 0x01) {
            index.load(archive);
            archive(uniq);
        } else {
            throw error_fmt("unknown file format version for index: {}", fileFormatVersion);
        }
    }
    auto config = index.config();


    fmt::print("  kmer mode:           {}\n", uint8_t(config.mode));
    if (config.mode == AdaptiveKmerIndex::KmerMode::Winnowing) {
        fmt::print("  window:           {}\n", config.window);
    } else if (config.mode == AdaptiveKmerIndex::KmerMode::Mod) {
        fmt::print("  kmer mod:            {}\n", config.modExp);
    } else {
        throw error_fmt("unknown kmer mode {}", uint8_t(config.mode));
    }

    timing.emplace_back("ld index", stopWatch.reset());

    // load fasta file
    auto reader = ivio::fasta::reader {{*cliQuery}};
    size_t totalSize{};
    size_t kmerLen{};
    auto ref_kmer = std::vector<std::vector<uint8_t>>{};
    size_t smallestKmer = std::numeric_limits<size_t>::max();
    size_t longestKmer{};
    size_t skipped{};
    [&]() {
        auto ref = std::vector<uint8_t>{};
        size_t recordNbr = 0;
        for (auto record : reader) {
            recordNbr += 1;
            totalSize += record.seq.size();
            ref.resize(record.seq.size());
            ivs::convert_char_to_rank<Alphabet>(record.seq, ref);
            if (auto pos = ivs::verify_rank(ref); pos) {
                throw error_fmt{"query '{}' ({}) has invalid character at position {} '{}'({:x})", record.id, recordNbr, *pos, record.seq[*pos], record.seq[*pos]};
            }

            [&]() {
                ref_kmer.emplace_back();
                if (config.mode == AdaptiveKmerIndex::KmerMode::Winnowing) {
                    for (auto v : ivs::winnowing_minimizer<Alphabet, /*DuplicatesAllowed=*/false>(ref, /*.k=*/config.kmerLen, /*.window=*/config.window)) {
                        if (auto iter = uniq.find(v); iter != uniq.end()) {
                            ref_kmer.back().emplace_back(iter->second);
                        } else {
                            ref_kmer.pop_back();
                            return;
                        }
                    }
                } else if (config.mode == AdaptiveKmerIndex::KmerMode::Mod) {
                    uint64_t mask = (1<<config.modExp)-1;
                    for (auto v : ivs::compact_encoding<Alphabet, /*UseCanonicalKmers=*/true>{ref, /*.k=*/config.kmerLen}) {
                        v = hash(v);
                        if ((v & mask) != 0) continue;
                        if (auto iter = uniq.find(v); iter != uniq.end()) {
                            ref_kmer.back().emplace_back(iter->second);
                        } else {
                            ref_kmer.pop_back();
                            return;
                        }
                    }
                } else {
                    throw error_fmt("unknown kmer mode: {}", uint8_t(config.mode));
                }
                if (ref_kmer.back().size() >= 6) {
                    smallestKmer = std::min(ref_kmer.back().size(), smallestKmer);
                    longestKmer = std::max(ref_kmer.back().size(), longestKmer);
                    kmerLen += ref_kmer.back().size();
                    if (!cliNoReverse) {
                        ref_kmer.emplace_back(ref_kmer.back());
                        std::ranges::reverse(ref_kmer.back());
                    }
                } else {
                    skipped += 1;
                    ref_kmer.pop_back();
                    if (!cliNoReverse) {
                        skipped += 1;
                    }
                }
            }();
        }
    }();
    fmt::print("skipped {} of {} queries\n", skipped, skipped + ref_kmer.size());
    fmt::print("avg kmer len: {}\n", kmerLen * 1.0/ ref_kmer.size());
    fmt::print("smallest/longest kmer len: {}/{}\n", smallestKmer, longestKmer);
    fmt::print("index uniq {}\n", uniq.size());

    if (ref_kmer.empty()) {
        throw error_fmt{"query file {} was empty - abort\n", *cliQuery};
    }

    {
        auto fwdQueries = ref_kmer.size() / (cliNoReverse?1:2);
        auto bwdQueries = ref_kmer.size() - fwdQueries;
        fmt::print("fwd queries: {}\n"
                   "bwd queries: {}\n",
                   fwdQueries, bwdQueries);
    }
    timing.emplace_back("ld queries", stopWatch.reset());

    auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
    for (size_t qidx{0}; qidx < ref_kmer.size(); ++qidx) {
        index.search(ref_kmer[qidx], [&](size_t refid, size_t refpos) {
            results.emplace_back(qidx, refid, refpos);
        });
    }
    timing.emplace_back("search", stopWatch.reset());

    auto finishTime = std::chrono::steady_clock::now();
    {
        auto ofs = fopen(cliOutput->c_str(), "w");
        for (auto const& [queryId, seqId, pos] : results) {
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
    fmt::print("  queries per second:  {:> 10.0f}q/s\n", ref_kmer.size() / totalTime);
    fmt::print("  number of hits:      {:>10}\n", results.size());
}
}
