#include "error_fmt.h"
#include "utils.h"

#include <cereal/archives/binary.hpp>
#include "cereal/types/unordered_map.hpp"
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
auto cli = clice::Argument{ .arg    = "kmer-search",
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
    using Alphabet = ivs::d_dna5;
    constexpr size_t Sigma = Alphabet::size();
    constexpr size_t KmerSigma = 256;

    auto timing = std::vector<std::tuple<std::string, double>>{};

    auto stopWatch = StopWatch();

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


    // load index file
    //using Table = fmindex_collection::occtable::interleaved32::OccTable<Sigma>;
    using Table = fmindex_collection::occtable::interleavedEPRV7::OccTable<Sigma>;
    using KmerTable = fmindex_collection::occtable::interleavedEPRV7_32::OccTable<KmerSigma>;


    if (!std::filesystem::exists(*cliIndex)) {
        throw error_fmt{"no valid index path at {}", *cliIndex};
    }

    auto index = fmindex_collection::BiFMIndex<Table, fmindex_collection::DenseCSA>{fmindex_collection::cereal_tag{}};
    auto kmerIndex = fmindex_collection::BiFMIndex_32<KmerTable, fmindex_collection::CSA_32>{fmindex_collection::cereal_tag{}};
    size_t kmer{};
    size_t window{};
    auto uniq = std::unordered_map<size_t, uint32_t>{};
    {
        auto ifs     = std::ifstream{*cliIndex, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        archive(index, kmerIndex, kmer, window, uniq);
    }
    timing.emplace_back("ld index", stopWatch.reset());

    // load fasta file
    auto reader = ivio::fasta::reader {{*cliQuery}};
    size_t totalSize{};
    size_t kmerLen{};
    auto ref = std::vector<std::vector<uint8_t>>{};
    auto ref_kmer = std::vector<std::vector<uint32_t>>{};
    auto uniq2 = uniq;
    size_t smallestKmer = std::numeric_limits<size_t>::max();
    size_t longestKmer{};
    for (auto record : reader) {
        totalSize += record.seq.size();
        ref.emplace_back(ivs::convert_char_to_rank<ivs::d_dna5>(record.seq));
        auto seq = ivs::convert_char_to_rank<ivs::dna5>(record.seq);
        if (!ivs::verify_rank(seq)) {
            throw std::runtime_error{"something went wrong"};
        }

        [&]() {
            ref_kmer.emplace_back();
            for (auto v : ivs::winnowing_minimizer<ivs::dna5>(seq, /*.k=*/kmer, /*.window=*/window)) {
                if (auto iter = uniq.find(v); iter != uniq.end()) {
                    ref_kmer.back().emplace_back(iter->second);
                } else {
                    ref_kmer.pop_back();
                    ref.pop_back();
                    return;
                }
            }
            smallestKmer = std::min(ref_kmer.back().size(), smallestKmer);
            longestKmer = std::max(ref_kmer.back().size(), longestKmer);
            kmerLen += ref_kmer.back().size();

            if (!cliNoReverse) {
                ref.emplace_back(ivs::reverse_complement_rank<ivs::d_dna5>(ref.back()));
                ref_kmer.emplace_back(ref_kmer.back());
                std::ranges::reverse(ref_kmer.back());
            }
        }();
    }
    fmt::print("avg kmer len: {}\n", kmerLen * 1.0/ ref_kmer.size());
    fmt::print("smallest/longest kmer len: {}/{}\n", smallestKmer, longestKmer);
    fmt::print("index uniq {}, query uniq {}\n", uniq2.size(), uniq.size());

    if (ref.empty()) {
        throw error_fmt{"query file {} was empty - abort\n", *cliQuery};
    }

    {
        auto fwdQueries = ref.size() / (cliNoReverse?1:2);
        auto bwdQueries = fwdQueries * (cliNoReverse?0:1);
        fmt::print("fwd queries: {}\n"
                   "bwd queries: {}\n",
                   fwdQueries, bwdQueries);
    }
    timing.emplace_back("ld queries", stopWatch.reset());

    auto k = *cliNumErrors;

    auto generator = [&]() {
        auto iter = search_schemes::generator::all.find(*cliGenerator);
        if (iter == search_schemes::generator::all.end()) {
            throw error_fmt{"unknown search scheme generetaror \"{}\"", *cliGenerator};
        }
        return iter->second;
    }();

    auto loadSearchScheme = [&](int minK, int maxK) {
        auto len = ref[0].size();
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

    auto loadKmerSearchScheme = [&](int minK, int maxK, int len) {
        auto oss = generator(minK, maxK, /*unused*/0, /*unused*/0);
        if (!cliDynGenerator) {
            oss = search_schemes::expand(oss, len);
        } else {
            oss = search_schemes::expandDynamic(oss, len, Sigma, index.size());
        }
        for (auto& os : oss) {
            for (size_t i{0}; i < os.pi.size(); ++i) {
                if (os.u[i] > 0) {
                    os.pi.erase(os.pi.begin()+i, os.pi.end());
                    os.l.erase(os.l.begin()+i, os.l.end());
                    os.u.erase(os.u.begin()+i, os.u.end());
                    break;
                }
            }
        }
//        fmt::print("node count: {}\n", search_schemes::nodeCount(oss, Sigma));
//        fmt::print("expected node count: {}\n", search_schemes::expectedNodeCount(oss, Sigma, index.size()));
        return oss;
    };


    auto resultCursors = std::vector<std::tuple<size_t, LeftBiFMIndexCursor<decltype(index)>, size_t>>{};
    auto res_cb = [&](size_t queryId, auto cursor, size_t errors) {
        resultCursors.emplace_back(queryId, cursor, errors);
    };
//    if (*cliSearchMode == SearchMode::All) {
        auto search_schemes  = loadSearchScheme(0, k);
        auto reordered_list  = search_ng21::prepare_reorder(search_schemes);

        auto kmer_search_schemes_by_len = std::vector<decltype(search_schemes)>{};
        auto kmer_reordered_list_by_len = std::vector<decltype(reordered_list)>{};
        kmer_search_schemes_by_len.resize(longestKmer+1);
        kmer_reordered_list_by_len.resize(longestKmer+1);
        for (size_t i{smallestKmer}; i <= longestKmer; ++i){
            kmer_search_schemes_by_len[i] = loadKmerSearchScheme(0, k, i);
            kmer_reordered_list_by_len[i] = search_ng21::prepare_reorder(kmer_search_schemes_by_len[i]);
        }

        timing.emplace_back("searchScheme", stopWatch.reset());

        size_t totalKmerHits{};
        for (size_t qidx{0}; qidx < ref.size(); ++qidx) {
            auto& kmer_search_schemes = kmer_search_schemes_by_len[ref_kmer[qidx].size()];
            auto& kmer_reordered_list = kmer_reordered_list_by_len[ref_kmer[qidx].size()];

            for (size_t i{0}; i < reordered_list.size(); ++i) {
                auto& kmer_reordered     = kmer_reordered_list[i];
                auto& kmer_search_scheme = kmer_search_schemes[i];
                bool found{false};
                auto& kmer_search = kmer_reordered;
                #if 1
                for (size_t k {0}; k < kmer_search.size(); ++k) {
                    kmer_search[k].rank = ref_kmer[qidx][kmer_search_scheme.pi[k]];
                    assert(kmer_search[k].rank != 0);
                }
                search_ng21::Search{kmerIndex, kmer_search, [&](auto cur, size_t e) {
                    totalKmerHits += cur.count();
                    found = true;
                    return std::true_type{};
                }}.run();
                //continue;
                #else
                found = true;
                #endif

                if (found) {
                    auto& search        = reordered_list[i];
                    auto& search_scheme = search_schemes[i];
                    for (size_t k {0}; k < search.size(); ++k) {
                        search[k].rank = ref[qidx][search_scheme.pi[k]];
                    }
                    search_ng21::Search{index, search, [&](auto cur, size_t e) {
                        res_cb(qidx, cur, e);
                        return std::false_type{};
                    }}.run();
                }
            }
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
    fmt::print("  queries per second:  {:> 10.0f}q/s\n", ref.size() / totalTime);
    fmt::print("  number of hits:      {:>10}\n", results.size());
    fmt::print("  number of kmerHits   {:>10}\n", totalKmerHits);
}
}
