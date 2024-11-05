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
#include <string>
#include <unordered_set>

using namespace fmindex_collection;

namespace {
void app();
auto cli = clice::Argument{ .args   = "uni-search",
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


auto cliNoReverse = clice::Argument{ .parent = &cli,
                                     .args   = "--no-reverse",
                                     .desc   = "do not search for reversed complements",
};

void app() {
    using Alphabet = ivs::d_dna5;
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
    if (queries.empty()) {
        throw error_fmt{"query file {} was empty - abort\n", *cliQuery};
    }
    timing.emplace_back("ld queries", stopWatch.reset());

    fmt::print(
        "config:\n"
        "  query:               {}\n"
        "  index:               {}\n"
        "  reverse complements: {}\n"
        "  output path:         {}\n",
        *cliQuery, *cliIndex, !cliNoReverse, *cliOutput);


    {
        auto fwdQueries = queries.size() / (cliNoReverse?1:2);
        auto bwdQueries = queries.size() - fwdQueries;
        fmt::print("fwd queries: {}\n"
                   "bwd queries: {}\n",
                   fwdQueries, bwdQueries);
    }

    using Table = fmindex_collection::occtable::Interleaved_32<Sigma>;
//    using Table = fmindex_collection::occtable::EprV7<Sigma>;

    if (!std::filesystem::exists(*cliIndex)) {
        throw error_fmt{"no valid index path at {}", *cliIndex};
    }

    auto index = fmindex_collection::FMIndex<Table, fmindex_collection::DenseCSA>{};
    {
        auto ifs     = std::ifstream{*cliIndex, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        archive(index);
    }
    timing.emplace_back("ld index", stopWatch.reset());


    auto resultCursors = std::vector<std::tuple<size_t, FMIndexCursor<decltype(index)>>>{};
    for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
        auto const& query = queries[qidx];
        auto cursor = fmindex_collection::search_no_errors::search(index, query);
        resultCursors.emplace_back(qidx, cursor);
    }

    timing.emplace_back("search", stopWatch.reset());

    auto results = std::vector<std::tuple<size_t, size_t, size_t>>{};
    for (auto const& [queryId, cursor] : resultCursors) {
        for (auto [seqId, pos] : LocateLinear{index, cursor}) {
            results.emplace_back(queryId, seqId, pos);
        }
    }

    timing.emplace_back("locate", stopWatch.reset());

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
    fmt::print("  queries per second:  {:> 10.0f}q/s\n", queries.size() / totalTime);
    fmt::print("  number of hits:      {:>10}\n", results.size());
}
}
