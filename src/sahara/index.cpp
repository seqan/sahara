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
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <string>

namespace {
void app();
auto cli = clice::Argument{ .args   = "index",
                            .desc   = "construct an index over a given input file",
                            .value  = std::filesystem::path{},
                            .cb     = app,
};

auto cliIgnoreUnknown = clice::Argument{ .parent = &cli,
                                         .args   = "--ignore_unknown",
                                         .desc   = "ignores unknown nuclioteds in input data and replaces them with 'N'",
};

auto cliUseDna4 = clice::Argument{ .parent = &cli,
                                   .args   = "--dna4",
                                   .desc   = "use dna 4 alphabet, replace 'N' with random ACG or T",
};



template <typename Alphabet>
void createIndex() {
    constexpr size_t Sigma = Alphabet::size();

    fmt::print("constructing an index for {}\n", *cli);
    using Table = fmindex_collection::occtable::Interleaved_32<Sigma>;

    auto timing = std::vector<std::tuple<std::string, double>>{};
    auto stopWatch = StopWatch();

    // load fasta file
    size_t totalSize{};
    auto ref = std::vector<std::vector<uint8_t>>{};
    for (auto record : ivio::fasta::reader {{*cli}}) {
        totalSize += record.seq.size();
        ref.emplace_back(ivs::convert_char_to_rank<Alphabet>(record.seq));
        if (cliIgnoreUnknown) {
            if (cliUseDna4) {
                for (auto& v : ref.back()) {
                    if (ivs::verify_rank(v)) continue;
                    v = Alphabet::char_to_rank('A') + rand() % 4;
                }
            } else {
                for (auto& v : ref.back()) {
                    if (ivs::verify_rank(v)) continue;
                    v = Alphabet::char_to_rank('N');
                }
            }
        }
        if (auto pos = ivs::verify_rank(ref.back()); pos) {
            throw error_fmt{"ref '{}' ({}) has invalid character '{}' (0x{:02x}) at position {}", record.id, ref.size(), record.seq[*pos], record.seq[*pos], *pos};
        }
    }
    if (ref.empty()) {
        throw error_fmt{"reference file {} was empty - abort\n", *cli};
    }


    fmt::print("config:\n");
    fmt::print("  file: {}\n", *cli);
    fmt::print("  sigma: {}\n", Sigma);
    fmt::print("  references: {}\n", ref.size());
    fmt::print("  totalSize: {}\n", totalSize);

    timing.emplace_back("ld queries", stopWatch.reset());

    // create index
    auto index = fmindex_collection::BiFMIndex<Table, fmindex_collection::DenseCSA>{ref, /*samplingRate*/16, /*threadNbr*/1};

    timing.emplace_back("index creation", stopWatch.reset());

    // save index
    auto indexPath = cli->string() + ".idx";
    if (cliUseDna4) {
        indexPath = cli->string() + ".dna4.idx";
    }
    auto ofs       = std::ofstream{indexPath, std::ios::binary};
    auto archive   = cereal::BinaryOutputArchive{ofs};
    archive(Sigma);
    archive(index);
    ofs.close();

    timing.emplace_back("saving to disk", stopWatch.reset());

    fmt::print("stats:\n");
    double totalTime{};
    for (auto const& [key, time] : timing) {
        fmt::print("  {:<20} {:> 10.2f}s\n", key + " time:", time);
        totalTime += time;
    }
    fmt::print("  total time:          {:> 10.2f}s\n", totalTime);

}


void app() {
    if (cliUseDna4) {
        createIndex<ivs::d_dna4>();
    } else {
        createIndex<ivs::d_dna5>();
    }
}
}
