// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "utils/StopWatch.h"
#include "utils/error_fmt.h"
#include "VarIndex.h"

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <clice/clice.h>
#include <fmindex-collection/suffixarray/DenseCSA.h>
#include <fmindex-collection/fmindex-collection.h>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <mmser/mmser.h>
#include <string>

namespace {
void app();
auto cli = clice::Argument {
    .args   = "index",
    .desc   = "construct an index over a given input file",
    .value  = std::filesystem::path{},
    .cb     = app,
};

auto cliIgnoreUnknown = clice::Argument {
    .parent = &cli,
    .args   = "--ignore_unknown",
    .desc   = "ignores unknown nuclioteds in input data and replaces them with 'N'",
};

auto cliIndexType = clice::Argument {
    .parent  = &cli,
    .args    = "--index_type",
    .desc    = "type of the index (implementation detail)",
    .value   = std::string{"ibv16"},
    .mapping = {{
        {"ibv16", "ibv16"},
        {"mbv64_64", "mbv64_64"},
        {"mbv512_64", "mbv512_64"},
        {"fbv64_64", "fbv64_64"},
        {"fbv512_64", "fbv512_64"},
    }},
};

auto cliIndexTypePaired = clice::Argument {
    .parent = &cliIndexType,
    .args   = "--paired",
    .desc   = "some types like fbv*_* have a specialzed 'paired' variant",
};

auto cliIndexTypeKStep = clice::Argument {
    .parent = &cliIndexType,
    .args   = "--k-step",
    .desc   = "enable additional k-step functionality, steps of 1 turns this function off",
    .value  = size_t{1},
};

auto cliIndexNoDelim = clice::Argument {
    .parent = &cliIndexType,
    .args   = "--no-delim",
    .desc   = "index type can also be build without delimiter, this introduces a false positives but decreases the alphabet size",
};

auto cliUseDna4 = clice::Argument {
    .parent = &cli,
    .args   = "--dna4",
    .desc   = "use dna 4 alphabet, replace 'N' with random ACG or T",
};

auto cliUseDna2 = clice::Argument {
    .parent = &cli,
    .args   = "--dna2",
    .desc   = "use dna 2 alphabet, replace 'N' with random ACG or T and reduce AT->S and CG->W",
};

auto cliIncludeReverse = clice::Argument {
    .parent = &cli,
    .args   = "--include-reverse",
    .desc   = "Includes the reverse text to the index",
};

auto cliThreads = clice::Argument {
    .parent = &cli,
    .args   = {"-t", "--threads"},
    .desc   = "number of threads to build the index",
    .value  = size_t{1},
};
auto cliSamplingRate = clice::Argument {
    .parent = &cli,
    .args   = {"-s", "--sampling_rate"},
    .desc   = "sampling rate of the fm index",
    .value  = size_t{16},
};
auto cliOutputFormat = clice::Argument {
    .parent  = &cli,
    .args    = {"--of", "--output_format"},
    .desc    = "cerealization technique used cereal or mmser",
    .value   = std::string{"mmser"},
    .mapping = {{
        {"cereal", "cereal"},
        {"mmser", "mmser"},
    }},
};

template <typename Alphabet>
void createIndex() {
    constexpr size_t Sigma = Alphabet::size();

    fmt::print("constructing an index for {}\n", *cli);

    auto timing = std::vector<std::tuple<std::string, double>>{};
    auto stopWatch = StopWatch();

    // load fasta file
    size_t totalSize{};
    auto ref = std::vector<std::vector<uint8_t>>{};
    for (auto record : ivio::fasta::reader {{*cli}}) {
        totalSize += record.seq.size();
        ref.emplace_back(ivs::convert_char_to_rank<Alphabet>(record.seq));
        if (cliIgnoreUnknown) {
            if (cliUseDna2) {
                for (auto& v : ref.back()) {
                    if (ivs::verify_rank(v)) continue;
                    v = Alphabet::char_to_rank('S') + rand() % 2;
                }
            } else if (cliUseDna4) {
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
    fmt::print("  threads: {}\n", *cliThreads);
    fmt::print("  samplingRate: {}\n", *cliSamplingRate);
    fmt::print("  index type: {}\n", *cliIndexType);
    fmt::print("    paired: {}\n", static_cast<bool>(cliIndexTypePaired));
    fmt::print("    use delimiter: {}\n", static_cast<bool>(cliIndexNoDelim));
    fmt::print("    k-step: {}\n", *cliIndexTypeKStep);
    fmt::print("  include reverse text: {}\n", static_cast<bool>(cliIncludeReverse));

    timing.emplace_back("ld queries", stopWatch.reset());

    // create index
    auto index = VarIndex<Alphabet>{};
    auto indexType = *cliIndexType;
    if (cliIndexTypeKStep) {
        indexType = fmt::format("{}_{}step", indexType, *cliIndexTypeKStep);
    }
    if (cliIndexNoDelim) {
        indexType += "-nd";
    }
    if (cliIncludeReverse) {
        indexType += "-rev";
    }
    if (cliIndexTypePaired) {
        indexType = "p" + indexType;
    }
    index.emplace(indexType, ref, *cliSamplingRate, *cliThreads);
    index.samplingRate = *cliSamplingRate;

    timing.emplace_back("index creation", stopWatch.reset());

    // save index
    auto indexPath = fmt::format("{}.{}.{}.idx", cli->string(), indexType, Sigma);
    fmt::print("  output path: {}\n", indexPath);

    if (*cliOutputFormat == "mmser") {
        mmser::saveFile(indexPath + ".mmser", index);
        timing.emplace_back("saving to disk via mmser", stopWatch.reset());
    } else if (*cliOutputFormat == "cereal") {
        auto ofs       = std::ofstream{indexPath, std::ios::binary};
        auto archive   = cereal::BinaryOutputArchive{ofs};
        archive(index);
        ofs.close();
        timing.emplace_back("saving to disk via cereal", stopWatch.reset());
    }

    fmt::print("stats:\n");
    double totalTime{};
    for (auto const& [key, time] : timing) {
        fmt::print("  {:<20} {:> 10.2f}s\n", key + " time:", time);
        totalTime += time;
    }
    fmt::print("  total time:          {:> 10.2f}s\n", totalTime);
}

void app() {
/*    if      (cliUseDna2 && cliIndexNoDelim)  createIndex<ivs::dna2>();
    else if (cliUseDna2 && !cliIndexNoDelim) createIndex<ivs::d_dna2>();
    else*/ if (cliUseDna4 && cliIndexNoDelim)  createIndex<ivs::dna4>();
/*    else if (cliUseDna4 && !cliIndexNoDelim) createIndex<ivs::d_dna4>();
    else if (cliIndexNoDelim)                createIndex<ivs::dna5>();
    else                                     createIndex<ivs::d_dna5>();*/
}
}
