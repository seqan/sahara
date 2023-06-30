#include "utils.h"

#include "cereal/types/unordered_map.hpp"
#include <clice/clice.h>
#include <fmindex-collection/fmindex-collection.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <string>

namespace {
void app();
auto cli = clice::Argument{ .arg    = "kmer-index",
                            .desc   = "construct an index over a given input file",
                            .value  = std::filesystem::path{},
                            .cb     = app,
};

auto cliKmer = clice::Argument{ .parent = &cli,
                                .arg    = "--kmer",
                                .desc   = "splitting the text into kmers",
                                .value  = size_t{1}
};

auto cliWindow = clice::Argument{ .parent = &cli,
                                .arg    = "--window",
                                .desc   = "using windows",
                                .value  = size_t{1}
};

void app() {
    using Alphabet = ivs::d_dna5;
    constexpr size_t Sigma = Alphabet::size();
    constexpr size_t KmerSigma = 256;

    fmt::print("constructing an index for {}\n", *cli);

    auto timing = std::vector<std::tuple<std::string, double>>{};
    auto stopWatch = StopWatch();

    // load fasta file
    auto reader = ivio::fasta::reader {{*cli}};
    size_t totalSize{};
    size_t kmerLen{};
    auto ref = std::vector<std::vector<uint8_t>>{};
    auto ref_kmer = std::vector<std::vector<uint32_t>>{};
    auto uniq = std::unordered_map<size_t, uint32_t>{};
    for (auto record : reader) {
        totalSize += record.seq.size();
        ref.emplace_back(ivs::convert_char_to_rank<ivs::d_dna5>(record.seq));
        auto seq = ivs::convert_char_to_rank<ivs::dna5>(record.seq);

        ref_kmer.emplace_back();
        for (auto v : ivs::winnowing_minimizer<ivs::dna5>(seq, /*.k=*/*cliKmer, /*.window=*/*cliWindow)) {
            if (auto iter = uniq.find(v); iter != uniq.end()) {
                ref_kmer.back().emplace_back(iter->second);
            } else {
                if (uniq.size() >= KmerSigma) throw std::runtime_error{fmt::format("to many different kmers {} >= KmerSigma, doesn't fit into OccTable", uniq.size(), KmerSigma)};
                uniq[v] = uniq.size()+1;
                ref_kmer.back().emplace_back(uniq[v]);
            }
        }
        kmerLen += ref_kmer.back().size();
    }


    fmt::print("config:\n");
    fmt::print("  file:            {}\n", *cli);
    fmt::print("  sigma:           {:>10}\n", Sigma);
    fmt::print("  references:      {:>10}\n", ref.size());
    fmt::print("  totalSize:       {:>10}\n", totalSize);
    fmt::print("  different kmers: {:>10}\n", uniq.size());
    fmt::print("  kmer-seq-len:    {:>10}\n", kmerLen);

    timing.emplace_back("ld queries", stopWatch.reset());

    // create index
    //using Table = fmindex_collection::occtable::interleaved32::OccTable<Sigma>;
    using Table = fmindex_collection::occtable::interleavedEPRV7::OccTable<Sigma>;
    auto index = fmindex_collection::BiFMIndex<Table, fmindex_collection::DenseCSA>{ref, /*samplingRate*/16, /*threadNbr*/1};

    timing.emplace_back("index creation", stopWatch.reset());

    // create kmer-index
    using KmerTable = fmindex_collection::occtable::interleavedEPRV7_32::OccTable<KmerSigma>;
    auto kmerIndex = fmindex_collection::BiFMIndex_32<KmerTable, fmindex_collection::CSA_32>{ref_kmer, /*samplingRate*/65536, /*threadNbr*/1};

    timing.emplace_back("kmer-index creation", stopWatch.reset());

    // save index
    auto indexPath = cli->string() + ".kmer.idx";
    auto ofs       = std::ofstream{indexPath, std::ios::binary};
    auto archive   = cereal::BinaryOutputArchive{ofs};
    archive(index, kmerIndex, *cliKmer, *cliWindow, uniq);
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
}
