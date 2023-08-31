#include "AdaptiveKmerIndex.h"
#include "hash.h"
#include "utils/StopWatch.h"
#include "utils/error_fmt.h"

#include <cereal/types/unordered_map.hpp>
#include <clice/clice.h>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <string>

namespace {
void app();
auto cli = clice::Argument{ .args   = "kmer-index",
                            .desc   = "construct an index over a given input file",
                            .value  = std::filesystem::path{},
                            .cb     = app,
};

auto cliKmer = clice::Argument{ .parent = &cli,
                                .args   = "--kmer",
                                .desc   = "splitting the text into kmers",
                                .value  = size_t{1}
};

auto cliKmerMode = clice::Argument{ .parent = &cli,
                                    .args   = "--kmer_mode",
                                    .desc   = "valid modes are: winnowing and mod",
                                    .value  = AdaptiveKmerIndex::KmerMode::Winnowing,
                                    .mapping = {{{"winnowing", AdaptiveKmerIndex::KmerMode::Winnowing}, {"mod", AdaptiveKmerIndex::KmerMode::Mod}}},
};

auto cliWindow = clice::Argument{ .parent = &cli,
                                .args   = "--window",
                                .desc   = "using windows (only valid for '--kmer_mode winnowing' mode",
                                .value  = size_t{1}
};

auto cliMod = clice::Argument{ .parent = &cli,
                               .args   = "--mod",
                               .desc   = "take every 'mod' element (only valid for '--kmer_mode mod' mode",
                               .value  = size_t{4}
};

auto cliIgnoreUnknown = clice::Argument{ .parent = &cli,
                                         .args   = "--ignore_unknown",
                                         .desc   = "ignores unknown nuclioteds in input data and replaces them with 'N'",
};

void app() {
    using Alphabet = ivs::d_dna5;
    constexpr size_t Sigma = Alphabet::size();
    constexpr size_t KmerSigma = 128;

    fmt::print("constructing an index for {}\n", *cli);

    auto timing = std::vector<std::tuple<std::string, double>>{};
    auto stopWatch = StopWatch();

    // load fasta file
    auto reader = ivio::fasta::reader {{*cli}};
    size_t totalSize{};
    size_t kmerLen{};
    auto ref_kmer = std::vector<std::vector<uint8_t>>{};
    auto uniq = std::unordered_map<size_t, uint8_t>{};
    auto ref = std::vector<uint8_t>{};
    size_t recordNbr = 0;
    for (auto record : reader) {
        recordNbr += 1;
        totalSize += record.seq.size();
        ref.resize(record.seq.size());
        ivs::convert_char_to_rank<Alphabet>(record.seq, ref);
        if (auto pos = ivs::verify_rank(ref); pos) {
            if (cliIgnoreUnknown) {
                ref[*pos] = Alphabet::char_to_rank('N');
            } else {
                throw error_fmt{"ref '{}' ({}) has invalid character at position {} '{}'({:x})", record.id, recordNbr, *pos, record.seq[*pos], record.seq[*pos]};
            }
        }


        ref_kmer.emplace_back();
        if (*cliKmerMode == AdaptiveKmerIndex::KmerMode::Winnowing) {
            for (auto v : ivs::winnowing_minimizer<Alphabet, /*DuplicatesAllowed=*/false>(ref, /*.k=*/*cliKmer, /*.window=*/*cliWindow)) {
                if (auto iter = uniq.find(v); iter != uniq.end()) {
                    ref_kmer.back().emplace_back(iter->second);
                } else {
                    uniq[v] = uniq.size()+1;
                    ref_kmer.back().emplace_back(uniq[v]);
                }
            }
        } else if (*cliKmerMode == AdaptiveKmerIndex::KmerMode::Mod) {
            uint64_t mask = (1<<*cliMod)-1;
            for (auto v : ivs::compact_encoding<Alphabet>(ref, /*.k=*/*cliKmer)) {
                v = hash(v);
                if ((v & mask) != 0) continue;
                if (auto iter = uniq.find(v); iter != uniq.end()) {
                    ref_kmer.back().emplace_back(iter->second);
                } else {
                    uniq[v] = uniq.size()+1;
                    ref_kmer.back().emplace_back(uniq[v]);
                }
            }
        } else {
            throw error_fmt("unknown kmer mode: {}", uint8_t(*cliKmerMode));
        }
        kmerLen += ref_kmer.back().size();
    }
    if (uniq.size() >= KmerSigma) throw error_fmt{"to many different kmers {} >= {}, doesn't fit into OccTable", uniq.size(), KmerSigma};


    fmt::print("config:\n");
    fmt::print("  file:            {}\n", *cli);
    fmt::print("  sigma:           {:>10}\n", Sigma);
    fmt::print("  references:      {:>10}\n", ref_kmer.size());
    fmt::print("  totalSize:       {:>10}\n", totalSize);
    if (*cliKmerMode == AdaptiveKmerIndex::KmerMode::Winnowing) {
        fmt::print("  kmerMode:        {:>10}\n", "winnowing");
        fmt::print("  windowSize       {:>10}\n", *cliWindow);
    } else if (*cliKmerMode == AdaptiveKmerIndex::KmerMode::Mod) {
        fmt::print("  kmerMode:        {:>10}\n", "mod");
        fmt::print("  modFactor        {:>10}\n", fmt::format("2^{}", *cliMod));
    } else {
        throw std::runtime_error("missing code path for unknown kmer mode type");
    }
    fmt::print("  different kmers: {:>10}\n", uniq.size());
    fmt::print("  kmer-seq-len:    {:>10}\n", kmerLen);

    timing.emplace_back("ld queries", stopWatch.reset());

    // create kmer-index
    auto index = AdaptiveKmerIndex {{
        .mode         = *cliKmerMode,
        .kmerLen      = *cliKmer,
        .window       = *cliWindow,
        .modExp       = *cliMod,
        .largestValue = uniq.size()
        }, std::move(ref_kmer)};


    timing.emplace_back("index creation", stopWatch.reset());

    // save index
    auto indexPath = cli->string() + ".kmer.idx";
    auto ofs       = std::ofstream{indexPath, std::ios::binary};
    auto archive   = cereal::BinaryOutputArchive{ofs};
    auto fileFormatVersion = uint32_t{0x01}; // Saving as format v0x01
    archive(fileFormatVersion);
    index.save(archive);
    archive(uniq);
    ofs.close();

    timing.emplace_back("saving to disk", stopWatch.reset());

    fmt::print("stats:\n");
    double totalTime{};
    for (auto const& [key, time] : timing) {
        fmt::print("  {:<25} {:> 10.2f}s\n", key + " time:", time);
        totalTime += time;
    }
    fmt::print("  total time:               {:> 10.2f}s\n", totalTime);

}
}
