#include "error_fmt.h"
#include "utils/StopWatch.h"

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <clice/clice.h>
#include <fmindex-collection/DenseCSA.h>
#include <fmindex-collection/fmindex-collection.h>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>
#include <string>

namespace {
void app();
auto cli = clice::Argument{ .arg    = "index",
                            .desc   = "construct an index over a given input file",
                            .value  = std::filesystem::path{},
                            .cb     = app,
};

void app() {
    using Alphabet = ivs::d_dna5;
    constexpr size_t Sigma = Alphabet::size();

    fmt::print("constructing an index for {}\n", *cli);
    using Table = fmindex_collection::occtable::interleaved32::OccTable<Sigma>;

    auto timing = std::vector<std::tuple<std::string, double>>{};
    auto stopWatch = StopWatch();

    // load fasta file
    size_t totalSize{};
    auto ref = std::vector<std::vector<uint8_t>>{};
    for (auto record : ivio::fasta::reader {{*cli}}) {
        totalSize += record.seq.size();
        ref.emplace_back(ivs::convert_char_to_rank<Alphabet>(record.seq));
        if (!ivs::verify_rank(ref.back())) {
            throw error_fmt{"reference '{}' ({}) has invalid characters", record.id, ref.size()};
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
    auto ofs       = std::ofstream{indexPath, std::ios::binary};
    auto archive   = cereal::BinaryOutputArchive{ofs};
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
}
