#include "utils.h"

#include <clice/clice.h>
#include <fmindex-collection/fmindex-collection.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <string>

namespace {
void app();
auto cli = clice::Argument{ .arg    = "index",
                            .desc   = "construct an index over a given input file",
                            .value  = std::filesystem::path{},
                            .cb     = app,
};

void app() {
    constexpr size_t Sigma = 5;


    fmt::print("constructing an index for {}\n", *cli);
    using Table = fmindex_collection::occtable::interleaved32::OccTable<Sigma>;

    auto timing = std::vector<std::tuple<std::string, double>>{};
    auto stopWatch = StopWatch();

    // load fasta file
    auto [ref, refInfo] = loadQueries<Table::Sigma>(*cli, /*reverse*/false);
    size_t totalSize{};
    for (auto const& r : ref) totalSize += r.size();

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
