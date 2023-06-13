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

    fmt::print("constructing an index over {}\n", *cli);
    using Table = fmindex_collection::occtable::interleaved32::OccTable<Sigma>;

    // load fasta file
    auto [ref, refInfo] = loadQueries<Table::Sigma>(*cli, false);

    // create index
    auto index = fmindex_collection::BiFMIndex<Table, fmindex_collection::DenseCSA>{ref, /*samplingRate*/16, /*threadNbr*/1};

    // save index
    auto indexPath = cli->string() + "." + Table::extension() + ".dense.index";
    auto ofs       = std::ofstream{indexPath, std::ios::binary};
    auto archive   = cereal::BinaryOutputArchive{ofs};
    archive(index);

}
}
