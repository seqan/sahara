#pragma once

#include "utils/utils.h"

#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/DenseCSA.h>
#include <ivio/ivio.h>
#include <ivsigma/ivsigma.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <string>
#include <vector>

struct Query {
    std::string name;
    bool reverse;
    Query(std::string _name, bool _reverse)
        : name{std::move(_name)}
        , reverse{_reverse}
        {}
};
template <size_t Sigma>
auto loadQueries(std::filesystem::path path, bool reverse) {
    std::vector<std::vector<uint8_t>> queries;
    std::vector<Query> queryInfos;
    if (path.empty() || !std::filesystem::exists(path)) {
        return std::make_tuple(queries, queryInfos);
    }
    auto reader = ivio::fasta::reader {{path}};
    for (auto record : reader) {
        queries.emplace_back(ivs::convert_char_to_rank<ivs::d_dna5>(record.seq));
        if (reverse) {
            queries.emplace_back(ivs::reverse_complement_rank<ivs::d_dna5>(queries.back()));
        }
    }
    return std::make_tuple(queries, queryInfos);
}

template <typename CSA, typename Table>
auto loadIndex(std::string path, size_t samplingRate, size_t threadNbr) {
    auto sw = StopWatch{};
    auto indexPath = path + ".idx";
    if (!std::filesystem::exists(indexPath)) {
        auto [ref, refInfo] = loadQueries<Table::Sigma>(path, false);
        auto index = fmindex_collection::BiFMIndex<Table>{ref, samplingRate, threadNbr};
        // save index here
        auto ofs     = std::ofstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(index);
        return index;
    } else {
        auto ifs     = std::ifstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        auto index = fmindex_collection::BiFMIndex<Table>{fmindex_collection::cereal_tag{}};
        archive(index);
        std::cout << "loading took " << sw.peek() << "s\n";
        return index;
    }
}

template <typename CSA, typename Table>
auto loadDenseIndex(std::string path, size_t samplingRate, size_t threadNbr) {
    auto sw = StopWatch{};
    auto indexPath = path + ".idx";
    if (!std::filesystem::exists(indexPath)) {
        auto [ref, refInfo] = loadQueries<Table::Sigma>(path, false);
        auto index = fmindex_collection::BiFMIndex<Table, fmindex_collection::DenseCSA>{ref, samplingRate, threadNbr};
        // save index here
        auto ofs     = std::ofstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(index);
        return index;
    } else {
        auto ifs     = std::ifstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        auto index = fmindex_collection::BiFMIndex<Table, fmindex_collection::DenseCSA>{fmindex_collection::cereal_tag{}};
        archive(index);
        std::cout << "loading took " << sw.peek() << "s\n";
        return index;
    }
}
