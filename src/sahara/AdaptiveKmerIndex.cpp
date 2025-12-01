// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include "AdaptiveKmerIndex.h"
#include "utils/error_fmt.h"

#include <cereal/types/array.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/locate.h>
#include <fmindex-collection/search/all.h>
#include <fmindex-collection/string/PairedFlattenedBitvectors2L.h>
#include <fmindex-collection/suffixarray/DenseCSA.h>
#include <fstream>
#include <unordered_map>

static constexpr size_t KmerSigma = 128;

struct AdaptiveKmerIndex::Pimpl {
    Config config;

    std::unordered_map<size_t, uint8_t> denseMap; // converting kmer hash values to uniq dense values

    // create kmer-index
    template <size_t Sigma>
    using Index = fmc::FMIndex<Sigma, fmc::string::PairedFlattenedBitvectors_512_64k>;
    std::variant<Index<3>, Index<4>, Index<5>, Index<6>, Index<16>, Index<32>, Index<64>, Index<128>> index{Index<3>{}};

    void initIndex() {
        size_t maxValue = config.largestValue;
        if (maxValue < 3) index = Index<3>{};
        else if (maxValue <   4) index = Index<  4>{};
        else if (maxValue <   5) index = Index<  5>{};
        else if (maxValue <   6) index = Index<  6>{};
        else if (maxValue <  16) index = Index< 16>{};
        else if (maxValue <  32) index = Index< 32>{};
        else if (maxValue <  64) index = Index< 64>{};
        else if (maxValue < 128) index = Index<128>{};
    }
};

AdaptiveKmerIndex::AdaptiveKmerIndex()
    : pimpl{std::make_unique<AdaptiveKmerIndex::Pimpl>()}
{}


AdaptiveKmerIndex::AdaptiveKmerIndex(Config config, std::vector<std::vector<uint8_t>> _text)
    : AdaptiveKmerIndex{} {
    if (config.largestValue > 128) {
        throw error_fmt("text with values above 128 is not allowed (requested largest value: {})", config.largestValue);
    }
    pimpl->config = config;
    pimpl->initIndex();
    std::visit([&]<typename Index>(Index& index) {
        index = Index{std::move(_text), /*sampingRate=*/16, /*threadNbr=*/1};
    }, pimpl->index);
}

AdaptiveKmerIndex::~AdaptiveKmerIndex() = default;
auto AdaptiveKmerIndex::config() const -> Config {
    return pimpl->config;
}

void AdaptiveKmerIndex::load(cereal::BinaryInputArchive& archive) {
    archive(pimpl->config.largestValue);
    pimpl->initIndex();
    std::visit([&](auto& index) {
        archive(index);
    }, pimpl->index);
    archive(pimpl->config.kmerLen, pimpl->config.mode);
    if (pimpl->config.mode == KmerMode::Winnowing) {
        archive(pimpl->config.window);
    } else if (pimpl->config.mode == KmerMode::Mod) {
        archive(pimpl->config.modExp);
    } else {
        throw error_fmt("unknown kmer mode {}", uint8_t(pimpl->config.mode));
    }

}

void AdaptiveKmerIndex::save(cereal::BinaryOutputArchive& archive) const {
    archive(pimpl->config.largestValue);
    std::visit([&](auto const& index) {
        archive(index);
    }, pimpl->index);
    archive(pimpl->config.kmerLen, pimpl->config.mode);
    if (pimpl->config.mode == KmerMode::Winnowing) {
        archive(pimpl->config.window);
    } else if (pimpl->config.mode == KmerMode::Mod) {
        archive(pimpl->config.modExp);
    } else {
        throw error_fmt("missing code path for unknown kmer mode type {}", uint8_t(pimpl->config.mode));
    }
}

void AdaptiveKmerIndex::search(std::span<uint8_t const> _query, std::function<void(size_t refid, size_t refpos)> const& _report) const {
     std::visit([&](auto const& index) {
        auto cursor = fmc::search_no_errors::search(index, _query);
        for (auto [seqId, seqPos, offset] : fmc::LocateLinear{index, cursor}) {
            _report(seqId, seqPos + offset);
        }
   }, pimpl->index);
}
