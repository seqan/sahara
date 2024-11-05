// SPDX-FileCopyrightText: 2006-2024, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cereal/archives/binary.hpp>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <span>
#include <vector>

struct AdaptiveKmerIndex {
    // Different Kmer selection schemas
    enum class KmerMode : uint8_t {
        Winnowing = 0, // Using winnowing minimizers
        Mod = 1,       // Using mod mers
    };

    struct Config {
        KmerMode mode;
        size_t   kmerLen;                             // Length of a single kmer
        size_t   window;                              // Window size (only for mode == KmerMode::Winnowing)
        size_t   modExp;                              // Exponent for the shift value `2^modExp` (only for mode == KmerMode::Mod)
        size_t   largestValue;                        // Largest Value in the reference text
    };

private:
    struct Pimpl;
    std::unique_ptr<Pimpl> pimpl;

public:
    AdaptiveKmerIndex();
    AdaptiveKmerIndex(Config config, std::vector<std::vector<uint8_t>> _text);
    ~AdaptiveKmerIndex();

    auto config() const -> Config;
    void load(cereal::BinaryInputArchive& archive);
    void save(cereal::BinaryOutputArchive& archive) const;

    void search(std::span<uint8_t const> _query, std::function<void(size_t refid, size_t refpos)> const& _report) const;
};
