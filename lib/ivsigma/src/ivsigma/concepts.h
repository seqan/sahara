// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <cstdint>
#include <concepts>

namespace ivs {

template <typename Alphabet>
concept alphabet_c = requires(Alphabet) {
    { Alphabet::char_to_rank(char{}) } -> std::same_as<uint8_t>;
    { Alphabet::rank_to_char(uint8_t{}) } -> std::same_as<char>;
    { Alphabet::normalize_char(char{}) } -> std::same_as<char>;
};

template <typename Alphabet>
concept alphabet_with_complement_c = requires(Alphabet) {
    { Alphabet::char_to_rank(char{}) } -> std::same_as<uint8_t>;
    { Alphabet::rank_to_char(uint8_t{}) } -> std::same_as<char>;
    { Alphabet::normalize_char(char{}) } -> std::same_as<char>;
    { Alphabet::complement_rank(uint8_t{}) } -> std::same_as<uint8_t>;
    { Alphabet::complement_char(char{}) } -> std::same_as<char>;
};


}
