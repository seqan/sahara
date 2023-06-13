// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "simple_alphabet.h"

#include <cstddef>

namespace ivs {

template <size_t Size, char Offset>
struct ranged_alphabet {
    /*! \brief Conerts a single char value to a rank value
     *
     * \param c char value
     * \return corresponding rank value
     */
    static constexpr auto char_to_rank(char c) noexcept -> uint8_t {
        auto index = static_cast<uint8_t>(c);
        if (c < Offset or c-Offset > Size) {
            return 255;
        }
        return c - Offset;
    }

    /*! \brief Converts a single rank value to a char value
     *
     * \param v a rank value
     * \return corresponding (normalized) char value
     */
    static constexpr auto rank_to_char(uint8_t v) noexcept -> char {
        if (v > Size) return '\0';
        return v + Offset;
    }
};


using pthred42       = ranged_alphabet<42, '!'>;
using pthred63       = ranged_alphabet<63, '!'>;
using pthred68solexa = ranged_alphabet<68, ';'>;
using pthred94       = ranged_alphabet<94, '!'>;

}
