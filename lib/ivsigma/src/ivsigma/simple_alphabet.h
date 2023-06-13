// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "concepts.h"

#include <array>
#include <cassert>
#include <cstdint>
#include <tuple>

namespace ivs {

template <uint8_t V, char H, char... Tail>
struct rank_char_mapping {
    static constexpr uint8_t rank = V;
    static constexpr char symb = H;
    static constexpr std::array<char, sizeof...(Tail)+1> all_symb = {H, Tail...};
};



template <typename... values>
struct simple_alphabet {
private:
    //! Table representing conversion from rank to char
    static constexpr std::array<char, 256> rank_to_char_table{[]() {
        auto table = std::array<char, 256>{};
        table.fill('\0');
        auto add_values_to_table = [&]<typename V>(V) {
            table[V::rank] = V::symb;
        };
        (add_values_to_table(values{}), ...);
        return table;
    }()};

    //! Table representing conversion from char to rank
    static constexpr std::array<uint8_t, 256> char_to_rank_table{[]() {
        auto table = std::array<uint8_t, 256>{};
        table.fill(255);
        auto add_values_to_table = [&]<typename V>(V) {
            for (auto c : V::all_symb) {
                table[c] = V::rank;
            }
        };
        (add_values_to_table(values{}), ...);
        return table;
    }()};

    //! Table normalizing a char
    static constexpr std::array<char, 256> normalize_char_table{[]() {
        auto table = std::array<char, 256>{};
        table.fill('\0');
        auto add_values_to_table = [&]<typename V>(V) {
            for (auto c : V::all_symb) {
                table[c] = V::symb;
            }
        };
        (add_values_to_table(values{}), ...);
        return table;
    }()};

    //! Table computing the complement in rank space (rank in, rank out)
    static constexpr std::array<uint8_t, 256> rank_complement_table{[]() {
        auto table = std::array<uint8_t, 256>{};
        table.fill(255);
        auto add_values_to_table = [&]<typename V>(V) {
            table[V::rank] = char_to_rank_table[V::complement];
        };
        (add_values_to_table(values{}), ...);
        return table;
    }()};

public:
    /*! \brief Conerts a single char value to a rank value
     *
     * \param c char value
     * \return corresponding rank value
     */
    static constexpr auto char_to_rank(char c) noexcept -> uint8_t {
        auto index = static_cast<uint8_t>(c);
        return char_to_rank_table[index];
    }

    /*! \brief Converts a single rank value to a char value
     *
     * \param v a rank value
     * \return corresponding (normalized) char value
     */
    static constexpr auto rank_to_char(uint8_t v) noexcept -> char {
        return rank_to_char_table[v];
    }

    /*! \brief Normalizes a single char
     *
     * \param c char value
     * \return corresponding normalized char value
     */
    static constexpr auto normalize_char(char c) noexcept -> uint8_t {
        auto index = static_cast<uint8_t>(c);
        return normalize_char_table[index];
    }
};


template <uint8_t V, char H, char Compl, char... Tail>
struct rank_char_mapping_with_compl {
    static constexpr uint8_t rank = V;
    static constexpr char symb = H;
    static constexpr char complement = Compl;
    static constexpr std::array<char, sizeof...(Tail)+1> all_symb = {H, Tail...};
};

template <typename... values>
struct alphabet_with_compl {
private:
    //! Table representing conversion from rank to char
    static constexpr std::array<char, 256> rank_to_char_table{[]() {
        auto table = std::array<char, 256>{};
        table.fill('\0');
        auto add_values_to_table = [&]<typename V>(V) {
            table[V::rank] = V::symb;
        };
        (add_values_to_table(values{}), ...);
        return table;
    }()};

    //! Table representing conversion from char to rank
    static constexpr std::array<uint8_t, 256> char_to_rank_table{[]() {
        auto table = std::array<uint8_t, 256>{};
        table.fill(255);
        auto add_values_to_table = [&]<typename V>(V) {
            for (auto c : V::all_symb) {
                table[c] = V::rank;
            }
        };
        (add_values_to_table(values{}), ...);
        return table;
    }()};

    //! Table normalizing a char
    static constexpr std::array<char, 256> normalize_char_table{[]() {
        auto table = std::array<char, 256>{};
        table.fill('\0');
        auto add_values_to_table = [&]<typename V>(V) {
            for (auto c : V::all_symb) {
                table[c] = V::symb;
            }
        };
        (add_values_to_table(values{}), ...);
        return table;
    }()};

    //! Table computing the complement in rank space (rank in, rank out)
    static constexpr std::array<uint8_t, 256> rank_complement_table{[]() {
        auto table = std::array<uint8_t, 256>{};
        table.fill(255);
        auto add_values_to_table = [&]<typename V>(V) {
            table[V::rank] = char_to_rank_table[V::complement];
        };
        (add_values_to_table(values{}), ...);
        return table;
    }()};

public:
    /*! \brief Normalizes a single char
     *
     * \param c char value
     * \return corresponding normalized char value
     */
    static constexpr auto normalize_char(char c) noexcept -> uint8_t {
        auto index = static_cast<uint8_t>(c);
        return normalize_char_table[index];
    }

    /*! \brief Converts a single char value to a rank value
     *
     * \param c char value
     * \return corresponding rank value
     */
    static constexpr auto char_to_rank(char c) noexcept -> uint8_t {
        auto index = static_cast<uint8_t>(c);
        return char_to_rank_table[index];
    }

    /*! \brief Converts a single rank value to a char value
     *
     * \param v a rank value
     * \return corresponding (normalized) char value
     */
    static constexpr auto rank_to_char(uint8_t v) noexcept -> char {
        return rank_to_char_table[v];
    }

    /*! \brief computes complement based on ranks
     *
     * \param r a rank value
     * \return corresponding complement
     */
    static constexpr auto complement_rank(uint8_t r) noexcept -> uint8_t {
        return rank_complement_table[r];
    }

    /*! \brief computes complement based on char
     *
     * \param c a char value
     * \return corresponding complement
     */
    static constexpr auto complement_char(char c) noexcept -> char {
        return rank_to_char(complement_rank(char_to_rank(c)));
    }
};

}
