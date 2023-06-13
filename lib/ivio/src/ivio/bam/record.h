// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <cstddef>
#include <optional>
#include <ranges>
#include <span>
#include <string_view>
#include <vector>

namespace ivio::bam {

struct record_view {
    constexpr static auto rank_to_char = std::string_view{"=ACMGRSVTWYHKDBN"};
    constexpr static auto char_to_rank = []() {
        auto r = std::array<uint8_t, 256>{};
        r.fill(15);
        for (size_t i{0}; i < rank_to_char.size(); ++i) {
            r[rank_to_char[i]] = i;
        }
        return r;
    }();

    struct compact_seq {
        std::span<uint8_t const> data; // Data compressed as described in bam
        size_t                   size; // Number of elements stored in the sequence
        auto operator[](size_t i) const -> uint8_t {
            auto c = data[i/2];
            if (i % 2 == 0) {
                c >>= 4;
            }
            return c & 0xf;
        }

        struct iter {
            using value_type = uint8_t;
            compact_seq const* seq;
            size_t i{};
            auto operator<=>(iter const& _other) const = default;
            auto operator*() const {
                return seq->operator[](i);
            }
            auto operator++() -> iter& {
                ++i;
                return *this;
            }
        };
        friend auto begin(compact_seq const& seq) {
            return iter{&seq, 0};
        }
        friend auto end(compact_seq const& seq) {
            return iter{&seq, seq.size};
        }
    };
    int32_t                     refID;
    int32_t                     pos;
    uint8_t                     mapq;
    uint16_t                    bin;
    uint16_t                    flag;
    int32_t                     next_refID;
    int32_t                     next_pos;
    int32_t                     tlen;
    std::string_view            read_name;
    std::span<uint8_t const>    cigar;
    compact_seq                 seq;
    std::span<uint8_t const>    qual;
};

struct record {
    int32_t                     refID;
    int32_t                     pos;
    uint8_t                     mapq;
    uint16_t                    bin;
    uint16_t                    flag;
    int32_t                     next_refID;
    int32_t                     next_pos;
    int32_t                     tlen;
    std::string                 read_name;
    std::vector<uint8_t>        cigar;
    std::vector<uint8_t>        seq;
    size_t                      seq_len;
    std::vector<uint8_t>        qual;

    record() = default;
    record(record_view v)
        : refID      {v.refID}
        , pos        {v.pos}
        , mapq       {v.mapq}
        , bin        {v.bin}
        , flag       {v.flag}
        , next_refID {v.next_refID}
        , next_pos   {v.next_pos}
        , tlen       {v.tlen}
        , read_name  {v.read_name}
        , cigar      {begin(v.cigar), end(v.cigar)}
        , seq        {begin(v.seq.data), end(v.seq.data)}
        , seq_len    {v.seq.size}
        , qual       {begin(v.qual), end(v.qual)}
    {}
    operator record_view() const {
        return record_view {
            refID,
            pos,
            mapq,
            bin,
            flag,
            next_refID,
            next_pos,
            tlen,
            read_name,
            cigar,
            {seq, seq_len},
            qual,
        };
    }
};

}
