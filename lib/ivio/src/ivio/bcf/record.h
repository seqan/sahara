// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace ivio::bcf {

struct record_view {
    int32_t                  chromId;
    int32_t                  pos;
    int32_t                  rlen;
    std::optional<float>     qual;
    uint16_t                 n_info;
    uint16_t                 n_allele;
    uint32_t                 n_sample; // Only the first 24bit can be used
    uint8_t                  n_fmt;
    std::string_view         id;
    std::string_view         ref;
    std::span<uint8_t const> alt;
    std::span<uint8_t const> filter;
    std::span<uint8_t const> info;
    std::span<uint8_t const> format;
};

struct record {
    int32_t                  chromId;
    int32_t                  pos;
    int32_t                  rlen;
    std::optional<float>     qual;
    uint16_t                 n_info;
    uint16_t                 n_allele;
    uint32_t                 n_sample; // Only the first 24bit can be used
    uint8_t                  n_fmt;
    std::string              id;
    std::string              ref;
    std::vector<uint8_t>     alt;
    std::vector<uint8_t>     filter;
    std::vector<uint8_t>     info;
    std::vector<uint8_t>     format;

    record() = default;
    record(record_view v)
        : chromId   {v.chromId}
        , pos       {v.pos}
        , rlen      {v.rlen}
        , qual      {v.qual}
        , n_info    {v.n_info}
        , n_allele  {v.n_allele}
        , n_sample  {v.n_sample}
        , n_fmt     {v.n_fmt}
        , id        {v.id}
        , ref       {v.ref}
        , alt       {begin(v.alt), end(v.alt)}
        , filter    {begin(v.filter), end(v.filter)}
        , info      {begin(v.info), end(v.info)}
        , format    {begin(v.format), end(v.format)}
    {}
    operator record_view() const {
        return record_view {
            chromId,
            pos,
            rlen,
            qual,
            n_info,
            n_allele,
            n_sample,
            n_fmt,
            id,
            ref,
            alt,
            filter,
            info,
            format,
        };
    }
};


}
