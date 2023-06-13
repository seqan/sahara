// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>

namespace ivio::vcf {

struct record;

struct record_view {
    std::string_view            chrom;
    int32_t                     pos;
    std::string_view            id;
    std::string_view            ref;
    std::string_view            alts;
    std::optional<float>        qual;
    std::string_view            filters;
    std::string_view            infos;
    std::string_view            formats;
    std::string_view            samples;

    operator record() const;
    auto operator<=>(record_view const&) const = default;
};

struct record {
    std::string                 chrom;
    int32_t                     pos;
    std::string                 id;
    std::string                 ref;
    std::string                 alts;
    std::optional<float>        qual;
    std::string                 filters;
    std::string                 infos;
    std::string                 formats;
    std::string                 samples;

    operator record_view() const&;
    auto operator<=>(record const&) const = default;
};

// Implementation of the convert operators
inline record_view::operator record() const {
    return {
        .chrom    = std::string{chrom},
        .pos      = pos,
        .id       = std::string{id},
        .ref      = std::string{ref},
        .alts     = std::string{alts},
        .qual     = qual,
        .filters  = std::string{filters},
        .infos    = std::string{infos},
        .formats  = std::string{formats},
        .samples  = std::string{samples},
    };
}
inline record::operator record_view() const& {
    return {
        .chrom    = chrom,
        .pos      = pos,
        .id       = id,
        .ref      = ref,
        .alts     = alts,
        .qual     = qual,
        .filters  = filters,
        .infos    = infos,
        .formats  = formats,
        .samples  = samples,
    };
}

}

// Specialization to describe their common types
template <>
struct std::common_type<ivio::vcf::record, ivio::vcf::record_view> {
    using type = ivio::vcf::record;
};

template <>
struct std::common_type<ivio::vcf::record_view, ivio::vcf::record> {
    using type = ivio::vcf::record;
};
