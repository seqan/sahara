// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../detail/reader_base.h"
#include "../faidx/record.h"
#include "record.h"

#include <filesystem>
#include <optional>
#include <tuple>
#include <variant>

namespace ivio::fasta {

struct reader : public reader_base<reader> {
    using record      = fasta::record;
    using record_view = fasta::record_view;

    struct config {
        // Source: file or stream
        std::variant<std::filesystem::path, std::reference_wrapper<std::istream>> input;

        // This is only relevant if a stream is being used
        bool compressed{};
    };

public:
    reader(config const& config_);
    ~reader();

    auto next() -> std::optional<record_view>;
    void close();

    auto tell() const -> size_t;
    void seek(size_t offset);

    auto tell_faidx() const -> faidx::record;
    void seek_faidx(faidx::record const& offset);
};

}
