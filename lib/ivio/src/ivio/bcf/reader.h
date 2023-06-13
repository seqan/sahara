// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../detail/reader_base.h"
#include "header.h"
#include "record.h"

#include <filesystem>
#include <optional>
#include <tuple>
#include <variant>
#include <vector>
#include <unordered_map>

namespace ivio::bcf {

struct reader : public reader_base<reader> {
    using record      = bcf::record;
    using record_view = bcf::record_view;

    bcf::header header_;

    struct config {
        // Source: file or stream
        std::variant<std::filesystem::path, std::reference_wrapper<std::istream>> input;

        // Value of 0 will run with a sequential implementation, other values will spawn new threads
        size_t threadNbr = 0;
    };

public:
    reader(config const& config_);
    ~reader();

    auto header() const -> bcf::header const& { return header_; }
    auto next() -> std::optional<record_view>;
    void close();
};

}
