// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../detail/reader_base.h"
#include "record.h"

#include <filesystem>
#include <optional>
#include <tuple>
#include <variant>
#include <vector>

namespace ivio::sam {

struct reader : public reader_base<reader> {
    using record      = sam::record;
    using record_view = sam::record_view;

    struct config {
        // Source: file or stream
        std::variant<std::filesystem::path, std::reference_wrapper<std::istream>> input;
    };

public:
    reader(config const& config_);
    ~reader();

    std::vector<std::string> header;

    auto next() -> std::optional<record_view>;
    void close();
};

}
