// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "../detail/writer_base.h"
#include "record.h"

#include <filesystem>
#include <functional>
#include <ostream>
#include <variant>

namespace ivio::fasta {

struct writer : writer_base<writer> {
    struct config {
        // Source: file or stream
        std::variant<std::filesystem::path, std::reference_wrapper<std::ostream>> output;

        // This is only relevant if a stream is being used
        bool compressed{};

        size_t length{80}; // Break after 80 characters
    };

    writer(config config);
    ~writer();

    void write(record_view record);
    void close();
};

}
