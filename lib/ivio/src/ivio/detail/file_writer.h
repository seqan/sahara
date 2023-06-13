// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "concepts.h"

#include <filesystem>
#include <fstream>

namespace ivio {

class file_writer {
protected:
    std::ofstream ofs;

public:
    file_writer(std::filesystem::path path)
        : ofs(path, std::ios_base::out | std::ios_base::binary)
    {}

    file_writer() = delete;
    file_writer(file_writer const&) = delete;
    file_writer(file_writer&& _other) noexcept = default;
    ~file_writer() = default;

    auto operator=(file_writer const&) -> file_writer& = delete;
    auto operator=(file_writer&&) -> file_writer& = delete;

    auto write(std::span<char const> buffer) -> size_t {
        ofs.write(buffer.data(), buffer.size());
        return buffer.size();
    }

    void close() {
        ofs.close();
    }
};

static_assert(writer_c<file_writer>);
}
