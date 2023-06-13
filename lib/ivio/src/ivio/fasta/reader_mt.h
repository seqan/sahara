// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "record.h"

#include <filesystem>
#include <functional>
#include <istream>
#include <memory>
#include <optional>
#include <string_view>
#include <variant>

namespace ivio::fasta {

struct reader_mt_config {
    // Source: file or stream
    std::variant<std::filesystem::path, std::reference_wrapper<std::istream>> input;

    // This is only relevant if a stream is being used
    bool compressed{};
};

struct reader_mt_pimpl;
struct reader_mt {
    using record_view = fasta::record_view;
    struct iter;

private:
    std::unique_ptr<reader_mt_pimpl> pimpl;

public:
    std::vector<std::tuple<std::string, std::string>> header;
    std::vector<std::string> genotypes;

    reader_mt(reader_mt_config config);
    ~reader_mt();

    auto next() -> std::optional<record_view>;

    friend auto begin(reader_mt& reader) -> iter;
    friend auto end(reader_mt& reader) {
        return nullptr;
    }
};

struct reader_mt::iter {
    reader_mt& _reader;
    std::optional<record_view> nextItem = _reader.next();

    auto operator*() const -> record_view {
       return *nextItem;
    }
    auto operator++() -> iter& {
        nextItem = _reader.next();
        return *this;
    }
    auto operator!=(std::nullptr_t) const {
        return nextItem.has_value();
    }
};


}
