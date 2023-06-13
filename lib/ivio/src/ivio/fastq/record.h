// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <string>
#include <string_view>

namespace ivio::fastq {

struct record_view {
    std::string_view id;
    std::string_view seq;
    std::string_view id2;
    std::string_view qual;
};

struct record {
    std::string id;
    std::string seq;
    std::string id2;
    std::string qual;

    record() = default;
    record(record_view v)
        : id   {v.id}
        , seq  {v.seq}
        , id2  {v.id2}
        , qual {v.qual}
    {}

    operator record_view() const {
        return record_view {
            id,
            seq,
            id2,
            qual
        };
    }
};

}
