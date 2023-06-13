// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "buffered_reader.h"
#include "file_reader.h"
#include "mmap_reader.h"
#include "stream_reader.h"

#include <zlib.h>
#include <ranges>

namespace ivio {

struct zlib_reader {
    VarBufferedReader reader;

    z_stream stream {
        .next_in = Z_NULL,
        .avail_in = 0,
        .total_out = 0,
        .zalloc = Z_NULL,
        .zfree = Z_NULL,
        .opaque = Z_NULL,
    };

    zlib_reader(VarBufferedReader reader)
        : reader{std::move(reader)}
    {
        if (inflateInit2(&stream, 16 + MAX_WBITS) != Z_OK) {
            throw std::runtime_error{"error"};
        }
    }
    zlib_reader(zlib_reader&& _other)
        : reader{std::move(_other.reader)}
    {
        if (inflateInit2(&stream, 16 + MAX_WBITS) != Z_OK) {
            throw std::runtime_error{"error"};
        }
    }

    ~zlib_reader() {
        inflateEnd(&stream);
    }

    size_t read(std::ranges::sized_range auto&& range) {
        while(true) {
            auto [ptr, avail_in] = reader.read(range.size());

            stream.next_in  = (unsigned char*)(ptr);
            stream.avail_in = avail_in;
            stream.avail_out = range.size();
            stream.next_out  = (unsigned char*)&*std::begin(range);
            auto ret = inflate(&stream, Z_NO_FLUSH);
            auto diff = avail_in - stream.avail_in;

            reader.dropUntil(diff);

            auto producedBytes = (size_t)(stream.next_out - (unsigned char*)&*std::begin(range));
            if (producedBytes > 0) {
                return producedBytes;
            }
            if (ret == Z_STREAM_END) {
                return 0;
            }
        }
    }
};

static_assert(Readable<zlib_reader>);

}
