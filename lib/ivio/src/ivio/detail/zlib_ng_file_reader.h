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

#include "zlib-ng.h"
#include <ranges>

namespace ivio {

template <BufferedReadable Reader>
struct zlib_ng_reader_impl {
    Reader file;

    zng_stream stream {
        .next_in = Z_NULL,
        .avail_in = 0,
        .total_out = 0,
        .zalloc = Z_NULL,
        .zfree = Z_NULL,
        .opaque = Z_NULL,
    };

    template <typename T>
    zlib_ng_reader_impl(T&& name)
        : file(std::forward<T>(name))
    {
        if (zng_inflateInit2(&stream, 16 + MAX_WBITS) != Z_OK) {
            throw "error";
        }
    }
    zlib_ng_reader_impl(zlib_ng_reader_impl&& _other)
        : file{std::move(_other.file)}
    {
        if (zng_inflateInit2(&stream, 16 + MAX_WBITS) != Z_OK) {
            throw "error";
        }
    }

    ~zlib_ng_reader_impl() {
        zng_inflateEnd(&stream);
    }

    size_t read(std::ranges::sized_range auto&& range) {
        while(true) {
            auto [ptr, avail_in] = file.read(range.size());

            stream.next_in  = (unsigned char*)(ptr);
            stream.avail_in = avail_in;
            stream.avail_out = range.size();
            stream.next_out  = (unsigned char*)&*std::begin(range);
            auto ret = zng_inflate(&stream, Z_NO_FLUSH);
            auto diff = avail_in - stream.avail_in;

            file.dropUntil(diff);

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

using zlib_ng_file_reader   = zlib_ng_reader_impl<buffered_reader<>>;
using zlib_ng_mmap_reader   = zlib_ng_reader_impl<mmap_reader>;
using zlib_ng_stream_reader = zlib_ng_reader_impl<buffered_reader<>>;

static_assert(Readable<zlib_ng_file_reader>);
static_assert(Readable<zlib_ng_mmap_reader>);
static_assert(Readable<zlib_ng_stream_reader>);


}
