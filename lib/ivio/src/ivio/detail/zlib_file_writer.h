// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "buffered_writer.h"
#include "file_writer.h"
#include "stream_writer.h"

#include <zlib.h>
#include <ranges>

namespace ivio {

template <writer_c writer>
struct zlib_writer_impl {
    writer file;

    z_stream stream {
        .next_in = Z_NULL,
        .avail_in = 0,
        .total_out = 0,
        .zalloc = Z_NULL,
        .zfree = Z_NULL,
        .opaque = Z_NULL,
    };

    zlib_writer_impl(writer&& name)
        : file{std::move(name)}
    {
        if (deflateInit2(&stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
            throw std::runtime_error{"error initializing zlib/deflateInit2"};
        }
    }
    zlib_writer_impl(zlib_writer_impl&& _other)
        : file{std::move(_other.file)}
    {
        if (deflateInit2(&stream, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
            throw std::runtime_error{"error initializing zlib/deflateInit2"};
        }
    }

    ~zlib_writer_impl() {
        close();
    }

    auto write(std::span<char> buffer) -> size_t {
        auto outBuffer = std::array<char, 2<<16>{};

        stream.next_in = (unsigned char*)buffer.data();
        stream.avail_in = buffer.size();
        stream.next_out  = (unsigned char*)&outBuffer[0];
        stream.avail_out = outBuffer.size();

        auto ret = deflate(&stream, Z_NO_FLUSH);
        if (ret != Z_OK && ret != Z_STREAM_END) {
            throw std::runtime_error{"error deflating data with zlib"};
        }
        file.write({&outBuffer[0], outBuffer.size() - stream.avail_out});
        return buffer.size() - stream.avail_in;
    }

    void close() {
        if (stream.next_in) {
            auto outBuffer = std::array<char, 2<<16>{};

            do {
                stream.next_out  = (unsigned char*)&outBuffer[0];
                stream.avail_out = outBuffer.size();

                auto ret = deflate(&stream, Z_FINISH);
                if (ret != Z_OK && ret != Z_STREAM_END) {
                    throw std::runtime_error{"error deflating data with zlib"};
                }
                file.write({&outBuffer[0], outBuffer.size() - stream.avail_out});
            } while (stream.avail_in);
            file.close();
            stream.next_in = nullptr;
        }
        deflateEnd(&stream);
    }
};

using zlib_file_writer   = zlib_writer_impl<file_writer>;
using zlib_stream_writer = zlib_writer_impl<stream_writer>;

static_assert(writer_c<zlib_file_writer>);
static_assert(writer_c<zlib_stream_writer>);


}
