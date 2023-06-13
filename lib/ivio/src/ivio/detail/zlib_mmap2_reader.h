// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "mmap_reader.h"

#include <zlib.h>
#include <ranges>

namespace ivio {

struct mmap_queue {
    // this mmap_queue only sufficient for 4 Terabytes of data
    size_t filesize = 1'024ul*1'024ul*1'024ul*1'024ul*4ul;
    char* buffer {[&]() {
        auto ptr = (char*)mmap(nullptr, filesize, PROT_READ | PROT_WRITE, MAP_ANONYMOUS | MAP_PRIVATE | MAP_NORESERVE, -1, 0);
        return ptr;
    }()};
    size_t inPos{};

    mmap_queue() = default;
    mmap_queue(mmap_queue&& _other) noexcept
        : filesize{_other.filesize}
        , buffer  {_other.buffer}
        , inPos   {_other.inPos}
    {
        _other.buffer = nullptr;
    }

    ~mmap_queue() {
        if (buffer == nullptr) return;
        munmap((void*)buffer, filesize);
    }


    void dropUntil(size_t i) {
        i = i + inPos;
        if (i < 1'024ul * 1'024ul) {
            inPos = i;
            return;
        }
        auto mask = std::numeric_limits<size_t>::max() - 4095;
        auto diff = i & mask;
        munmap((void*)buffer, diff);
        buffer = buffer + diff;
        filesize -= diff;
        inPos = i - diff;
    }
};


struct zlib_mmap2_reader : protected mmap_reader {
    z_stream stream {};

    mmap_queue queue;

    zlib_mmap2_reader(char const* name)
        : mmap_reader(name)
    {
        if (inflateInit2(&stream, 16 + MAX_WBITS) != Z_OK) {
            throw "error";
        }
        stream.avail_in = filesize_;
        stream.next_in  = (unsigned char*)buffer;
        stream.next_out = (unsigned char*)queue.buffer;
    }
    zlib_mmap2_reader(zlib_mmap2_reader&& _other)
        : mmap_reader{std::move(_other)}
    {
        if (inflateInit2(&stream, 16 + MAX_WBITS) != Z_OK) {
            throw "error";
        }
        stream.avail_in = filesize_;
        stream.next_in  = (unsigned char*)buffer;
        stream.next_out = (unsigned char*)queue.buffer;


    }
    zlib_mmap2_reader(zlib_mmap2_reader const&) = delete;
    ~zlib_mmap2_reader() {
        inflateEnd(&stream);
    }

    auto operator=(zlib_mmap2_reader const&) -> zlib_mmap2_reader& = delete;
    auto operator=(zlib_mmap2_reader&&) -> zlib_mmap2_reader& = delete;

    size_t filesize{};
    size_t largestSize{1'024ul*1'024ul*4};

    auto readMore() -> bool {
        while(true) {
            if (filesize < largestSize) {
                stream.avail_out = largestSize-filesize;
            } else {
                stream.avail_out += largestSize*0.25;
            }
            auto available = stream.avail_out;
            auto ret = inflate(&stream, Z_NO_FLUSH);
            mmap_reader::dropUntil((char*)stream.next_in - buffer);

            auto producedBytes = available - stream.avail_out;
            filesize += producedBytes;
            largestSize = std::max(largestSize, filesize);
            if (producedBytes > 0) {
                return true;
            }
            if (ret == Z_STREAM_END) {
                return false;
            }
        }
    }

    size_t readUntil(char c, size_t lastUsed) {
        while (true) {
            auto pos = std::string_view{queue.buffer, filesize}.find(c, lastUsed);
            if (pos != std::string_view::npos) {
                return pos;
            }
            if (!readMore()) {
                return filesize;
            }
        }
    }

    auto read(size_t ct) -> std::tuple<char const*, size_t> {
        return {nullptr, 0}; //!TODO this is broken
    }


    void dropUntil(size_t i) {
        queue.dropUntil(i);
        filesize = filesize - i;
    }

    bool eof(size_t i) const {
        return filesize == i;
    }

    auto string_view(size_t start, size_t end) -> std::string_view {
        return std::string_view{queue.buffer+start, queue.buffer+end};
    }
};

static_assert(BufferedReadable<zlib_mmap2_reader>);
}
