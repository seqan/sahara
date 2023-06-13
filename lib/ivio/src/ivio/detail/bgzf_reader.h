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

#include <algorithm>
#include <ranges>
#include <zlib.h>

namespace ivio {

/*!\brief Convert a little_endian_to_host byte-order
 * \param  in   The input value to convert.
 * \returns the converted value in host byte-order.
 */
auto little_endian_to_host(std::integral auto const in) noexcept {
    if constexpr (sizeof(in) == 1) {
        return in;
    } else if constexpr (sizeof(in) == 2) {
        return le16toh(in);
    } else if constexpr (sizeof(in) == 4) {
        return le32toh(in);
    } else if constexpr (sizeof(in) == 8) {
        return le64toh(in);
    } else {
        []<bool flag=false>() {
            static_assert(flag, "Conversion only for 1, 2, 4 and 8 byte types possible");
        }();
    }
}


inline constexpr char _magic_bgzf_header[] =
    // GZip header
    "\x1f\x8b\x08"
    // FLG[MTIME         ] XFL OS [XLEN  ]
    "\x04\x00\x00\x00\x00\x00\xff\x06\x00"
    // B   C [SLEN  ][BSIZE ]
    "\x42\x43\x02\x00\x00\x00";

inline constexpr auto magic_bgzf_header = std::string_view{_magic_bgzf_header, sizeof(_magic_bgzf_header)};

template <typename T>
inline auto bgzfUnpack(char const* buffer) -> T {
    T v;
    std::uninitialized_copy(buffer, buffer + sizeof(v), reinterpret_cast<char*>(&v));
    if constexpr (std::integral<T>) {
        return little_endian_to_host(v);
    } else {
        return v;
    }
}

struct ZlibContext {
    inline static constexpr size_t BlockHeaderLength = magic_bgzf_header.size();
    inline static constexpr size_t BlockFooterLength = 8;

    z_stream strm{};

    ZlibContext() {
        constexpr auto GzipWindowBits = -15; // no zlib header
        auto status = inflateInit2(&strm, GzipWindowBits);
        if (status != Z_OK) {
            throw std::runtime_error{"BGZF inflateInit2() failed"};
        }
    }

    ZlibContext(ZlibContext const&) = delete;
    ZlibContext(ZlibContext&&) = delete;

    ~ZlibContext() noexcept(false) {
        auto status = inflateEnd(&strm);
        if (status != Z_OK) {
            throw std::runtime_error{"BGZF inflateEnd() failed"};
        }
    }

    auto operator=(ZlibContext const&) -> ZlibContext = delete;
    auto operator=(ZlibContext&&) -> ZlibContext = delete;


    void reset() {
        auto status = inflateReset(&strm);
        if (status != Z_OK) {
            throw std::runtime_error{"BGZF inflateReset() failed"};
        }
    }

    size_t decompressedSize(std::span<char const> in) const {
        return bgzfUnpack<uint32_t>(in.data() + in.size() - 4);
    }

    size_t decompressBlock(std::span<char const> in, std::span<char> out) {
        reset();

        if (in.size() < BlockFooterLength) {
            throw std::runtime_error{"BGZF block too short. " + std::to_string(in.size())};
        }

//        if (!detail::bgzf_compression::validate_header(std::span{srcBegin, srcLength})) {
//            throw io_error("Invalid BGZF block header.");
//        }
        strm.next_in   = (Bytef *)in.data();
        strm.next_out  = (Bytef *)out.data();
        strm.avail_in  = in.size() - BlockFooterLength;
        strm.avail_out = out.size();

        auto status = inflate(&strm, Z_FINISH);
        if (status != Z_STREAM_END) {
            throw std::runtime_error{"Inflation failed. Decompressed BGZF data is too big. " + std::to_string(strm.avail_in) + " " + std::to_string(strm.avail_out)};
        }

        // Compute and check checksum
        unsigned crc = crc32(crc32(0, nullptr, 0), (Bytef *)out.data(), out.size() - strm.avail_out);
        unsigned ecrc = bgzfUnpack<uint32_t>(in.data() + in.size() - 8);
        if (ecrc != crc)
            throw std::runtime_error{"BGZF wrong checksum." + std::to_string(ecrc) + " " + std::to_string(crc)};

        // Check uncompressed data length
        auto dlen = bgzfUnpack<uint32_t>(in.data() + in.size() - 4);
        if (dlen != out.size() - strm.avail_out)
            throw std::runtime_error{"BGZF size mismatch."};

        return out.size() - strm.avail_out;

    }
};

struct bgzf_reader {
    VarBufferedReader reader;
    ZlibContext       zlibCtx;

    bgzf_reader(VarBufferedReader reader)
        : reader{std::move(reader)}
    {}

    bgzf_reader(bgzf_reader const&) = delete;
    bgzf_reader(bgzf_reader&& _other)
        : reader{std::move(_other.reader)}
    {}

    auto operator=(bgzf_reader const&) -> bgzf_reader& = delete;
    auto operator=(bgzf_reader&&) -> bgzf_reader& = delete;


    ~bgzf_reader() = default;

    size_t read(std::ranges::sized_range auto&& range) {
        while(true) {
            auto [ptr, avail_in] = reader.read(18);
            if (avail_in == 0) return 0;
            if (avail_in < 18) throw std::runtime_error{"failed reading (1)"};

            size_t compressedLen = bgzfUnpack<uint16_t>(ptr + 16) + 1u;
            std::tie(ptr, avail_in) = reader.read(compressedLen);
            if (avail_in < compressedLen) throw std::runtime_error{"failed reading (2)"};

            assert(range.size() >= (1<<16));

            size_t size = zlibCtx.decompressBlock({ptr+18, compressedLen-18}, {range.data(), range.size()});
            reader.dropUntil(compressedLen);
            return size;
        }
    }
};

static_assert(Readable<bgzf_reader>);
}
