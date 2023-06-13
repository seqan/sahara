// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "file_writer.h"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <cstring>
#include <ranges>
#include <zlib.h>

namespace ivio {

namespace bgzf_writer::detail {

/*!\brief Convert a host_to_little_endian byte-order
 * \param  in   The input value to convert.
 * \returns the converted value in host byte-order.
 */
auto host_to_little_endian(std::integral auto const in) noexcept {
    if constexpr (sizeof(in) == 1) {
        return in;
    } else if constexpr (sizeof(in) == 2) {
        return htole16(in);
    } else if constexpr (sizeof(in) == 4) {
        return htole32(in);
    } else if constexpr (sizeof(in) == 8) {
        return htole64(in);
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
inline auto bgzfPack(T v, char* buffer) -> size_t {
    if constexpr (std::integral<T>) {
        v = host_to_little_endian(v);
    }
    std::memcpy(buffer, reinterpret_cast<char*>(&v), sizeof(v));
    return sizeof(v);
}

struct ZlibContext {
    inline static constexpr size_t BlockHeaderLength = magic_bgzf_header.size();
    inline static constexpr size_t BlockFooterLength = 8;

    z_stream stream{};

    ZlibContext() {
        constexpr auto GzipWindowBits = -15; // no zlib header
        auto status = deflateInit2(&stream, 6, Z_DEFLATED, GzipWindowBits, 8, Z_DEFAULT_STRATEGY);
        if (status != Z_OK) {
            throw "BGZF deflateInit2() failed";
        }
    }

    ~ZlibContext() noexcept(false) {
        auto status = deflateEnd(&stream);
        if (status != Z_OK) {
            throw "BGZF inflateEnd() failed";
        }

    }

    void reset() {
        auto status = deflateReset(&stream);
        if (status != Z_OK) {
            throw "BGZF deflateReset() failed";
        }
    }

    size_t compressBlock(std::span<char const> in, std::span<char> out) {
        reset();

/*        if (in.size() < BlockFooterLength) {
            throw "BGZF block too short. " + std::to_string(in.size());
        }*/

//        if (!detail::bgzf_compression::validate_header(std::span{srcBegin, srcLength})) {
//            throw io_error("Invalid BGZF block header.");
//        }
        for (size_t i{0}; i < sizeof(_magic_bgzf_header); ++i) {
            out[i] = _magic_bgzf_header[i];
        }
        stream.next_in   = (Bytef *)in.data();
        stream.next_out  = (Bytef *)out.data()+18;
        stream.avail_in  = in.size();
        stream.avail_out = out.size()-18-BlockFooterLength;

        auto status = deflate(&stream, Z_FINISH);
        if (status != Z_STREAM_END) {
            throw "Deflation failed. compressed BGZF data is too big. " + std::to_string(stream.avail_in) + " " + std::to_string(stream.avail_out);
        }

        auto length = out.size() - stream.avail_out;
        bgzf_writer::detail::bgzfPack(static_cast<uint16_t>(length-1ul), &out[16]);

        uint32_t crc = crc32(crc32(0, nullptr, 0), (Bytef *)in.data(), in.size());

        bgzf_writer::detail::bgzfPack(static_cast<uint32_t>(crc), &out[length-8ul]);
        bgzf_writer::detail::bgzfPack(static_cast<uint32_t>(in.size()), &out[length-4ul]);


        // Compute and check checksum
        //unsigned crc = crc32(crc32(0, nullptr, 0), (Bytef *)out.data(), out.size() - stream.avail_out);
        //unsigned ecrc = bgzfUnpack<uint32_t>(in.data() + in.size() - 8);
        //if (ecrc != crc)
        //    throw "BGZF wrong checksum." + std::to_string(ecrc) + " " + std::to_string(crc);

        //// Check uncompressed data length
        //auto dlen = bgzfUnpack<uint32_t>(in.data() + in.size() - 4);
        //if (dlen != out.size() - stream.avail_out)
        //    throw "BGZF size mismatch.";

        return length;

    }
};

}

template <writer_c Writer>
struct bgzf_writer_impl {
    Writer      file;
    bgzf_writer::detail::ZlibContext zlibCtx;

    std::vector<char> buffer{};
    std::vector<char> outBuffer{};

    template <typename T>
    bgzf_writer_impl(T&& name)
        : file(std::forward<T>(name))
    {}

    bgzf_writer_impl(bgzf_writer_impl&& _other)
        : file{std::move(_other.file)}
    {}

    ~bgzf_writer_impl() {
        close();
    }

    static constexpr auto fullLength = 65280;
    auto write(std::span<char const> out) -> size_t {
        auto oldSize = buffer.size();
        buffer.resize(buffer.size() + out.size());
        std::ranges::copy(out, buffer.data() + oldSize);

        auto writeData = [&](std::span<char const> v) {
            outBuffer.resize(fullLength);
            auto length = zlibCtx.compressBlock(v, outBuffer);
            outBuffer.resize(length);

            // write to file
            file.write(outBuffer);
        };

        while (buffer.size() >= fullLength) {
            writeData({buffer.data(), buffer.data() + fullLength});

            // move left over data to the beginning
            std::memcpy(buffer.data(), buffer.data() + fullLength, buffer.size() - fullLength);//!TODO maybe with ranges::copy?
            buffer.resize(buffer.size() - fullLength);
        }
        return out.size();
    }

    void close() {
        assert(buffer.size() < fullLength);
        if (!buffer.empty()) {
            outBuffer.resize(fullLength);
            auto length = zlibCtx.compressBlock({buffer.data(), buffer.size()}, outBuffer);
            outBuffer.resize(length);

            // write to file
            file.write(outBuffer);

            buffer.clear();
        }
        file.close();
    }
};

using bgzf_file_writer   = bgzf_writer_impl<file_writer>;
//using bgzf_mmap_reader   = bgzf_reader_impl<mmap_reader>;
//using bgzf_stream_reader = bgzf_reader_impl<buffered_reader<stream_reader>>;
//
//static_assert(reader_c<bgzf_file_reader>);
//static_assert(reader_c<bgzf_mmap_reader>);
//static_assert(reader_c<bgzf_stream_reader>);
}
