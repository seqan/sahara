// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include "../detail/buffered_writer.h"
#include "../detail/file_writer.h"
#include "../detail/stream_writer.h"
#include "../detail/zlib_file_writer.h"
#include "writer.h"

#include <cassert>
#include <charconv>

template <>
struct ivio::writer_base<ivio::sam::writer>::pimpl {
    using Writers = std::variant<file_writer,
                                 buffered_writer<zlib_file_writer>,
                                 stream_writer,
                                 buffered_writer<zlib_stream_writer>
                                 >;

    ivio::sam::writer::config config;
    Writers writer;

    pimpl(std::filesystem::path output, bool)
        : writer {[&]() -> Writers {
            if (output.extension() == ".gz") {
                return buffered_writer{zlib_file_writer{file_writer{output}}};
            }
            return file_writer{output};
        }()}
    {}

    pimpl(std::ostream& output, bool compressed)
        : writer {[&]() -> Writers {
            if (compressed) {
                return buffered_writer{zlib_stream_writer{stream_writer{output}}};
            }
            return stream_writer{output};
        }()}
    {}
};

namespace ivio::sam {

writer::writer(config config_)
    : writer_base{std::visit([&](auto& p) {
        return std::make_unique<pimpl>(p, config_.compressed);
    }, config_.output)}
{
    // write header
    std::visit([&](auto& writer) {
        for (auto const& value : config_.header) {
            writer.write(value);
            writer.write(std::string_view{"\n"});
        }
    }, pimpl_->writer);
}

writer::~writer() = default;

void writer::write(record_view record) {
    assert(pimpl_);

    static thread_local auto ss = std::string{};
    ss = record.qname; ss += '\t';
    ss += std::to_string(record.flag); ss += '\t';
    ss += record.rname; ss += '\t';
    ss += std::to_string(record.pos); ss += '\t';
    ss += std::to_string(record.mapq); ss += '\t';
    ss += record.cigar; ss += '\t';
    ss += record.rnext; ss += '\t';
    ss += std::to_string(record.pnext); ss += '\t';
    ss += std::to_string(record.tlen); ss += '\t';
    ss += record.seq; ss += '\t';
    ss += record.qual; ss += '\t';
    ss += record.tags; ss += '\n';

    std::visit([&](auto& writer) {
       writer.write(ss);
    }, pimpl_->writer);
}

void writer::close() {
    pimpl_.reset();
}

static_assert(record_writer_c<writer>);


}
