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

namespace ivio {

template <>
struct writer_base<fasta::writer>::pimpl {
    using Writers = std::variant<file_writer,
                                 buffered_writer<zlib_file_writer>,
                                 stream_writer,
                                 buffered_writer<zlib_stream_writer>
                                 >;

    size_t contig_length;
    Writers writer;
    std::string buffer;
    pimpl(std::filesystem::path output, size_t contig_length, bool)
        : contig_length{contig_length}
        , writer {[&]() -> Writers {
            if (output.extension() == ".gz") {
                return buffered_writer{zlib_file_writer{file_writer{output}}};
            }
            return file_writer{output};
        }()}
    {}

    pimpl(std::ostream& output, size_t contig_length, bool compressed)
        : contig_length{contig_length}
        , writer {[&]() -> Writers {
            if (compressed) {
                return buffered_writer{zlib_stream_writer{stream_writer{output}}};
            }
            return stream_writer{output};
        }()}
    {}
};

}

namespace ivio::fasta {

writer::writer(config config_)
    : writer_base{std::visit([&](auto& p) {
        return std::make_unique<pimpl>(p, config_.length, config_.compressed);
    }, config_.output)}
{
    assert(config_.length > 0);
}

writer::~writer() = default;

void writer::write(record_view record) {
    assert(pimpl_);
    std::visit([&](auto& writer) {
        auto const& contig_length = pimpl_->contig_length;
        auto& buffer = pimpl_->buffer;
        buffer.clear();
        buffer = '>';
        buffer += record.id;
        buffer += '\n';
        auto seq = record.seq;
        while (seq.size() > contig_length) {
            buffer += seq.substr(0, contig_length);
            buffer += '\n';
            seq = seq.substr(contig_length);
        }
        if (seq.size() > 0) {
            buffer += seq;
            buffer += '\n';
        }
        writer.write(buffer);
    }, pimpl_->writer);
}

void writer::close() {
    pimpl_.reset();
}

static_assert(record_writer_c<writer>);

}
