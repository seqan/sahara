// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include "../detail/buffered_reader.h"
#include "../detail/file_reader.h"
#include "../detail/mmap_reader.h"
#include "../detail/stream_reader.h"
#include "../detail/zlib_file_reader.h"
#include "../detail/zlib_mmap2_reader.h"
#include "reader.h"

static_assert(std::ranges::range<ivio::fasta::reader>, "reader must be a range (unittest)");
static_assert(ivio::record_reader_c<ivio::fasta::reader>, "must fulfill the record_reader concept (unittest)");
static_assert(ivio::Seekable<ivio::fasta::reader>, "must fulfill the Seekable concept (unittest)");

namespace ivio {

template <>
struct reader_base<fasta::reader>::pimpl {
    VarBufferedReader ureader;
    size_t lastUsed{};
    std::string s;
    faidx::record faidx_record{{}, 0, 0, 0, 0};

    pimpl(std::filesystem::path file, bool)
        : ureader {[&]() -> VarBufferedReader {
            if (file.extension() == ".gz") {
                return zlib_reader{mmap_reader{file.c_str()}};
            }
            return mmap_reader{file.c_str()};
        }()}
    {}
    pimpl(std::istream& file, bool compressed)
        : ureader {[&]() -> VarBufferedReader {
            if (!compressed) {
                return stream_reader{file};
            } else {
                return zlib_reader{stream_reader{file}};
            }
        }()}
    {}
};
}

namespace ivio::fasta {

reader::reader(config const& config_)
    : reader_base{std::visit([&](auto& p) {
        return std::make_unique<pimpl>(p, config_.compressed);
    }, config_.input)}
{}

reader::~reader() = default;

auto reader::next() -> std::optional<record_view> {
    assert(pimpl_);

    auto& ureader  = pimpl_->ureader;
    auto& lastUsed = pimpl_->lastUsed;
    auto& s        = pimpl_->s;

    auto startId = ureader.readUntil('>', lastUsed);
    if (ureader.eof(startId)) return std::nullopt;
    ureader.dropUntil(startId+1);

    auto endId = ureader.readUntil('\n', 0);
    if (ureader.eof(endId)) return std::nullopt;

    auto startSeq = endId+1;

    // convert into dense string representation
    s.clear();
    {
        auto s2 = startSeq;
        do {
            auto s1 = s2;
            s2 = ureader.readUntil('\n', s1);
            s += ureader.string_view(s1, s2);
            if (!ureader.eof(s2)) {
                s2 += 1;
            }
        } while (!ureader.eof(s2) and ureader.string_view(s2, s2+1)[0] != '>');
        lastUsed = s2;
    }


    return record_view {
        .id  = ureader.string_view(0, endId),
        .seq = s,
    };
}

void reader::close() {
    pimpl_.reset();
}

auto reader::tell() const -> size_t {
    assert(pimpl_);

    auto& ureader  = pimpl_->ureader;
    return ureader.tell();
}

void reader::seek(size_t offset) {
    assert(pimpl_);

    auto& ureader  = pimpl_->ureader;

    ureader.seek(offset);
    ureader.dropUntil(0);
    pimpl_->lastUsed = 0;
}

auto reader::tell_faidx() const -> faidx::record {
    assert (pimpl_);
    return pimpl_->faidx_record;
}

void reader::seek_faidx(faidx::record const& faidx) {
    assert (pimpl_);

    auto& ureader  = pimpl_->ureader;

    size_t numberOfLines = faidx.length / faidx.linebases;

    if (faidx.length % faidx.linebases != 0) { // last line is not full
        numberOfLines += 1;
    }

    size_t numberOfNewLineCharacters = numberOfLines * (faidx.linewidth - faidx.linebases);
    size_t offset = faidx.offset + numberOfNewLineCharacters + faidx.length;
    if (!ureader.eof(offset - ureader.tell())) {
        seek(offset);
    } else {
        seek(0);
    }
}

}
