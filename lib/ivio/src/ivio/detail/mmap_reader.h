// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "file_reader.h"

#include <array>
#include <cassert>
#include <filesystem>

namespace ivio {

class mmap_reader {
protected:
    file_reader reader;
    size_t filesize_; // size of the file from inPos to the end
    char const* buffer;
    size_t inPos{};

public:
    mmap_reader(std::filesystem::path path)
        : reader{path}
        , filesize_{reader.filesize()}
        , buffer{[&]() {
            auto ptr = (char const*)mmap(nullptr, filesize_, PROT_READ, MAP_PRIVATE, reader.getFileHandler(), 0);
            return ptr;
        }()}
    {
        assert(buffer);
    }

    mmap_reader() = delete;
    mmap_reader(mmap_reader const&) = delete;
    mmap_reader(mmap_reader&& _other) noexcept
        : reader{std::move(_other.reader)}
        , filesize_{_other.filesize_}
        , buffer{_other.buffer}
        , inPos{_other.inPos}
    {
        assert(buffer);
        _other.buffer = nullptr;
    }
    auto operator=(mmap_reader const&) -> mmap_reader& = delete;
    auto operator=(mmap_reader&&) -> mmap_reader& = delete;

    ~mmap_reader() {
        if (buffer == nullptr) return;
        munmap((void*)buffer, filesize_);
    }

    size_t readUntil(char c, size_t lastUsed) {
        lastUsed += inPos;
        assert(lastUsed <= filesize_);
        auto pos = std::string_view{buffer, filesize_}.find(c, lastUsed);
        if (pos != std::string_view::npos) {
            return pos - inPos;
        }
        return filesize_ - inPos;
    }

    auto read(size_t) -> std::tuple<char const*, size_t> {
        return {buffer+inPos, filesize_-inPos};
    }

    void dropUntil(size_t i) {
        i += inPos;
        assert(i <= filesize_);
        if (i < 1'024ul * 1'024ul) {
            inPos = i;
            return;
        }

        auto mask = std::numeric_limits<size_t>::max() - 4095;
        auto diff = (i & mask);
        assert(diff <= filesize_);
        munmap((void*)buffer, diff);
        buffer = buffer + diff;
        filesize_ -= diff;
        inPos = i - diff;
    }

    bool eof(size_t i) const {
        i += inPos;
        assert(i <= filesize_);
        return filesize_ == i;
    }

    auto string_view(size_t start, size_t end) -> std::string_view {
        start += inPos;
        end   += inPos;
        assert(start <= filesize_);
        assert(end <= filesize_);
        return std::string_view{buffer+start, buffer+end};
    }

    auto tell() const -> size_t {
        return inPos + reader.filesize() - filesize_;
    }

    void seek(size_t offset) {
        if (offset >= tell()) { // Seeking forwards
            inPos += offset - tell();
            return;
        }
        // Seeking backwards requires remapping the file
        munmap((void*)buffer, filesize_);
        filesize_ = reader.filesize();
        buffer = (char const*)mmap(nullptr, filesize_, PROT_READ, MAP_PRIVATE, reader.getFileHandler(), 0);
        inPos = offset;
    }
};

static_assert(BufferedReadable<mmap_reader>);
static_assert(Seekable<mmap_reader>);

}
