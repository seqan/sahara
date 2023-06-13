// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "concepts.h"

#include <functional>
#include <memory>
#include <vector>

namespace ivio {

template <typename Reader, size_t minV = (1<<12)>
class buffered_reader {
    Reader reader;
    std::vector<char> buf = []() { auto vec = std::vector<char>{}; vec.reserve(minV); return vec; }();
    int inPos{};

public:
    buffered_reader(Reader reader)
        : reader{std::move(reader)}
    {}
    buffered_reader(buffered_reader const&) = delete;
    buffered_reader(buffered_reader&& _other) = default;

private:
    auto readMore() -> bool {
        size_t lastSize = buf.size();
        if (buf.capacity() - buf.size() >= minV) {
            buf.resize(buf.capacity());
        } else {
            buf.resize(buf.capacity()*2);
        }

        auto bytes_read = reader.read(std::span{buf.data() + lastSize, buf.size() - lastSize});
        buf.resize(lastSize + bytes_read);

        return bytes_read != 0;
    }

public:
    size_t readUntil(char c, size_t lastUsed) {
        while (true) {
            auto pos = std::string_view{buf.data(), buf.size()}.find(c, lastUsed + inPos);
            if (pos != std::string_view::npos) {
                return pos - inPos;
            }

            if (!readMore()) {
                return buf.size() - inPos;
            }
        }
    }

    auto read(size_t ct) -> std::tuple<char const*, size_t> {
        while (buf.size()-inPos < ct) {
            if (!readMore()) break;
        }
        return {buf.data()+inPos, buf.size()-inPos};
    }

    void dropUntil(size_t i) {
        i = i + inPos;
        if (i < minV) {
            inPos = i;
            return;
        }
        std::copy(begin(buf)+i, end(buf), begin(buf));
        buf.resize(buf.size()-i);
        inPos = 0;
    }

    bool eof(size_t i) const {
        return i+inPos == buf.size();
    }

    auto string_view(size_t start, size_t end) -> std::string_view {
        return std::string_view{buf.data()+start+inPos, buf.data()+end+inPos};
    }

    auto tell() const -> size_t requires Seekable<Reader> {
        return reader.tell();
    }
    void seek(size_t offset) requires Seekable<Reader> {
        return reader.seek(offset);
    }

};
template <size_t minV, typename Reader>
auto make_buffered_reader(Reader&& reader) {
    return buffered_reader<Reader, minV>{std::forward<Reader>(reader)};
}

struct VarBufferedReader {
    VarBufferedReader() = default;
    VarBufferedReader(VarBufferedReader const&) = delete;
    VarBufferedReader(VarBufferedReader&&) = default;

    template <BufferedReadable T>
    VarBufferedReader(T&& t) {
        auto sptr = std::make_shared<T>(std::forward<T>(t));
        readUntil = [sptr] (char c, size_t s) {
            return sptr->readUntil(c, s);
        };
        dropUntil = [sptr] (size_t s) {
            return sptr->dropUntil(s);
        };
        read = [sptr] (size_t s) {
            return sptr->read(s);
        };
        eof = [sptr] (size_t s) {
            return sptr->eof(s);
        };
        string_view = [sptr] (size_t s, size_t e) {
            return sptr->string_view(s, e);
        };
        tell = [sptr] () -> size_t {
            if constexpr (Seekable<T>) {
                return sptr->tell();
            }
            throw std::runtime_error("this file format does not support tell/seek(1)");
        };
        seek = [sptr](size_t offset) {
            if constexpr (Seekable<T>) {
                sptr->seek(offset);
                return;
            }
            throw std::runtime_error("this file format does not support tell/seek(2)");
        };
    }
    template <Readable T>
    VarBufferedReader(T&& t)
        : VarBufferedReader{buffered_reader{std::forward<T>(t)}}
    {}

    auto operator=(VarBufferedReader const&) -> VarBufferedReader& = delete;
    auto operator=(VarBufferedReader&&) -> VarBufferedReader& = default;

    std::function<size_t(char, size_t)>                    readUntil;
    std::function<void(size_t)>                            dropUntil;
    std::function<std::tuple<char const*, size_t>(size_t)> read;
    std::function<bool(size_t)>                            eof;
    std::function<std::string_view(size_t, size_t)>        string_view;
    std::function<size_t()>                                tell;
    std::function<void(size_t)>                            seek;
};
}
