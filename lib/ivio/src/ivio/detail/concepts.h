// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace ivio {

/* \brief describes a Reader that supports `read`
 * - read: reads data into a given buffer
 */
template <typename T>
concept Readable = requires(T t) {
    { t.read(std::span<char>{}) } -> std::same_as<size_t>;
};


/* \brief a Reader that supports the functions:
 * - readUntil: reads from a given point until a certain character is found or eof
 * - dropUntil: allows the internal buffer to drop all memory until a certain point
 * - eof: check if a value is at the end of the file
 * - string_view: returns a view to a buffer range as returned by readUntil
 */
template <typename T>
concept BufferedReadable = requires(T t) {
    { t.readUntil(char{}, size_t{}) } -> std::same_as<size_t>;
    { t.dropUntil(size_t{}) } -> std::same_as<void>;
    { t.read(size_t{}) } -> std::same_as<std::tuple<char const*, size_t>>;
    { t.eof(size_t{}) } -> std::same_as<bool>;
    { t.string_view(size_t{}, size_t{}) } -> std::same_as<std::string_view>;
};

/* \brief a Reader that is seekable
 * This means that random access is possible.
 * tell(): returns the current position
 * seek(): moves the internal position to a certain position
 */
template <typename T>
concept Seekable = requires(T t) {
    { t.tell() } -> std::same_as<size_t>;
    { t.seek(size_t{}) } -> std::same_as<void>;
};

/* \brief a Reader that reads record by record
 * - next: returns an optional, which is a record or std::nullopt if not available
 */
template <typename T>
concept record_reader_c = requires(T t) {
    { t.next() };
    { t.close() };
};


template <typename T>
concept writer_c = requires(T t) {
    { t.write(std::declval<std::span<char>>()) } -> std::same_as<size_t>;
    { t.close() };
//    { t.write(std::declval<char const*>()) };
};

template <typename T>
concept record_writer_c = requires(T t) {
    { t.write({}) };
    { t.close() };
};

}
