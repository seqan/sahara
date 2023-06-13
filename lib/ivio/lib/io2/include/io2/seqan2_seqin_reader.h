#pragma once

#include "Input.h"
#include "alphabet_seqan223.h"
#include "common.h"
#include "iterator.h"
#include "typed_range.h"

#include <filesystem>
#include <seqan/seq_io.h>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/utility/range/to.hpp>
#include <string_view>


namespace io2::seqan2_seqin_io {

template <typename AlphabetS3, typename QualitiesS3>
struct record;

/**\brief A view onto a single record
 *
 * This record represents a single entry in the file.
 */
template <typename AlphabetS3, typename QualitiesS3>
struct record_view {
    using sequence_view  = typed_range<AlphabetS3>;
    using qualities_view = typed_range<QualitiesS3>;
    using record         = seqan2_seqin_io::record<AlphabetS3, QualitiesS3>;

    std::string_view id{};
    sequence_view    seq{};
    qualities_view   qual{};
};

/**\brief A copy of a fastq_io record
 */
template <typename AlphabetS3, typename QualitiesS3>
struct record {
    using sequence_t  = std::vector<AlphabetS3>;
    using qualities_t = std::vector<QualitiesS3>;
    using record_view = seqan2_seqin_io::record_view<AlphabetS3, QualitiesS3>;

    std::string id{};
    sequence_t  seq{};
    qualities_t qual{};

    record(record_view v)
        : id{v.id}
        , seq{v.seq | seqan3::ranges::to<std::vector>()}
        , qual{v.qual | seqan3::ranges::to<std::vector>()}
    {}
    record() = default;
    record(record const&) = default;
    record(record&&) = default;
    record& operator=(record const&) = default;
    record& operator=(record&&) = default;
};

/** A reader to read sequence files like fasta, fastq, genbank, embl
 *
 * Usage:
 *    auto reader = io2::seq_io::reader {
 *       .input    = _file,                   // accepts string and streams
 *       .alphabet = io2::type<seqan3::dna5>, // default dna5
 *   };
 */
template <typename AlphabetS3,
          typename QualitiesS3,
          typename ExtensionAndFormat>
struct reader {
    using record_view = seqan2_seqin_io::record_view<AlphabetS3, QualitiesS3>;
    using record      = seqan2_seqin_io::record<AlphabetS3, QualitiesS3>;

    // configurable from the outside
    io2::Input<seqan::SeqFileIn, ExtensionAndFormat::format> input;
    AlphabetS3                   alphabet_type{};
    QualitiesS3                  qualities_type{};

    static auto extensions() -> std::vector<std::string> {
        return ExtensionAndFormat::extensions();
    }

    static bool validExt(std::filesystem::path const& p) {
        return io2::validExtension(p, extensions());
    }


    // internal variables
    // storage for one record
    class {
        friend reader;
        seqan::CharString id;
        seqan::String<detail::AlphabetAdaptor<AlphabetS3>>  seq;
        seqan::String<detail::AlphabetAdaptor<QualitiesS3>> qual;

        record_view return_record;
    } storage;

    auto next() -> record_view const* {
        if (input.atEnd()) return nullptr;
        input.readRecord(storage.id, storage.seq, storage.qual);

        storage.return_record = record_view {
            .id   = detail::convert_to_view(storage.id),
            .seq  = detail::convert_to_seqan3_view(storage.seq),
            .qual = detail::convert_to_seqan3_view(storage.qual),
        };
        return &storage.return_record;
    }

    using iterator = detail::iterator<reader, record_view, record>;
    auto end() const {
        return iterator{.reader = nullptr};
    }

    friend auto begin(reader& _reader) {
        return iterator{.reader = &_reader};
    }
    friend auto end(reader const& _reader) {
        return _reader.end();
    }
};

}
