#pragma once

#include "common.h"
#include "alphabet_seqan223.h"
#include "Output.h"
#include "seq_io_reader.h"
#include "typed_range.h"

#include <seqan/seq_io.h>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>

namespace io2::seqan2_seqout_io {

template <typename AlphabetS3,
          typename QualitiesS3,
          typename ExtensionAndFormat>
struct writer {
    // configurable from the outside
    Output<seqan::SeqFileOut, ExtensionAndFormat::format> output{};
    AlphabetS3  alphabet_type{};
    QualitiesS3 qualities_type{};

    struct record {
        typed_range<char>        id{};
        typed_range<AlphabetS3>  seq{};
        typed_range<QualitiesS3> qual{};
    };

    static auto extensions() -> std::vector<std::string> {
        return ExtensionAndFormat::extensions();
    }

    static bool validExt(std::filesystem::path const& p) {
        return io2::validExtension(p, extensions());
    }


    void write(record _record) {
        auto id   = detail::convert_to_seqan2_string(_record.id);
        auto seq  = detail::convert_to_seqan2_alphabet(_record.seq);
        auto qual = detail::convert_to_seqan2_qualities(_record.qual);

        writeRecord(*output.fileOut, id, seq, qual);
    }

    template <typename record_like>
    void write(record_like const& _record) {
        write(record{
            .id   = _record.id,
            .seq  = _record.seq,
            .qual = _record.qual,
        });
    }
};

}
