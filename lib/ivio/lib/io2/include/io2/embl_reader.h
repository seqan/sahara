#pragma once

#include "seqan2_seqin_reader.h"

namespace io2::embl_io {

template <typename AlphabetS3, typename QualitiesS3>
using record_view = seqan2_seqin_io::record_view<AlphabetS3, QualitiesS3>;

template <typename AlphabetS3, typename QualitiesS3>
using record      = seqan2_seqin_io::record<AlphabetS3, QualitiesS3>;

struct ExtensionAndFormat {
    constexpr static auto format = seqan::Embl{};

    static auto extensions() -> std::vector<std::string> {
        static auto list = std::vector<std::string> {
            ".embl",
        };
        return list;
    }
};

template <typename AlphabetS3 = seqan3::dna5,
          typename QualitiesS3 = seqan3::phred42>
using reader = seqan2_seqin_io::reader<AlphabetS3, QualitiesS3, ExtensionAndFormat>;

}
