#pragma once

#include "seqan2_seqout_writer.h"
#include "fasta_reader.h"

namespace io2::fasta_io {

struct ExtensionAndFormat {
    constexpr static auto format = seqan::Fastq{};

    static auto extensions() -> std::vector<std::string> {
        static auto list = std::vector<std::string> {
            ".fasta",
            ".fa",
        };
        return list;
    }
};


template <typename AlphabetS3 = seqan3::dna5,
          typename QualitiesS3 = seqan3::phred42>
using writer = seqan2_seqout_io::writer<AlphabetS3, QualitiesS3, ExtensionAndFormat>;

}
