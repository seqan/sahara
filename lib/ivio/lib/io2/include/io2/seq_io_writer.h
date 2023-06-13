#pragma once

#include "fastq_writer.h"
#include "fasta_writer.h"

#include "common.h"
#include "alphabet_seqan223.h"
#include "Output.h"
#include "seq_io_reader.h"
#include "typed_range.h"

#include <seqan/seq_io.h>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>

namespace io2::seq_io {

template <typename AlphabetS3 = seqan3::dna5,
          typename QualitiesS3 = seqan3::phred42>
struct writer {
    // configurable from the outside
    std::filesystem::path output{};
    AlphabetS3 alphabet_type{};

    fasta_io::writer<AlphabetS3, QualitiesS3> fasta{};
    fastq_io::writer<AlphabetS3, QualitiesS3> fastq{};


    struct record {
        typed_range<char>        id{};
        typed_range<AlphabetS3>  seq{};
        typed_range<QualitiesS3> qual{};
    };

    void write(record _record) {
        writeImpl(_record);
    }

    template <typename record_like>
    void write(record_like const& _record) {
        write(record{
            .id   = _record.id,
            .seq  = _record.seq,
            .qual = _record.qual,
        });
    }

    std::function<void(record _record)> writeImpl{};
    void* ctor = [this]() {
        if (validExtension(output, decltype(fasta)::extensions())) {
            fasta.output = output;
            writeImpl = [this](record _record) {
                fasta.write(_record);
            };
        } else if (validExtension(output, decltype(fastq)::extensions())) {
            fastq.output = output;
            writeImpl = [this](record _record) {
                fasta.write(_record);
            };
        } else if (!output.empty()) {
            throw std::runtime_error("unknown file format");
        }
        return nullptr;
    }();

};

}
