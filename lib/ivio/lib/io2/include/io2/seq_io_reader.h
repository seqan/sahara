#pragma once

#include "fasta_reader.h"
#include "fastq_reader.h"
#include "embl_reader.h"
#include "genbank_reader.h"


namespace io2::seq_io {

template <typename AlphabetS3, typename QualitiesS3>
using record_view = seqan2_seqin_io::record_view<AlphabetS3, QualitiesS3>;

template <typename AlphabetS3, typename QualitiesS3>
using record      = seqan2_seqin_io::record<AlphabetS3, QualitiesS3>;


/** A reader to read sequence files like fasta, fastq, genbank, embl
 *
 * Usage:
 *    auto reader = io2::seq_io::reader {
 *       .input     = _file,                     // accepts string and streams
 *       .alphabet  = io2::type<seqan3::dna5>,   // default dna5
 *       .qualities = io2::type<seqan::phred42>, // default phred42
 *   };
 */
template <typename AlphabetS3 = seqan3::dna5,
          typename QualitiesS3 = seqan3::phred42>
struct reader {
    using record_view = seq_io::record_view<AlphabetS3, QualitiesS3>;
    using record      = seq_io::record<AlphabetS3, QualitiesS3>;

    // configurable from the outside
//    io2::Input<seqan::SeqFileIn> input;
    std::filesystem::path input{};
    AlphabetS3            alphabet_type{};
    QualitiesS3           qualities_type{};


    fasta_io::reader<AlphabetS3>              fasta{};
    fastq_io::reader<AlphabetS3, QualitiesS3> fastq{};
    fastq_io::reader<AlphabetS3, QualitiesS3> embl{};
    fastq_io::reader<AlphabetS3, QualitiesS3> genbank{};

    static auto extensions() -> std::vector<std::string> {
        static auto list = [&]() {
            auto l = std::vector<std::string>{};
            for (auto e : decltype(fasta)::extensions()) l.push_back(e);
            for (auto e : decltype(fastq)::extensions()) l.push_back(e);
            for (auto e : decltype(embl)::extensions()) l.push_back(e);
            for (auto e : decltype(genbank)::extensions()) l.push_back(e);
            return l;
        }();
        return list;
    }

    static bool validExt(std::filesystem::path const& p) {
        return io2::validExtension(p, extensions());
    }


    // internal variables
    // storage for one record
    class {
        friend reader;
        record_view return_record;
    } storage{};


    std::function<record_view const*()> next{};
    void* ctor = [this]() {
        if (validExtension(input, decltype(fasta)::extensions())) {
            fasta.input = input;
            next = [this]() -> record_view const* {
                auto r = fasta.next();
                if (r == nullptr) return nullptr;
                storage.return_record = record_view {
                    .id  = r->id,
                    .seq = r->seq,
                };
                return &storage.return_record;
            };
        } else if (validExtension(input, decltype(fastq)::extensions())) {
            fastq.input = input;
            next = [this]() -> record_view const* {
                auto r = fastq.next();
                if (r == nullptr) return nullptr;
                storage.return_record = record_view {
                    .id   = r->id,
                    .seq  = r->seq,
                    .qual = r->qual,
                };
                return &storage.return_record;
            };
        } else if (validExtension(input, decltype(embl)::extensions())) {
            embl.input = input;
            next = [this]() -> record_view const* {
                auto r = fastq.next();
                if (r == nullptr) return nullptr;
                storage.return_record = record_view {
                    .id   = r->id,
                    .seq  = r->seq,
                    .qual = r->qual,
                };
                return &storage.return_record;
            };
        } else if (validExtension(input, decltype(genbank)::extensions())) {
            genbank.input = input;
            next = [this]() -> record_view const* {
                auto r = fastq.next();
                if (r == nullptr) return nullptr;
                storage.return_record = record_view {
                    .id   = r->id,
                    .seq  = r->seq,
                    .qual = r->qual,
                };
                return &storage.return_record;
            };
        } else if (!input.empty()) {
            throw std::runtime_error("unknown file format");
        }
        return nullptr;
    }();

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
