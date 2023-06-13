#pragma once

#include "Input.h"
#include "alphabet_seqan223.h"
#include "common.h"
#include "iterator.h"
#include "typed_range.h"

#include <filesystem>
#include <optional>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/utility/range/all.hpp>
#include <string_view>


namespace io2::vcf_io {

enum class format {
    Vcf,
    Bcf,
};

/**
 * \noapi
 */
void convert_format(format _format, auto&& cb) {
    switch(_format) {
    case format::Vcf: cb(seqan::Vcf()); break;
    case format::Bcf: cb(seqan::Bcf()); break;
    }
}

/**\brief A view onto a single record
 *
 * This record represents a single entry in the file.
 */
template <typename AlphabetS3>
struct record_view {
    // views for string types
    using sequence_view  = decltype(detail::convert_string_to_seqan3_view<AlphabetS3>({}));

    int32_t                             rID{};
    int32_t                             beginPos{};
    std::string_view                    id{};
    sequence_view                       ref{};
    sequence_view                       alt{};
    float                               qual{};
    std::string_view                    filter{};
    std::string_view                    info{};
    std::string_view                    format{};
    sized_typed_range<std::string_view> genotypeInfos{};
};

/**\brief A copy of a seq_io record
 */
template <typename AlphabetS3>
struct record {
    // views for string types
    using sequence_t  = std::vector<AlphabetS3>;

    int32_t                  rID{};
    int32_t                  beginPos{};
    std::string              id{};
    sequence_t               ref{};
    sequence_t               alt{};
    float                    qual{};
    std::string              filter{};
    std::string              info{};
    std::string              format{};
    std::vector<std::string> genotypeInfos{};


    record(record_view<AlphabetS3> v)
        : rID{v.rID}
        , beginPos{v.beginPos}
        , id{v.id}
        , ref{v.ref | seqan3::ranges::to<std::vector>()}
        , alt{v.alt | seqan3::ranges::to<std::vector>()}
        , qual{v.qual}
        , filter{v.filter}
        , info{v.info}
        , format{v.format}
        , genotypeInfos{v.genotypeInfos | seqan3::ranges::to<std::vector<std::string>>()}
    {}

    record() = default;
    record(record const&) = default;
    record(record&&) = default;
    record& operator=(record const&) = default;
    record& operator=(record&&) = default;
};


/** A reader to vcf and bcf files
 *
 * Usage:
 *    auto reader = io2::vcf_io::reader {
 *       .input    = _file,                   // accepts string and streams
 *       .alphabet = io2::type<seqan3::dna5>, // default dna5
 *   };
 */
template <typename AlphabetS3 = seqan3::dna5>
struct reader {
    using record_view = vcf_io::record_view<AlphabetS3>;
    using record      = vcf_io::record<AlphabetS3>;

    // configurable from the outside
    io2::Input<seqan::VcfFileIn> input;
    [[no_unique_address]] detail::empty_class<AlphabetS3>  alphabet{};

    // internal variables
    /* storage for one record
     * \noapi
     */
    struct {
        seqan::VcfRecord seqan2_record;
        record_view      return_record;

        seqan::VcfHeader seqan2_header;
        std::decay_t<decltype(sampleNames(context(seqan::VcfFileIn{})))> seqan2_sampleNames;
    } storage{};

    /* A fake constructor
     *
     * We are using this function as a constructor, to keep the designated initializer
     * syntax intact.
     *
     * \noapi
     */
    [[no_unique_address]] detail::empty_class<nullptr_t> _fakeConstructor = [this]() {
        seqan::readHeader(storage.seqan2_header, *input.fileIn);
        storage.seqan2_sampleNames = sampleNames(context(*input.fileIn));
        return nullptr;
    }();




    /** accessing header information
     * \experimentalapi
     */
    struct {
        reader& reader_;

        auto samples() -> sized_typed_range<std::string_view> {
            auto& names = reader_.storage.seqan2_sampleNames;

            return std::ranges::subrange{begin(names), begin(names) + length(names)}
                | std::views::transform([](seqan::String<char> const& s) {
                    return detail::convert_to_view(s);
                });
        }
        auto entries() const -> sized_typed_range<std::tuple<std::string_view, std::string_view>> {
            return {std::ranges::subrange{begin(reader_.storage.seqan2_header), begin(reader_.storage.seqan2_header) + length(reader_.storage.seqan2_header)}
                | std::views::transform([](seqan::VcfHeaderRecord const& v) {
                    return std::make_tuple(detail::convert_to_view(v.key),
                                            detail::convert_to_view(v.value));
                })};
        }
        auto entries(std::string const& _filter) const -> typed_range<std::string_view> {
            return {std::ranges::subrange{begin(reader_.storage.seqan2_header), begin(reader_.storage.seqan2_header) + length(reader_.storage.seqan2_header)}
                | std::views::filter([_filter](seqan::VcfHeaderRecord const& v) {
                    return toCString(v.key) == _filter;
                }) | std::views::transform([](seqan::VcfHeaderRecord const& v) {
                    return detail::convert_to_view(v.value);
                })};
        }
    } header {*this};

    auto next() -> record_view const* {
        if (input.atEnd()) return nullptr;
        input.readRecord(storage.seqan2_record);


        auto const& r = storage.seqan2_record; // shorter name
        storage.return_record = record_view {
            .rID      = r.rID,
            .beginPos = r.beginPos,
            .id       = detail::convert_to_view(r.id),
            .ref      = detail::convert_string_to_seqan3_view<AlphabetS3>(r.ref),
            .alt      = detail::convert_string_to_seqan3_view<AlphabetS3>(r.alt),
            .qual     = r.qual,
            .filter   = detail::convert_to_view(r.filter),
            .info     = detail::convert_to_view(r.info),
            .format   = detail::convert_to_view(r.format),
            .genotypeInfos = {r.genotypeInfos | std::views::transform([](auto& s) { return detail::convert_to_view(s); })},
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
