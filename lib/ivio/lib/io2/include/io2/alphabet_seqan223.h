#pragma once

#include <filesystem>
#include <functional>
#include <optional>
#include <ranges>
#include <seqan/seq_io.h>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

namespace io2::detail {

/**\brief Count the number of bits to represents the given number of states.
 *
 * \param nrStates
 * \return the number of bits required to represents at least nrStates different states
 *
 * \noapi
 */
constexpr unsigned numberOfBits(unsigned nrStates) {
    return nrStates < 2 ? nrStates : 1+numberOfBits(nrStates >> 1);
}


/**\brief Tag for seqan3 alphabets used in seqan2.
 *
 * \tparam AlphabetS3 specifies the seqan3 alphabet. e.g.: seqan3::dna4
 *
 * \noapi
 */
template <typename AlphabetS3>
struct AlphabetAdaptor_ {};

/**\brief The alphabet adaptor, making seqan3 alphabets available in seqan2.
 *
 * \tparam AlphabetS3 specifies the seqan3 alphabet. e.g.: seqan3::dna4
 *
 * \noapi
 */
template <typename AlphabetS3>
using AlphabetAdaptor = seqan::SimpleType<uint8_t, AlphabetAdaptor_<AlphabetS3>>;


/**\brief Provides look up tables for seqan3 alphabets
 *
 * \tparam AlphabetS3 specifies the seqan3 alphabet. e.g.: seqan3::dna4
 *
 * \noapi
 */
template <typename AlphabetS3>
struct Tables {
    /**\brief Array that maps from a char to a rank.
     * For example in seqan::dna4 'A' -> 0
     */
    static constexpr std::array<uint8_t, 256> toRank = []() {
        auto map = std::array<uint8_t, 256>{};
        for (size_t i{0}; i<map.size(); ++i) {
            map[i] = seqan3::assign_char_to((char)i, AlphabetS3{}).to_rank();
        }
        return map;
    }();

    /**\brief Array that maps from rank to char.
     * For example in seqan::dna4 0 -> 'A'
     */
    static constexpr std::array<uint8_t, AlphabetS3::alphabet_size> toChar = []() {
        auto map = std::array<uint8_t, AlphabetS3::alphabet_size>{};
        for (size_t i{0}; i < map.size(); ++i) {
            map[i] = seqan3::assign_rank_to(i, AlphabetS3{}).to_char();
        }
        return map;
    }();
};

}

namespace seqan {

/**\brief Converts a char to a rank according to AlphabetAdaptor/seqan3.
 *
 * Required by seqan2.
 *
 * \tparam AlphabetS3 specifies the seqan3 alphabet. e.g.: seqan3::dna4
 *
 * \noapi
 */
template <typename AlphabetS3>
constexpr inline void assign(io2::detail::AlphabetAdaptor<AlphabetS3> & target, char c_source) {
    target.value = io2::detail::Tables<AlphabetS3>::toRank[(size_t) c_source];
}

/**\brief Converts a rank to a char according to AlphabetAdaptor/seqan3.
 *
 * Required by seqan2.
 *
 * \tparam AlphabetS3 specifies the seqan3 alphabet. e.g.: seqan3::dna4
 *
 * \noapi
 */
template <typename AlphabetS3>
constexpr inline void assign(char & c_target, io2::detail::AlphabetAdaptor<AlphabetS3> const & source) {
    c_target = io2::detail::Tables<AlphabetS3>::toChar[(size_t)source.value];
}

}
