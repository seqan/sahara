#pragma once

#include "alphabet_seqan223.h"

#include <ranges>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <span>
#include <string_view>


namespace io2::detail {

/**\brief An empty base class for type deduction.
 *
 * It is meant to be used in combination with the attribute
 * [[no_unique_address]].
 * See io2::type for an example.
 *
 * \tparam t a captured type.
 *
 * \noapi
 */
template <typename t>
struct empty_class {
    using type = t;
    constexpr empty_class() noexcept = default;
    constexpr empty_class(t const &) noexcept {}
};


/**\brief Creates a view onto a seqan3::String<t>.
 *
 * \tparam t Type of the string elements.
 * \param v  Desired string which to create a view for.
 * \return   A view that points to v.
 *
 * \noapi
 */
template <typename t>
auto convert_to_view(seqan::String<t> const& v) {
    auto ptr = &*begin(v, seqan::Standard());
    return std::span<t const>(ptr, length(v));
}

/**\brief Creates a view onto a seqan::String<char>.
 *
 * \param v Desired string which to create a view for.
 * \return  A string_view that points to v.
 *
 * \noapi
 */
inline auto convert_to_view(seqan::String<char> const& v) {
    return std::string_view{&*begin(v, seqan::Standard()), length(v)};
}


/**\brief Creates a view onto a seqan::String<AlphabetAdaptor<·>>.
 *
 * \tparam AlphabetS3 A seqan3 alphabet, e.g.: seqan3::dna4.
 * \param v           Desired string which to create a view for.
 * \return            A view that points to v.
 *
 * \noapi
 */
template <typename AlphabetS3>
auto convert_to_seqan3_view(seqan::String<detail::AlphabetAdaptor<AlphabetS3>> const& v) {
    return convert_to_view(v) | std::views::transform([](auto const& v) {
        return seqan3::assign_rank_to(v.value, AlphabetS3{});
    });
}

/**\brief Creates a view onto a seqan::String<seqan::CigarElement>.
 *
 * \param v Desired string which to create a view for.
 * \return  A view over seqan3::cigar elements that points to v.
 *
 * \noapi
 */
inline auto convert_to_seqan3_view(seqan::String<seqan::CigarElement<>> const& v) {
    return convert_to_view(v) | std::views::transform([](auto const& v) {
        using namespace seqan3::literals;
        seqan3::cigar::operation letter{};
        letter.assign_char(v.operation);
        return seqan3::cigar{v.count, letter};
    });
}

/**\brief Creates a view onto a seqan::IupacString
 *
 * \param v Desired string which to create a view for.
 * \return  A view over AlpahbetS3 elements that points to v.
 *
 * \noapi
 */
template <typename AlphabetS3>
inline auto convert_iupac_to_seqan3_view(seqan::IupacString const& v) {
    return convert_to_view(v) | std::views::transform([](auto const& v) {
        return seqan3::assign_char_to(static_cast<char>(v), AlphabetS3{});
    });
}


/**\brief Creates a view onto a seqan::CharString
 *
 * \tparam AlphabetS3 A seqan3 alphabet, e.g.: seqan3::dna4.
 * \param v           Desired string which to create a view for.
 * \return            A view that points to v.
 *
 * \noapi
 */
template <typename AlphabetS3>
auto convert_string_to_seqan3_view(seqan::CharString const& v) {
    return convert_to_view(v) | std::views::transform([](auto const& v) {
        return seqan3::assign_char_to(v, AlphabetS3{});
    });
}

/**\brief Creates a seqan2 alphabet string over a seqan3 alphabet range
 *
 * \param rng A seqan3 alphabet range.
 * \return    A seqan::String over AlphabetAdaptor<·> elements.
 *
 * \noapi
 */
template <std::ranges::range rng_t>
auto convert_to_seqan2_alphabet(rng_t&& rng) {
    using AlphabetS3 = std::decay_t<decltype(*rng.begin())>;
    seqan::String<io2::detail::AlphabetAdaptor<AlphabetS3>> v;

    auto view = rng | std::views::transform([](auto c) {
            auto t = io2::detail::AlphabetAdaptor<AlphabetS3>{};
            t.value = c.to_rank();
            return t;
        });

    if constexpr (requires() {
        { std::ranges::size(rng) } -> std::same_as<size_t>;
    }) {
        resize(v, std::ranges::size(rng), seqan::Exact());
        std::ranges::copy(view, begin(v));
    } else {
        for (auto c : view) {
            appendValue(v, c);
        }
    }
    return v;
}

/**\brief Creates a seqan2 alphabet string over a seqan3 alphabet range
 *
 * \param rng A seqan3 alphabet range.
 * \return    A seqan::String over seqan::Iupac.
 *
 * \noapi
 */
template <std::ranges::range rng_t>
auto convert_to_seqan2_iupac_alphabet(rng_t&& rng) {
    seqan::String<seqan::Iupac> v;

    auto view = rng | std::views::transform([](auto c) {
            auto t = seqan::Iupac{};
            t = c.to_char();
            return t;
        });

    if constexpr (requires() {
        { std::ranges::size(rng) } -> std::same_as<size_t>;
    }) {
        resize(v, std::ranges::size(rng), seqan::Exact());
        std::ranges::copy(view, begin(v));
    } else {
        for (auto c : view) {
            appendValue(v, c);
        }
    }
    return v;
}


/**\brief Creates a seqan2 alphabet string over a seqan3 alphabet range
 *
 * \param rng A seqan3 alphabet range.
 * \return    A seqan::String over AlphabetAdaptor<·> elements.
 *
 * \noapi
 */
template <std::ranges::range rng_t>
auto convert_to_seqan2_qualities(rng_t&& rng) {
    return convert_to_seqan2_alphabet(std::forward<rng_t>(rng));
}

/**\brief Creates a seqan2 cigar string over a seqan3 cigar range
 *
 * \param rng A seqan3 cigar range.
 * \return    A seqan::String over AlphabetAdaptor<seqan:CigarElement<>> elements.
 *
 * \noapi
 */
template <std::ranges::range rng_t>
auto convert_to_seqan2_cigar(rng_t&& rng) {
    seqan::String<seqan::CigarElement<>> v;
    auto view = rng | std::views::transform([](auto c) {
        auto t = seqan::CigarElement{};
        t.operation = get<1>(c).to_char();
        t.count     = get<0>(c);
        return t;
    });

    if constexpr (requires() {
        { std::ranges::size(rng) } -> std::same_as<size_t>;
    }) {
        resize(v, std::ranges::size(rng), seqan::Exact());
        std::ranges::copy(view, begin(v));
    } else {
        for (auto c : view) {
            appendValue(v, c);
        }
    }
    return v;
}

/**\brief Creates a seqan2 string over any range with elements of type 'char'
 *
 * \param rng Any range with values of type 'char'
 * \return    A seqan::String over 'char'.
 *
 * \noapi
 */
template <std::ranges::range rng_t>
auto convert_to_seqan2_string(rng_t&& rng) {
    seqan::String<char> v;
    if constexpr (requires() {
        { std::ranges::size(rng) } -> std::same_as<size_t>;
    }) {
        resize(v, std::ranges::size(rng));
        std::ranges::copy(rng, begin(v));
    } else {
        for (auto c : rng) {
            appendValue(v, c);
        }
    }
    return v;
}

}

namespace io2 {

/**\brief A helper variable to declare types
 *
 * \tparam t the type being deduced.
 *
 * for example:
 * ```
 * template <typename t>
 * struct A {
 *     [[no_unique_address]] io2::detail::empty_class<t> value_type;
 * };
 *
 * // This allows to declare `t` in two ways:
 * A<int> a{};
 *
 * // Or the fancy way:
 * A a {
 *  .value_type = io2::type<int>;
 * }
 * ```
 */
template <typename t>
inline constexpr auto type = detail::empty_class<t>{};



bool validExtension(std::filesystem::path const& _path, std::vector<std::string> const& _list) {
    for (auto e : _list) {
        auto s = _path.filename().string();
        if (s.size() < e.size()) continue;
        s = s.substr(s.size() - e.size());
        if (s == e) {
            return true;
        }
    }
    return false;
}

}
