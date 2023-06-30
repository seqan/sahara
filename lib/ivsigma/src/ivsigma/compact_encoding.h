#pragma once

#include "concepts.h"

#include <cstdint>
#include <cmath>
#include <span>

namespace ivs::detail {

/**
 * Computes a integer processing pow function
 */
template <size_t Sigma>
auto myPow(size_t exp) -> size_t {
    if (exp == 0) return 1;
    if (exp == 1) return Sigma;
    if (Sigma == 0) return 0;
    if (Sigma == 1) return 1;
    auto v = myPow<Sigma>(exp/2);
    v = v*v;
    if (exp % 2 == 1) {
        v = v * Sigma;
    }
    return v;
}

/**
 * computes a compact encoding, limited to 64bit
 * - call 'nextRight(,)' to receive the next hash value.
 * - or call 'nextLeft(,)' to receive the next hash value.
 * for initialization call 'next*(0, x)'
 */
template <size_t Sigma>
struct compact_encoding_gadget {
    size_t const k;
    size_t const maxExp{myPow<Sigma>(k-1)};
    size_t hash{};

    size_t value() const {
//        return hash ^ 0x8F3F73B5CF1C9ADE;
        return hash;
    }
    void nextRight(size_t removeValue, size_t insertValue) {
        hash = (hash - removeValue*maxExp)*Sigma + insertValue;
    }
    void nextLeft(size_t removeValue, size_t insertValue) {
        hash = (hash - removeValue)/Sigma + insertValue*maxExp;
    }
};

template <alphabet_c Alphabet>
struct compact_encoding {
    std::span<uint8_t const> values;
    size_t const k;
    size_t const seed;

   compact_encoding(std::span<uint8_t const> _values, size_t _k, size_t _seed = 0)
        : values{_values}
        , k{_k}
        , seed{_seed}
    {}

    auto size() const -> size_t {
        if (values.size() < k) return 0;
        return values.size() - k + 1;
    }

    struct iterator {
        compact_encoding const* ptr;

        compact_encoding_gadget<Alphabet::size()> fwdHash;
        compact_encoding_gadget<Alphabet::size()> bwdHash;
        size_t minHash{};
        size_t pos{};

        iterator(compact_encoding const& hash)
            : ptr{&hash}
            , fwdHash{ptr->k}
            , bwdHash{ptr->k}
        {
            if (ptr->size() > 0) {
                for (size_t i{0}; i < ptr->k; ++i) {
                    auto addValue = ptr->values[i];
                    fwdHash.nextRight(0, addValue);
                    if constexpr (alphabet_with_complement_c<Alphabet>) {
                        bwdHash.nextLeft(0, Alphabet::complement_rank(addValue));
                    }
                }
            }
            pos = ptr->k;
            if constexpr (alphabet_with_complement_c<Alphabet>) {
                minHash = std::min(fwdHash.value() ^ ptr->seed, bwdHash.value() ^ ptr->seed);
            } else {
                minHash = fwdHash.value() ^ ptr->seed;
            }
        }
        auto operator*() const -> size_t {
            return minHash;
        }
        auto operator++() -> iterator& {
            auto rmValue  = ptr->values[pos-ptr->k];
            auto addValue = ptr->values[pos];
            fwdHash.nextRight(rmValue, addValue);
            if constexpr (alphabet_with_complement_c<Alphabet>) {
                bwdHash.nextLeft(Alphabet::complement_rank(rmValue), Alphabet::complement_rank(addValue));
                minHash = std::min(fwdHash.value() ^ ptr->seed, bwdHash.value() ^ ptr->seed);
            } else {
                minHash = fwdHash.value() ^ ptr->seed;
            }
            pos += 1;
            return *this;
        }
        bool operator==(nullptr_t) const {
            return pos > ptr->values.size();
        }
    };

    friend auto begin(compact_encoding const& index) -> iterator {
        return iterator{index};
    }
    friend auto end(compact_encoding const&) -> std::nullptr_t {
        return nullptr;
    }
};
}

namespace ivs {

template <alphabet_c Alphabet>
using compact_encoding = detail::compact_encoding<Alphabet>;

}
