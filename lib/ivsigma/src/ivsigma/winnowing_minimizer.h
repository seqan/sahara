#pragma once

#include "compact_encoding.h"

#include <deque>

namespace ivs {

template <alphabet_c Alphabet, bool DuplicatesAllowed=true, bool UseCanonicalKmers=true>
struct winnowing_minimizer {
    using Encoding = compact_encoding<Alphabet, UseCanonicalKmers>;

    Encoding hash;
    size_t   window{};

    winnowing_minimizer(std::span<uint8_t const> _values, size_t _k, size_t _window, size_t _seed = 0)
        : hash{_values, _k, _seed}
        , window{_window} {
    }
    auto size() const -> size_t {
        if (hash.size() < window) return 0;
        return hash.size() - window + 1;
    }

    struct iterator {
        winnowing_minimizer const* ptr;

        Encoding::iterator                     iter;
        std::deque<std::pair<size_t, size_t>>  values{};
        size_t                                 pos{};

        iterator(winnowing_minimizer const& minimizer)
            : ptr{&minimizer}
            , iter{begin(ptr->hash)}
        {
            if (ptr->size() == 0) return;

            values.emplace_back(0, *iter);
            while (pos+1 < ptr->window) {
                ++iter;
                ++pos;
                // pop at the end, until element smaller is found
                while (!values.empty() && values.back().second > *iter) {
                    values.pop_back();
                }
                values.emplace_back(pos, *iter);
            }
        }


        auto operator*() const -> size_t {
            return values.front().second;
        }

        auto operator++() -> iterator& {
            auto lastReportedPos  = values.front().first;
            auto lastReportedHash = values.front().second;
            auto abortPred = [&]() {
                if constexpr (DuplicatesAllowed) {
                    return lastReportedPos == values.front().first;
                } else {
                    return lastReportedHash == values.front().second;
                }
            };

            for (; abortPred(); ) {
                ++iter;
                ++pos;
                if (iter == end(ptr->hash)) {
                    values.clear();
                    return *this;
                }

                // drop at the beginning if outside of the window
                auto const& [_pos, _hash] = values.front();
                if (_pos + ptr->window <= pos) {
                    values.pop_front();
                }

                // pop at the end, until element smaller is found
                while (!values.empty() && values.back().second >= *iter) {
                    values.pop_back();
                }

                values.emplace_back(pos, *iter);
                if (ptr->window == 1) break;
            }
            return *this;
        }


        bool operator==(nullptr_t) const {
            return values.empty();
        }
    };

    friend auto begin(winnowing_minimizer const& index) -> iterator {
        return iterator{index};
    }
    friend auto end(winnowing_minimizer const&) -> std::nullptr_t {
        return nullptr;
    }
};

}
