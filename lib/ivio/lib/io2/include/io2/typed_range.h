#pragma once

#include "any_iter.h"

#include <ranges>

namespace io2 {

template <typename T>
struct typed_range {
    std::any range;

    detail::any_iter<T> begin_;
    detail::any_iter<T> end_;

    typed_range() = default;

    typed_range(typed_range const& _other)
        : begin_{_other.begin_}
        , end_{_other.end_}
    {}

    typed_range(typed_range&& _other)
        : range{std::move(_other.range)}
        , begin_{std::move(_other.begin_)}
        , end_{std::move(_other.end_)}
    {}
    template <typename rng_t>
    typed_range(rng_t const& rng)
        : begin_{std::ranges::begin(rng)}
        , end_{std::ranges::end(rng)}
    {}

    template <typename rng_t>
    typed_range(rng_t&& rng)
        : range{std::make_any<rng_t>(std::move(rng))}
        , begin_{std::ranges::begin(get<rng_t>())}
        , end_{std::ranges::end(get<rng_t>())}
    {}

    auto operator=(typed_range const& _other) -> typed_range& = delete;
/*        begin_ = _other.begin_;
        end_   = _other.end_;
        return *this;
    }*/

    auto operator=(typed_range&& _other) -> typed_range& {
        range  = std::move(_other.range);
        begin_ = std::move(_other.begin_);
        end_   = std::move(_other.end_);
        return *this;
    }
    template <typename rng_t>
    auto operator=(rng_t const& rng) -> typed_range& {
        begin_ = std::ranges::begin(rng);
        end_ = std::ranges::end(rng);
        return *this;
    }

    template <typename rng_t>
    auto operator=(rng_t&& rng) -> typed_range& {
        range  = std::make_any<rng_t>(std::move(rng));
        begin_ = std::ranges::begin(get<rng_t>());
        end_   = std::ranges::end(get<rng_t>());
        return *this;
    }


    template <typename rng_t>
    auto get() -> auto& {
        return any_cast<rng_t&>(range);
    }

    auto begin() const {
        return begin_;
    }
    auto end() const {
        return end_;
    }

    friend auto begin(typed_range const& tr) {
        return tr.begin_;
    };

    friend auto end(typed_range const& tr) {
        return tr.end_;
    }
};

template <typename T>
struct sized_typed_range : typed_range<T> {
    size_t size_{};

    sized_typed_range() = default;

    sized_typed_range(sized_typed_range const& _other)
        : typed_range<T>{static_cast<typed_range<T> const&>(_other)}
        , size_{_other.size_}
    {}

    sized_typed_range(sized_typed_range&& _other)
        : typed_range<T>{static_cast<typed_range<T>&&>(_other)}
        , size_{_other.size_}
    {}

    template <typename rng_t>
    sized_typed_range(rng_t&& rng)
        : typed_range<T>{std::move(rng)}
        , size_{std::ranges::size(typed_range<T>::template get<rng_t>())}
    {}

    template <typename rng_t>
    sized_typed_range(rng_t const& rng)
        : typed_range<T>{rng}
        , size_{std::ranges::size(rng)}
    {}

    auto operator=(sized_typed_range const& _other) -> sized_typed_range& = delete;

    auto operator=(sized_typed_range&& _other) -> sized_typed_range& {
        typed_range<T>::operator=(std::move(_other));
        size_ = _other.size_;
        return *this;
    }
    template <typename rng_t>
    auto operator=(rng_t const& rng) -> sized_typed_range& {
        typed_range<T>::operator=(rng);
        size_ = std::ranges::size(rng);
        return *this;
    }

    template <typename rng_t>
    auto operator=(rng_t&& rng) -> sized_typed_range& {
        typed_range<T>::operator=(std::move(rng));
        size_ = std::ranges::size(typed_range<T>::template get<rng_t>());
        return *this;
    }


    auto size() const {
        return size_;
    }
};


}
