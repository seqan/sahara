#pragma once

#include <any>
#include <functional>
#include <ranges>

namespace io2::detail {

template <typename T>
struct any_iter {
    using difference_type = std::ptrdiff_t;
    using value_type      = T;
    using reference_type  = T&;
    using rvalue_reference_type = T&&;

    std::any                                              iter;
    std::function<T(std::any const&)>                     dref;
    std::function<void(std::any&)>                        inc;
    std::function<bool(std::any const&, std::any const&)> equal;
    std::function<std::any(std::any const&)>              copy = [](std::any const&) { return std::any{}; };

    any_iter() = default;
    any_iter(any_iter const& _other)
        : iter{_other.copy(_other.iter)}
        , dref{_other.dref}
        , inc{_other.inc}
        , equal{_other.equal}
        , copy{_other.copy}
    {}
    any_iter(any_iter&&) = default;


    template <typename I>
    any_iter(I i) {
        iter = i;
        dref = [](std::any const& i) -> T {
            return *(std::any_cast<I const&>(i));
        };
        inc = [](std::any& i) {
            ++std::any_cast<I&>(i);
        };
        equal = [](std::any const& i, std::any const& _other) -> bool {
            return std::any_cast<I const&>(i) == std::any_cast<I const&>(_other);
        };
        copy = [](std::any const& i) -> std::any {
            return std::any_cast<I const&>(i);
        };
    }

    auto operator=(any_iter const& _other) -> any_iter& {
        iter  = _other.copy(_other.iter);
        dref  = _other.dref;
        inc   = _other.inc;
        equal = _other.equal;
        copy  = _other.copy;
        return *this;
    }
    auto operator=(any_iter&&) -> any_iter& = default;

    auto operator*() const {
        return dref(iter);
    }
    auto operator++() -> any_iter& {
        inc(iter);
        return *this;
    }
    auto operator++(int) -> any_iter& {
        inc(iter);
        return *this;
    }

    auto operator!=(any_iter const& _other) const {
        if (!equal and !_other.equal) return false;
        if (equal and !_other.equal) return false;
        if (!equal and _other.equal) return false;
        return !equal(iter, _other.iter);
    }
    auto operator==(any_iter const& _other) const {
        return !(*this != _other);
    }
};

}
