// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <cassert>
#include <memory>
#include <optional>

namespace ivio {

template <typename reader>
struct reader_base {
    struct iter {
        using record_view = typename reader::record_view;
        using difference_type = std::ptrdiff_t;
        using value_type = typename reader::record;

        reader* reader_{};
        std::optional<record_view> nextItem = [this]() -> std::optional<record_view> { if (reader_) return reader_->next(); return std::nullopt; }();

        auto operator*() const -> record_view {
           return *nextItem;
        }
        auto operator->() const -> record_view* {
            return &*nextItem;
        }
        auto operator++() -> iter& {
            nextItem = reader_->next();
            return *this;
        }
        auto operator++(int) {
            auto currentItem = std::optional<value_type>(*nextItem);
            nextItem = reader_->next();
            return currentItem;
        }
        auto operator==(iter const& _other) const {
            assert(reader_);
            assert(_other.reader_ == nullptr);
            return !nextItem.has_value();
        }
        //!WORKAROUND clang15 is requiring this
        auto operator!=(iter const& _other) const {
            return !(*this == _other);
        }

    };

    struct pimpl; //!WORKAROUND, this should be protected, but clang15 fails, see readLine
protected:

    std::unique_ptr<pimpl> pimpl_;
public:
    reader_base(std::unique_ptr<pimpl> pimpl_)
        : pimpl_{std::move(pimpl_)}
    {}
    friend auto begin(reader& reader_) {
        return iter{&reader_};
    }
    friend auto end(reader&) {
        return iter{nullptr};
    }
};

}
