#pragma once

#include <optional>

namespace io2::detail {
/** Iterator that enables iterating over the reader
 */
template <typename Reader, typename RecordView, typename Record=RecordView>
struct iterator {
    using difference_type = std::ptrdiff_t;
    using value_type      = Record;
    using reference_type  = RecordView&;
    using rvalue_reference_type = RecordView&&;
    Reader* reader;
    RecordView const* nextItem = reader?reader->next():nullptr;

    auto operator*() const -> RecordView const& {
       return *nextItem;
    }
    auto operator++() -> iterator& {
        nextItem = reader->next();
        return *this;
    }
    auto operator++(int) -> iterator& {
        nextItem = reader->next();
        return *this;
    }

    auto operator!=(iterator const& iter) const {
        return nextItem != iter.nextItem;
    }
    auto operator==(iterator const& iter) const {
        return nextItem == iter.nextItem;
    }

};

}
