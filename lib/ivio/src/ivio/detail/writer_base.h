// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <memory>
#include <optional>

namespace ivio {

template <typename writer>
struct writer_base {
protected:
    struct pimpl;
    std::unique_ptr<pimpl> pimpl_;
public:
    writer_base(std::unique_ptr<pimpl> pimpl_)
        : pimpl_{std::move(pimpl_)}
    {}
};

}
