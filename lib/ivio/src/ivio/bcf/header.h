// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <vector>
#include <string>

namespace ivio::bcf {

struct header {
    std::vector<std::tuple<std::string, std::string>> table;
    std::vector<std::string> genotypes;
};

}
