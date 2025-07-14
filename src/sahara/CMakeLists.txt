# SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

cmake_minimum_required (VERSION 3.14)

add_executable(sahara
    AdaptiveKmerIndex.cpp
    index.cpp
    kmer-index.cpp
    kmer-search.cpp
    main.cpp
    rbi-index.cpp
    rbi-search.cpp
    rbi-index-dna4.cpp
    rbi-search-dna4.cpp
    search.cpp
    uni-index.cpp
    uni-search.cpp
    search_scheme.cpp
    read_simulator.cpp
    columba_prepare.cpp
)

target_link_libraries(sahara
    fmindex-collection::fmindex-collection
    fmt::fmt
    ivio::ivio
    ivsigma::ivsigma
    clice::clice
    cereal::cereal
    xxhash
)

set_property(TARGET sahara PROPERTY CXX_STANDARD 20)
set_property(TARGET sahara PROPERTY CXX_EXTENSIONS FALSE)
