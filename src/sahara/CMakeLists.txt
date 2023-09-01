 # -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.14)

add_executable(sahara
    AdaptiveKmerIndex.cpp
    index.cpp
    kmer-index.cpp
    kmer-search.cpp
    main.cpp
    search.cpp
    uni-index.cpp
    uni-search.cpp
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
