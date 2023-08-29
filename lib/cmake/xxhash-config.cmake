# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.16)
if (TARGET xxHash::xxhash)
    return()
endif()

option(BUILD_SHARED_LIBS "Build shared libs" OFF) #optional
set(XXHASH_BUILD_ENABLE_INLINE_API TRUE) #optional
set(XXHASH_BUILD_XXHSUM OFF) #optional
add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/../xxHash/cmake_unofficial)
