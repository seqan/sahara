# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------

cmake_minimum_required (VERSION 3.12)
if (TARGET ivsigma::ivsigma)
    return()
endif()

add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/src/ivsigma ${CMAKE_CURRENT_BINARY_DIR}/ivsigma)
