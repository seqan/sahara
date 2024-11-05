# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.16)

project (sahara VERSION 1.0.0)

# Force language standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


find_package(OpenMP QUIET)
include(cmake/CPM.cmake)
CPMAddPackage("gh:SGSSGene/cpm.dependencies@1.0.0")
CPMLoadDependenciesFile("${CMAKE_CURRENT_SOURCE_DIR}/cpm.dependencies")

add_subdirectory(src/sahara)
