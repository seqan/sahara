# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.14)

find_package(ZLIB REQUIRED)

add_library(ivio bam/reader.cpp
                 bcf/reader.cpp
                 bcf/writer.cpp
                 fasta/reader.cpp
                 fasta/reader_mt.cpp
                 fasta/writer.cpp
                 faidx/reader.cpp
                 fastq/reader.cpp
                 sam/reader.cpp
                 sam/writer.cpp
                 vcf/reader.cpp
                 vcf/writer.cpp)

add_library(ivio::ivio ALIAS ivio)

target_include_directories(ivio
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/..>
    $<INSTALL_INTERFACE:include>
)
target_link_libraries(ivio
    ${ZLIB_LIBRARIES}
)
