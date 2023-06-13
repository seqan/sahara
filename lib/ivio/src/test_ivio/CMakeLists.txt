# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.8)

project(test_ivio)

# fmindex-collectionunittests
add_executable(${PROJECT_NAME}
    main.cpp
    fasta_reader.cpp
    fasta_writer.cpp
    sam_reader.cpp
    sam_writer.cpp
    vcf_reader.cpp
    vcf_writer.cpp
)
target_link_libraries(${PROJECT_NAME}
    Catch2::Catch2
    ivio::ivio
)


add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
