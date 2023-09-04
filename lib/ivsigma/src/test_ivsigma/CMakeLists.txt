# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.8)

project(test_ivsigma)

add_executable(${PROJECT_NAME}
    main.cpp
)
target_link_libraries(${PROJECT_NAME}
    ivsigma::ivsigma
    fmt
)
target_compile_options(${PROJECT_NAME} PUBLIC
    -U NDEBUG
)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
