# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------
cmake_minimum_required (VERSION 3.14)

# check if this is the main project
if (NOT DEFINED PROJECT_NAME)
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
    if(WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /EHsc")
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -pedantic")
    endif()
    set(ROOT_PROJECT TRUE)
else()
    set(ROOT_PROJECT FALSE)
endif ()


project(IVSigma)

add_subdirectory(src/ivsigma)

if (${ROOT_PROJECT})
    find_package (fmt REQUIRED)

#    add_subdirectory(lib/Catch2)
    enable_testing()


    enable_testing()
    add_subdirectory(docs/snippets)
    add_subdirectory(src/test_ivsigma)
endif()
