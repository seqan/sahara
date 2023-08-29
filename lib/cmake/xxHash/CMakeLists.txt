# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file.
# -----------------------------------------------------------------------------------------------------

option(BUILD_SHARED_LIBS "Build shared libs" OFF) #optional
...
set(XXHASH_BUILD_ENABLE_INLINE_API OFF) #optional
set(XXHASH_BUILD_XXHSUM OFF) #optional
add_subdirectory(</path/to/xxHash/cmake_unofficial/> </path/to/xxHash/build/> EXCLUDE_FROM_ALL)
...
target_link_libraries(MyTarget PRIVATE xxHash::xxhash)
