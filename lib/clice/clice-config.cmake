cmake_minimum_required (VERSION 3.12)
if (TARGET ivio::ivio)
    return()
endif()

add_library(clice INTERFACE)
add_library(clice::clice ALIAS clice)
target_sources(clice INTERFACE FILE_SET set0 TYPE HEADERS BASE_DIRS src)
target_include_directories(clice
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/src>
    $<INSTALL_INTERFACE:include>
)
include(GNUInstallDirs)
install (TARGETS clice)
install (DIRECTORY src/clice TYPE INCLUDE FILES_MATCHING PATTERN "*.h")
install (FILES clice-config.cmake DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/)
