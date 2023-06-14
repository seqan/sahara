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
endif ()


project(clice LANGUAGES CXX
            DESCRIPTION "clice -- command line interface argument parsing")

include(clice-config.cmake)

if (${ROOT_PROJECT})
    find_package (fmt REQUIRED)

    enable_testing()

    add_executable(clice-demo src/clice-demo/main.cpp)
    target_link_libraries(clice-demo
                            clice::clice
                            fmt::fmt)
    add_test(clice-demo clice-demo)
endif()
