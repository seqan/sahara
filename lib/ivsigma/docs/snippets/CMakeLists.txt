cmake_minimum_required (VERSION 3.14)

function(test_snippet file_base_name)
    add_executable(snippet_${file_base_name} ${file_base_name})
    target_link_libraries(snippet_${file_base_name} ivsigma::ivsigma fmt)

    add_test (NAME test_snippet_${file_base_name}
              COMMAND bash -c "$<TARGET_FILE:snippet_${file_base_name}> | diff - ${CMAKE_CURRENT_SOURCE_DIR}/${file_base_name}.out")
endfunction()

test_snippet("char_to_rank.cpp")
test_snippet("fasta_reader_example.cpp")
