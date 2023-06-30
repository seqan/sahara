cmake_minimum_required (VERSION 3.14)

function(test_snippet file_base_name)
    add_executable(snippet_${file_base_name} ${file_base_name})
    target_link_libraries(snippet_${file_base_name} ivsigma::ivsigma fmt)

    add_test (NAME test_snippet_${file_base_name}
              COMMAND bash -c "$<TARGET_FILE:snippet_${file_base_name}> | diff - ${CMAKE_CURRENT_SOURCE_DIR}/${file_base_name}.out")
endfunction()

test_snippet("char_to_rank.cpp")
test_snippet("compact_encoding.cpp")
test_snippet("complement.cpp")
test_snippet("fasta_reader_example.cpp")
test_snippet("normalize_char.cpp")
test_snippet("rank_to_char.cpp")
test_snippet("reverse_complement.cpp")
test_snippet("verify.cpp")
test_snippet("winnowing_minimizers.cpp")
