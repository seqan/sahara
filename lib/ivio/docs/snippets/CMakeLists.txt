cmake_minimum_required (VERSION 3.14)

function(test_snippet file_base_name params)
    add_executable(snippet_${file_base_name} ${file_base_name})
    target_link_libraries(snippet_${file_base_name} ivio::ivio fmt)

    add_test (NAME test_snippet_${file_base_name}
              COMMAND bash -c "$<TARGET_FILE:snippet_${file_base_name}> ${params} | diff - ${CMAKE_CURRENT_SOURCE_DIR}/${file_base_name}.out")
endfunction()

test_snippet("read_vcf.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/data/vcf.vcf")
test_snippet("fasta_example_01.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/data/demo.fa")
