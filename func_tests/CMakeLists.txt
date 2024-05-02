project(FuncTests)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/runtests.sh
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_1.fq.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_2.fq.gz
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_w.fq
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_x.fq
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_y.fq
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_1_L.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/func_tests/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_1_V.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/func_tests/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_1_R.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/func_tests/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_2_L.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/func_tests/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_2_V.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/func_tests/expected_output)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/expected_output/dummy_2_R.fa
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/func_tests/expected_output)
