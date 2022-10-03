add_executable(simplex-test simplex-test.cpp)
target_link_libraries(simplex-test PUBLIC simplex)

add_executable(process-file process-file.cpp)
target_link_libraries(process-file PUBLIC simplex)
