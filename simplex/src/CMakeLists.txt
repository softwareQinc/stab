add_library(simplex
    simplex.cpp
    A_matrix.cpp
    Q_matrix.cpp
    parse-stim.cpp)

target_include_directories(simplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set_property(TARGET simplex PROPERTY POSITION_INDEPENDENT_CODE ON)

option(SIMPLEX_DENSE "Use dense-matrix implementation" OFF)
if(${SIMPLEX_DENSE})
    add_compile_definitions(SIMPLEX_DENSE)
endif()
