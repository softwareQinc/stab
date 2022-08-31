cmake_minimum_required(VERSION 3.12)
project(stab)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 REQUIRED)
find_package(qpp REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/src)
file(GLOB SOURCES CONFIGURE_DEPENDS "src/*.cpp")

#### qasmtools_stab library
add_library(libqasmtools INTERFACE)
target_include_directories(libqasmtools INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qasmtools/include/>)

add_executable(stab main.cpp ${SOURCES})
target_link_libraries(stab Eigen3::Eigen libqasmtools)
target_link_libraries(stab ${QPP_LINK_DEPS} libqpp)

#### Unit testing
add_subdirectory(${CMAKE_SOURCE_DIR}/unit_tests)
