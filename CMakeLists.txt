cmake_minimum_required(VERSION 3.14)
project(stab)
enable_testing()

set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 3.4 REQUIRED) # because we're using Eigen::seq()
find_package(qpp REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/src)
file(GLOB SOURCES CONFIGURE_DEPENDS "src/*.cpp")

#### qasmtools_stab library
add_library(qasmtools INTERFACE)
target_include_directories(qasmtools INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qasmtools/include/>)

#### Export stab-common library symbols on Windows
if (WIN32)
    if (MSVC)
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    endif ()
endif ()
add_library(stab-common SHARED ${SOURCES} "benchmarking/stab_tests.cpp")
target_link_libraries(stab-common qasmtools Eigen3::Eigen ${QPP_LINK_DEPS} libqpp)

add_executable(stab main.cpp "benchmarking/stab_tests.cpp")
target_link_libraries(stab stab-common qasmtools ${QPP_LINK_DEPS} libqpp)

#### Unit testing
add_subdirectory(${CMAKE_SOURCE_DIR}/unit_tests)
