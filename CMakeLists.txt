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


#### stab-common library
if (WIN32)
    if (MSVC)
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    endif ()
endif ()
add_library(stab-common SHARED ${SOURCES})
target_link_libraries(stab-common Eigen3::Eigen)

add_executable(stab main.cpp)
target_link_libraries(stab stab-common qasmtools ${QPP_LINK_DEPS} libqpp)

#### Unit testing
add_subdirectory(${CMAKE_SOURCE_DIR}/unit_tests)
