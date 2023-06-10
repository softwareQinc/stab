cmake_minimum_required(VERSION 3.15)
project(stab)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

#### Eigen3
message(STATUS "Detecting Eigen3")
#### Location specified via an environment variable
set(LOCATION_SET_VIA_ENV FALSE)
if (DEFINED ENV{EIGEN3_INSTALL_DIR})
    set(EIGEN3_INSTALL_DIR_ENV $ENV{EIGEN3_INSTALL_DIR})
    set(LOCATION_SET_VIA_ENV TRUE)
endif ()
#### Location set via CMake variable, trumps all other settings
set(EIGEN3_INSTALL_DIR "" CACHE PATH "Path to Eigen3")
if (NOT ${EIGEN3_INSTALL_DIR} STREQUAL "")
    message(STATUS "Overriding automatic Eigen3 detection (EIGEN3_INSTALL_DIR CMake variable)")
    if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR})
        message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INSTALL_DIR})")
        include_directories(SYSTEM "${EIGEN3_INSTALL_DIR}")
    else ()
        message(FATAL_ERROR "Invalid path to Eigen3 installation")
    endif ()
    #### Location set via environment variable
elseif (LOCATION_SET_VIA_ENV)
    message(STATUS "Overriding automatic Eigen3 detection (EIGEN3_INSTALL_DIR environment variable)")
    if (IS_DIRECTORY ${EIGEN3_INSTALL_DIR_ENV})
        message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INSTALL_DIR_ENV})")
        include_directories(SYSTEM "${EIGEN3_INSTALL_DIR_ENV}")
    else ()
        message(FATAL_ERROR "Invalid path to Eigen3 installation")
    endif ()
    #### Try to find the location automatically
else ()
    find_package(Eigen3 3.0 QUIET NO_MODULE)
    if (NOT TARGET Eigen3::Eigen) # did not find Eigen3 automatically
        message(FATAL_ERROR "Eigen3 not detected! Please point EIGEN3_INSTALL_DIR to your Eigen3 location when building with CMake, for example
    cmake --build build -DEIGEN3_INSTALL_DIR=$HOME/eigen3
or set the EIGEN3_INSTALL_DIR environment variable to point to your Eigen3 installation, for example (UNIX/Linux)
    export EIGEN3_INSTALL_DIR=$HOME/eigen3")
        #### Eigen3
        message(STATUS "Detecting Eigen3")
    endif ()
    message(STATUS "Detecting Eigen3 - done (in ${EIGEN3_INCLUDE_DIR})")
    # Eigen3 header-only dependencies to be injected in the main CMakeLists.txt
    set(STAB_EIGEN3_LINK_DEPS Eigen3::Eigen)
endif ()

# Enable Quantum++ support
option(USE_QPP "Enable Quantum++ support" OFF)
if (${USE_QPP})
    find_package(qpp REQUIRED)
    add_compile_definitions(USE_QPP)
endif ()

include_directories(${PROJECT_SOURCE_DIR}/src)
file(GLOB SOURCES CONFIGURE_DEPENDS "src/*.cpp")

#include_directories(${PROJECT_SOURCE_DIR}/simplex/src/include)

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

# add_library(stab-common SHARED ${SOURCES} "benchmarking/stab_tests.cpp")
add_library(stab-common SHARED ${SOURCES})

if (${USE_QPP})
    target_link_libraries(stab-common ${QPP_LINK_DEPS} libqpp)
else ()
    target_link_libraries(stab-common ${STAB_EIGEN3_LINK_DEPS})
endif ()

#### Maximum number of qubits allowed in qpp::AffineState::to_ket()
set(MAX_QUBITS_STATE_VECTOR "16" CACHE STRING "Maximum number of qubits allowed in AffineState::to_ket()")
add_compile_definitions(MAX_QUBITS_STATE_VECTOR=${MAX_QUBITS_STATE_VECTOR})

add_executable(stab stab/main.cpp)
target_link_libraries(stab stab-common qasmtools)
if (${USE_QPP})
    target_link_libraries(stab ${QPP_LINK_DEPS} libqpp)
endif ()

# add_executable(bench_stab benchmarking/stab_tests.cpp)
# target_link_libraries(bench_stab stab-common)

#### Unit testing
add_subdirectory(${CMAKE_SOURCE_DIR}/unit_tests)
