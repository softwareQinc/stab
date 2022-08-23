cmake_minimum_required(VERSION 3.12)
project(stab)

set(CMAKE_CXX_STANDARD 17)
find_package(Eigen3 REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/src)
file(GLOB SOURCES CONFIGURE_DEPENDS "src/*.cpp")

add_executable(stab main.cpp ${SOURCES})
target_link_libraries(stab Eigen3::Eigen)

#### qasmtools library
target_include_directories(stab PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/qasmtools/include/>)