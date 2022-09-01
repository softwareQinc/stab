include(GoogleTest)
set(TARGET_NAME "unit_tests")

find_package(qpp REQUIRED)

#### Link Google Test dynamically if using MSVC
if (MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif ()

add_subdirectory(lib/googletest-release-1.12.1)

aux_source_directory(tests TEST_FILES)

add_executable(${TARGET_NAME} tests/main.cpp)
add_dependencies(unit_tests ${TARGET_NAME})

#### Build all tests in ${TEST_FILES}
foreach (FILE ${TEST_FILES})
    target_sources(${TARGET_NAME} PUBLIC ${FILE})
endforeach ()

target_link_libraries(${TARGET_NAME} PUBLIC gmock stab-common qasmtools Eigen3::Eigen ${QPP_LINK_DEPS} libqpp)

gtest_discover_tests(${TARGET_NAME})

#### If on Windows, copy the stab-common.dll into the unit_tests folder
if (WIN32)
    if (MSVC)
        add_custom_command(TARGET unit_tests POST_BUILD        # execute after the linking stage
                COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
                "${PROJECT_BINARY_DIR}/stab-common.dll"        # <--this is in-file
                $<TARGET_FILE_DIR:unit_tests>)                 # <--this is out-file path
    endif ()
endif ()