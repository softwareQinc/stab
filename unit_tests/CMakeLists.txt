include(GoogleTest)
set(TARGET_NAME "unit_tests")

#### Link Google Test dynamically if using MSVC
if (MSVC)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif ()

add_subdirectory(libs/googletest-release-1.12.1 EXCLUDE_FROM_ALL SYSTEM)

aux_source_directory(tests TEST_FILES)

add_executable(${TARGET_NAME} EXCLUDE_FROM_ALL tests/main.cpp)
add_dependencies(unit_tests ${TARGET_NAME})

#### Build all tests in ${TEST_FILES}
foreach (FILE ${TEST_FILES})
    target_sources(${TARGET_NAME} PUBLIC ${FILE})
endforeach ()

target_link_libraries(${TARGET_NAME} gmock stab-common qasmtools)
if (${USE_QPP})
    target_link_libraries(${TARGET_NAME} ${QPP_LINK_DEPS} libqpp)
endif ()

gtest_discover_tests(${TARGET_NAME})
