include(GNUInstallDirs)

include(FetchContent)
FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v2.9.1.tar.gz
  URL_HASH SHA256=c6160321dc98e6e1184cc791fbeadd2907bb4a0ce0e447f2ea4ff8ab56550913
)
FetchContent_MakeAvailable(pybind11)
pybind11_add_module(_simplex MODULE pysimplex.cpp)

target_link_libraries(_simplex PUBLIC simplex)

install(TARGETS _simplex DESTINATION .)

set(ignoreMe "${SKBUILD}")
