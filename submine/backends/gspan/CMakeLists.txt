cmake_minimum_required(VERSION 3.16)
project(gspan_cpp_bindings LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(FetchContent)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

# ---- pybind11 (vendored) ----

include(FetchContent)

# Try system / vcpkg / conda / installed pybind11 first
find_package(pybind11 CONFIG QUIET)

if (NOT pybind11_FOUND)
  message(STATUS "pybind11 not found via find_package; fetching with FetchContent...")

  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG        v2.12.0
  )
  FetchContent_MakeAvailable(pybind11)
endif()



set(GSPAN_SOURCES
  dfs.cpp
  graph.cpp
  gspan.cpp
  ismin.cpp
  misc.cpp
)

add_executable(gspan_cli
  main.cpp
  ${GSPAN_SOURCES}
)
target_include_directories(gspan_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

pybind11_add_module(gspan_cpp
  python_bindings.cpp
  ${GSPAN_SOURCES}
)
target_include_directories(gspan_cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(gspan_cpp PROPERTIES PREFIX "")

# Warnings
if (MSVC)
  target_compile_options(gspan_cli PRIVATE /W4)
  target_compile_options(gspan_cpp PRIVATE /W4)
else()
  target_compile_options(gspan_cli PRIVATE -Wall -Wextra -Wformat)
  target_compile_options(gspan_cpp PRIVATE -Wall -Wextra -Wformat)
endif()

# ---- Install into the wheel ----
install(TARGETS gspan_cpp
  LIBRARY DESTINATION "submine/algorithms"
  RUNTIME DESTINATION "submine/algorithms"
  ARCHIVE DESTINATION "submine/algorithms"
)
