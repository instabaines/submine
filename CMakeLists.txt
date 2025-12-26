message(STATUS "ROOT CMakeLists executing: installing Python package files")


cmake_minimum_required(VERSION 3.15)
project(submine LANGUAGES CXX)
# add_subdirectory(submine/backends/pybind11 EXCLUDE_FROM_ALL)

# addpybind
include(FetchContent)
FetchContent_Declare(pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.11.1
)
FetchContent_MakeAvailable(pybind11)

# Build the C++ extensions
add_subdirectory(submine/backends/sopagrami/cpp)
add_subdirectory(submine/backends/gspan)

# Install the *Python* package directory into the wheel.
# IMPORTANT: Source has NO trailing slash so the "submine" folder is preserved.
# DESTINATION "." means it lands at wheel root: submine/__init__.py, submine/core/..., etc.
install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/submine"
  DESTINATION "."
)
