message(STATUS "ROOT CMakeLists executing: installing Python package files")

cmake_minimum_required(VERSION 3.15)
project(submine LANGUAGES CXX)

# Build the C++ extension
add_subdirectory(submine/thirdparty/sopagrami/cpp)

# Install the *Python* package directory into the wheel.
# IMPORTANT: Source has NO trailing slash so the "submine" folder is preserved.
# DESTINATION "." means it lands at wheel root: submine/__init__.py, submine/core/..., etc.
install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/submine"
  DESTINATION "."
  FILES_MATCHING
    PATTERN "*.py"
)
