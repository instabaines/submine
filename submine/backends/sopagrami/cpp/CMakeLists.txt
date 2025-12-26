cmake_minimum_required(VERSION 3.15)
project(SopaGrami LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(alg STATIC src/alg.cpp)
target_sources(alg PRIVATE
    src/alg.cpp
    src/dump.cpp
)

target_include_directories(alg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(alg PROPERTIES POSITION_INDEPENDENT_CODE ON)

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

pybind11_add_module(sopagrami_cpp src/python_bindings.cpp)

target_include_directories(sopagrami_cpp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(sopagrami_cpp PRIVATE alg)
set_target_properties(sopagrami_cpp PROPERTIES PREFIX "")

# Install into scikit-build-core wheel
install(TARGETS sopagrami_cpp
  LIBRARY DESTINATION "submine/algorithms"
  RUNTIME DESTINATION "submine/algorithms"
  ARCHIVE DESTINATION "submine/algorithms"
)
