cmake_minimum_required(VERSION 3.15)
project(SopaGrami LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(alg STATIC src/alg.cpp)
target_include_directories(alg PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_target_properties(alg PROPERTIES POSITION_INDEPENDENT_CODE ON)


# add_subdirectory(
#   ${CMAKE_CURRENT_SOURCE_DIR}/../pybind11
#   ${CMAKE_CURRENT_BINARY_DIR}/pybind11_build
# )
# ---- pybind11 (reuse if provided by top-level; otherwise vendor locally) ----
if (NOT COMMAND pybind11_add_module)
  add_subdirectory(
    ${CMAKE_CURRENT_LIST_DIR}/../../../pybind11
    ${CMAKE_CURRENT_BINARY_DIR}/pybind11_build
  )
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
