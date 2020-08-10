# This is the CMakeLists file for the SolutioLabs repository.
cmake_minimum_required(VERSION 3.0.0)
if(COMMAND CMAKE_POLICY)
  cmake_policy(SET CMP0003 NEW)
endif()

project(SolutioLabs)

# DCMTK - needed by SolutioCpp
find_package(DCMTK REQUIRED)

# Set compiler flags
set(GCC_COVERAGE_COMPILE_FLAGS "-fopenmp -std=c++17 -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(ct_hu ${CMAKE_CURRENT_SOURCE_DIR}/ct_hu.cpp)
target_link_libraries(ct_hu solutio ${DCMTK_LIBRARIES})

add_executable(gamma_1d ${CMAKE_CURRENT_SOURCE_DIR}/gamma_1d.cpp)
target_link_libraries(gamma_1d solutio ${DCMTK_LIBRARIES})
