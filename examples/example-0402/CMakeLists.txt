set(ExampleName "example" ${IRON_MESH_READER_ROOT}/src/meshReader.F90)

cmake_minimum_required(VERSION 3.4 FATAL_ERROR)

project(${ExampleName} VERSION 1.0 LANGUAGES C Fortran)

list(APPEND CMAKE_PREFIX_PATH "${OPENCMISS_INSTALL_DIR}")

find_package(OpenCMISSLibs 1.3.0 REQUIRED COMPONENTS Iron CONFIG)

add_subdirectory(src/iron)

