cmake_minimum_required( VERSION 3.12)
project(landincr LANGUAGES Fortran C CXX )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_FORTRAN_STANDARD 08)
set(CMAKE_FORTRAN_STANDARD_REQUIRED ON)
set(CMAKE_FORTRAN_EXTENSIONS OFF)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU)$")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-0")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

find_package(NetCDF REQUIRED COMPONENTS CXX Fortran)
find_package(oops CONFIG REQUIRED HINTS "../../build" "../../install")
find_package(fv3jedi CONFIG REQUIRED HINTS "../../build" "../../install")
find_package(MPI REQUIRED Fortran)

add_subdirectory( sorc/snow )
add_subdirectory( sorc/fv3jedi_soil )
