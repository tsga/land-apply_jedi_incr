# (C) Copyright 2022 .
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


cmake_minimum_required( VERSION 3.12)
project( jediincr VERSION 2022.10 LANGUAGES Fortran )

# ecbuild integration
find_package(ecbuild 3.3.2 REQUIRED)
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include( jediincr_compiler_flags )

################################################################################
# Dependencies
################################################################################

find_package(OpenMP COMPONENTS C Fortran)
find_package(MPI REQUIRED COMPONENTS C Fortran)
find_package(NetCDF REQUIRED COMPONENTS Fortran )

#add_definitions ( -DCPP=0 )

set( JEDIINCR_LINKER_LANGUAGE Fortran )

# public_libs: snow and soil incr
set( PUBLIC_LIBS
    NetCDF::NetCDF_Fortran
    MPI::MPI_Fortran
    OpenMP::OpenMP_C
    OpenMP::OpenMP_Fortran
)

set( MODULE_DIR ${CMAKE_BINARY_DIR}/module )

################################################################################
# jediincr_snow
################################################################################

list ( APPEND sorc_files_snowinc
./sorc/NoahMPdisag_module.f90
)

ecbuild_add_library( TARGET   jediincr_snow
                     SOURCES  ${sorc_files_snowinc}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${JEDIINCR_LINKER_LANGUAGE}
                    )

target_link_libraries(jediincr_snow PUBLIC ${PUBLIC_LIBS})

#target_compile_options(jediincr PRIVATE  -r8 -init=huge -cpp, -DCCPP=.false. )

set_target_properties(jediincr_snow PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIR}/snow_inc)
target_include_directories(jediincr_snow PUBLIC $<BUILD_INTERFACE:${MODULE_DIR}/snow_inc>)

# jedi_incr_snow.exe
set ( snow_exe_files ./sorc/apply_incr_noahmp_snow.f90 )

ecbuild_add_executable( TARGET   apply_incr.exe
                        SOURCES  ${snow_exe_files}
                        LIBS     jediincr_snow
                        LINKER_LANGUAGE ${JEDIINCR_LINKER_LANGUAGE}
                       )

################################################################################
# jediincr_soil
################################################################################

list( APPEND sorc_files_soilinc
    ./sorc/soil_increments.f90
    ./sorc/set_soilveg_snippet.f90
    ./sorc/sflx_snippet.f90
)

ecbuild_add_library( TARGET          jediincr_soil
                     SOURCES         ${sorc_files_soilinc}
                     INSTALL_HEADERS LISTED
                     LINKER_LANGUAGE ${JEDIINCR_LINKER_LANGUAGE}
                    )

target_link_libraries(jediincr_soil PUBLIC ${PUBLIC_LIBS})

set_target_properties(jediincr_soil PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_DIR}/soil_inc)
target_include_directories(jediincr_soil PUBLIC $<BUILD_INTERFACE:${MODULE_DIR}/soil_inc>)

set( soil_exe_files ./sorc/apply_incr_noahmp_soil.f90 )

ecbuild_add_executable( TARGET          apply_soil_incr.x
                        SOURCES         ${soil_exe_files}
                        LIBS            jediincr_soil
                        LINKER_LANGUAGE ${JEDIINCR_LINKER_LANGUAGE}
                       )

################################################################################
# Finalise configuration
################################################################################

ecbuild_install_project( NAME jediincr )
ecbuild_print_summary()
