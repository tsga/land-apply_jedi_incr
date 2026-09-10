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
./NoahMPdisag_module.f90
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
set ( snow_exe_files ./apply_incr_noahmp_snow.f90 )

ecbuild_add_executable( TARGET   apply_snow_incr.exe
                        SOURCES  ${snow_exe_files}
                        LIBS     jediincr_snow
                        LINKER_LANGUAGE ${JEDIINCR_LINKER_LANGUAGE}
                       )

################################################################################
