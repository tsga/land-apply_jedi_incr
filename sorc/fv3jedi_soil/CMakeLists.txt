################################################################################
# fv3jedi_soil_incr
################################################################################

list( APPEND sorc_files_fv3jedisoilinc
    ./sorc/fv3jedi_soil/fv3jedi_add_land_increment.cc
    ./sorc/fv3jedi_soil/fv3jedi_add_land_increment.h
    ./sorc/fv3jedi_soil/soil_increments_cpp_interface.h
    ./sorc/fv3jedi_soil/soil_increments_c.f90
)

ecbuild_add_executable( TARGET          fv3jedi_soil_incr.x
                        SOURCES         ${sorc_files_fv3jedisoilinc} )
target_compile_features( fv3jedi_soil_incr.x PUBLIC cxx_std_17)
target_link_libraries( fv3jedi_soil_incr.x PUBLIC NetCDF::NetCDF_CXX NetCDF::NetCDF_Fortran oops atlas fv3jedi)

################################################################################
