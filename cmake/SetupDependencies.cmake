
include(cmake/SetupPaths.cmake)
include(cmake/SetupStdFilesystem.cmake)
include(cmake/SetupDependenciesFind.cmake)
include(cmake/SetupDependenciesCMake.cmake)
include(cmake/SetupDependenciesConan.cmake)


##################################################################
### Link all the things!                                       ###
##################################################################
add_library(gl-deps INTERFACE)

if(GL_ENABLE_OPENMP)
    target_link_libraries(gl-flags INTERFACE OpenMP::OpenMP_CXX)
endif()
if(GL_ENABLE_MPI)
    target_link_libraries(gl-deps INTERFACE MPI::MPI_CXX)
endif()
target_link_libraries(gl-deps INTERFACE h5pp::h5pp)
target_link_libraries(gl-flags INTERFACE Threads::Threads)

if(GL_ENABLE_OPENGL)
    target_link_libraries(gl-deps INTERFACE gl-opengl)
endif()
