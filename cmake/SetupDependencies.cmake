if(NOT TARGET flags)
    add_library(flags INTERFACE)
endif()

# Find and setup threads (used by OpenMP)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
target_link_libraries(Threads::Threads INTERFACE rt dl)

# Find OpenMP
find_package(OpenMP COMPONENTS CXX REQUIRED)
target_link_libraries(flags INTERFACE OpenMP::OpenMP_CXX)

# Find MPI
if(GL_ENABLE_MPI AND NOT TARGET MPI::MPI_C)
    find_package(MPI COMPONENTS C REQUIRED)
    target_link_libraries(flags INTERFACE MPI::MPI_C)
endif()


# Setup dependencies
include(cmake/SetupDependenciesCMake.cmake)
#include(cmake/SetupDependenciesConan.cmake)

# Link dependencies
if(NOT TARGET deps)
    add_library(deps INTERFACE)
endif()


find_package(h5pp REQUIRED)
target_link_libraries(deps INTERFACE h5pp::h5pp)


