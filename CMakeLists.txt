cmake_minimum_required(VERSION 3.12)
project(Villain_2component)
enable_language(C CXX)

### Set options for compilation
option(BUILD_SHARED_LIBS     "Enables shared library linking (recommended for MPI)"     ON)
option(GL_ENABLE_TESTS       "Enables CTest testing"                                    ON)
option(GL_ENABLE_H5PP        "h5pp is a C++ front-end to HDF5"                          ON)
option(GL_ENABLE_EIGEN3      "Eigen3 is a linear algebra template library"              ON)
option(GL_ENABLE_SPDLOG      "Spdlog is a fast logging library"                         ON)
option(GL_ENABLE_OPENMP      "Enable OpenMP flags such as -fopenmp"                     ON)
option(GL_ENABLE_MPI         "Enable MPI parallelization"                               ON)
option(GL_ENABLE_LTO         "Enable link time optimization"                            OFF)
option(GL_ENABLE_ASAN        "Enable runtime address sanitizer -fsanitize=address"      ON)
option(GL_AS_SUBPROJECT      "Try using this as a subproject with add_subdirectory()"   OFF)
option(GL_DEPS_IN_SUBDIR     "Install dependencies into CMAKE_INSTALL_PREFIX/<libname>" OFF)
option(GL_PRINT_INFO         "Prints additional information during CMake configure"     OFF)

if(GL_ENABLE_MPI AND NOT BUILD_SHARED_LIBS)
    message(WARNING "This project uses MPI which strongly encourages shared library linking")
endif()
#################################################################
### Make an "enum" for valid download methods                 ###
#################################################################
set(GL_VALID_PACKAGE_MANAGERS find cmake conan)
set(GL_PACKAGE_MANAGER cmake CACHE STRING "Package manager for external dependencies")
set_property(CACHE GL_PACKAGE_MANAGER  PROPERTY STRINGS ${GL_VALID_PACKAGE_MANAGERS})
if(NOT GL_PACKAGE_MANAGER IN_LIST GL_VALID_PACKAGE_MANAGERS)
    message(FATAL_ERROR "GL_PACKAGE_MANAGER must be one of ${GL_VALID_PACKAGE_MANAGERS}")
endif()


################################################################
### Get git version number                                   ###
### Generates a header gitversion/gitversion.h               ###
### Include it using #include <gitversion.h>                 ###
### Gives a namespace GIT:: with several git identifiers     ###
################################################################
include(cmake/gitversion.cmake)


############################################################
### Print extra info during CMake configure              ###
############################################################
include(cmake/PrintHostInfo.cmake)
include(cmake/PrintOptionInfo.cmake)


### Print compiler details
message("Compiling with: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")



###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
include(cmake/CompilerFlags.cmake)

### Add all source files recursively
file(GLOB_RECURSE SOURCES "source/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS} source/main.cpp source/main.h source/villain_MC.cpp source/villain_MC.h source/pt.cpp source/pt.h source/measures.cpp source/measures.h source/initialization.cpp source/initialization.h source/rng.cpp)



##################################################################
###  Force cmake to find .a library suffixes on static builds  ###
##################################################################
if(BUILD_SHARED_LIBS)
    set(LINK_TYPE SHARED)
else()
    set(LINK_TYPE STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
    #target_link_libraries  (project-settings INTERFACE -static)   ### Static linkage
endif()


# Setup search paths
include(cmake/SetupPaths.cmake)

# Find or install all the dependencies
include(cmake/SetupDependencies.cmake)

# Link all the things!
target_link_libraries(${PROJECT_NAME} PRIVATE deps flags)

# Enable link-time optimization
target_enable_lto(${PROJECT_NAME})

# Print summary of CMake configuration
include(cmake/PrintTargetInfo.cmake)
print_project_summary(${PROJECT_NAME})

### CTest testing
if(GL_ENABLE_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()


