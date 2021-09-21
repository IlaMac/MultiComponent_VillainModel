

### Enforce an out-of-source build, because c'mon...
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "CMakeTemplate requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
endif()

cmake_minimum_required(VERSION 3.15)
### Setup project name and language
project(Villain_2component CXX)

# Make sure we use GL's own find modules
list(INSERT CMAKE_MODULE_PATH 0  ${PROJECT_SOURCE_DIR}/cmake)


### Set options for compilation
option(BUILD_SHARED_LIBS     "Enables shared library linking (recommended for MPI)"                     ON)
option(GL_ENABLE_TESTS       "Enables CTest testing"                                                    ON)
option(GL_ENABLE_H5PP        "h5pp is a C++ front-end to HDF5"                                          ON)
option(GL_ENABLE_EIGEN3      "Eigen3 is a linear algebra template library"                              ON)
option(GL_ENABLE_SPDLOG      "Spdlog is a fast logging library"                                         ON)
option(GL_ENABLE_OPENMP      "Enable OpenMP flags such as -fopenmp"                                     ON)
option(GL_ENABLE_OPENGL      "Enable OpenGL visualizations"                                             ON)
option(GL_ENABLE_MPI         "Enable MPI parallelization"                                               ON)
option(GL_ENABLE_LTO         "Enable link time optimization"                                            OFF)
option(GL_ENABLE_ASAN        "Enable runtime address sanitizer -fsanitize=address"                      ON)
option(GL_AS_SUBPROJECT      "Try using this as a subproject with add_subdirectory()"                   OFF)
option(GL_PREFIX_ADD_PKGNAME "Install h5pp and dependencies into <CMAKE_INSTALL_PREFIX>/<PackageName>"  ON)
option(GL_PRINT_INFO         "Prints additional information during CMake configure"                     ON)


if(NOT GL_ENABLE_MPI)
    message(FATAL_ERROR "This project uses MPI and should be enabled by default.")
endif()

if(GL_ENABLE_MPI AND NOT BUILD_SHARED_LIBS)
    message(WARNING "This project uses MPI which strongly encourages shared library linking")
endif()
#################################################################
### Make an "enum" for valid package managers                 ###
#################################################################
set(GL_VALID_PACKAGE_MANAGERS find cmake find-or-cmake conan)
set(GL_PACKAGE_MANAGER cmake CACHE STRING "Package manager for external dependencies")
set_property(CACHE GL_PACKAGE_MANAGER PROPERTY STRINGS ${GL_VALID_PACKAGE_MANAGERS})
if(NOT GL_PACKAGE_MANAGER IN_LIST GL_VALID_PACKAGE_MANAGERS)
    message(FATAL_ERROR "GL_PACKAGE_MANAGER must be one of ${GL_VALID_PACKAGE_MANAGERS}. Got [${GL_PACKAGE_MANAGER}]")
endif()


# Used when GL is included as subproject (e.g., as Git submodule/subtree) in the source
# tree of a project that uses it. Users may set the non-cached variable GL_AS_SUBPROJECT
# before add_subdirectory(GL)
if(NOT GL_AS_SUBPROJECT)
    if("^${CMAKE_SOURCE_DIR}$" STREQUAL "^${PROJECT_SOURCE_DIR}$")
        set (AS_SUBPROJECT FALSE)
    else()
        message(STATUS "Detected usage as subproject")
        set (AS_SUBPROJECT TRUE)
    endif()
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


### Add all source files recursively
file(GLOB_RECURSE SOURCES "source/*.cpp" "OpenGL/*.cpp")
file(GLOB_RECURSE HEADERS "source/*.h"   "OpenGL/*.h")
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

if(GL_ENABLE_OPENGL)
    add_library(gl-opengl OBJECT)
    file(GLOB_RECURSE OPENGL_SOURCES "OpenGL/*.cpp")
    file(GLOB_RECURSE OPENGL_HEADERS "OpenGL/*.h")
    target_sources(gl-opengl PRIVATE ${OPENGL_SOURCES} ${OPENGL_HEADERS})
endif()


###########################################
###  Apply RELEASE/DEBUG compile flags  ###
###########################################
include(cmake/CompilerFlags.cmake)

### Create a helper target with all the options you need to compile
### This one is useful to pass to tests and examples so that they compile
### In the same way as our main project.
add_library(gl-flags INTERFACE IMPORTED)



#######################################
# Settings for sanitizers           ###
#######################################
if(GL_ENABLE_ASAN)
    target_compile_options(gl-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>: -fsanitize=address -fno-omit-frame-pointer>)
    target_link_options(gl-flags INTERFACE -fsanitize=address)
endif()



##################################################
###  Add optional RELEASE/DEBUG compile flags  ###
##################################################
# -fsanitize=leak -fsanitize=undefined -fsanitize=alignment  -fsanitize=leak -fsanitize=undefined -fsanitize=alignment
target_compile_options(gl-flags INTERFACE $<$<COMPILE_LANGUAGE:CXX>:-Wall>)
target_compile_options(gl-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELEASE>>:-g>)
target_compile_options(gl-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:DEBUG>>: -Wpedantic -Wextra>)
target_compile_options(gl-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:RELWITHDEBINFO>>:>)
target_compile_options(gl-flags INTERFACE $<$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:MINSIZEREL>>:>)
target_compile_features(gl-flags INTERFACE cxx_std_17)

if(GL_ENABLE_OPENGL)
    target_link_libraries(gl-opengl PUBLIC gl-flags)
endif()


# Fetch or find all the dependencies
include(cmake/SetupDependencies.cmake)

##################################################################
### Link all the things!                                       ###
##################################################################
target_link_libraries(${PROJECT_NAME} PRIVATE gl-deps gl-flags)  # Link everything to main project





# Print summary of CMake configuration
if (GL_PRINT_INFO)
    include(cmake/PrintTargetInfo.cmake)
    include(cmake/getExpandedTarget.cmake)
    expand_target_all_targets(gl-flags GL_FLAGS_EXPANDED)
    expand_target_all_targets(gl-deps GL_DEPS_EXPANDED)
    message(STATUS "| GL TARGET SUMMARY")
    message(STATUS "|--------------------")
    print_compiler_info("| ")
    foreach(tgt ${GL_FLAGS_EXPANDED};${GL_DEPS_EXPANDED})
        print_target_info(${tgt} "| ")
    endforeach()
    print_target_info(${PROJECT_NAME} "| ")
endif()




# Enable link time optimization
if(GL_ENABLE_LTO)
    ### Check if we can enable LTO
    include(CheckIPOSupported)
    check_ipo_supported(RESULT supported OUTPUT error)
    if(supported)
        message(STATUS "IPO / LTO supported, enabling")
        set_target_properties(${PROJECT_NAME} PROPERTIES INTERPROCEDURAL_OPTIMIZATION ON)
    else()
        message(STATUS "IPO / LTO not supported: <${error}>, will not enable")
    endif()
endif()


### CTest testing
if(GL_ENABLE_TESTS AND NOT GL_AS_SUBPROJECT)
    enable_testing()
    add_subdirectory(tests)
endif()


