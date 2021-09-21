
message(DEBUG "C compiler ${CMAKE_C_COMPILER}")
message(DEBUG "FC compiler ${CMAKE_Fortran_COMPILER}")
message(DEBUG "CXX compiler ${CMAKE_CXX_COMPILER}")

#####################################################
### Set the  microarchitecture for OpenBLAS       ###
#####################################################
cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
if($ENV{CI} OR $ENV{GITHUB_ACTIONS} OR GL_MICROARCH MATCHES "generic|Generic|GENERIC")
    set(MARCH -march=x86-64)
    set(MTUNE -mtune=generic)
elseif(DEFINED GL_MICROARCH)
    set(MARCH -march=${GL_MICROARCH})
    set(MTUNE -mtune=${GL_MICROARCH})
else()
    set(MARCH -march=haswell)
    set(MTUNE -mtune=native)
endif()


message(DEBUG "Using ${MARCH} ${MTUNE}")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${MARCH} ${MTUNE}")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -fno-strict-aliasing -Wall -Wextra -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-strict-aliasing -Wall -Wextra -Wpedantic -fstack-protector -D_FORTIFY_SOURCE=2 -fno-omit-frame-pointer") #-D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
#set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set these variables so that the same flags are used for building dependencies
set(CMAKE_CXX_FLAGS_INIT                 "${CMAKE_CXX_FLAGS_INIT} ${CMAKE_CXX_FLAGS}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELEASE_INIT         "${CMAKE_CXX_FLAGS_RELEASE_INIT} ${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_DEBUG_INIT           "${CMAKE_CXX_FLAGS_DEBUG_INIT} ${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING "" FORCE)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT  "${CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_INIT          "${CMAKE_EXE_LINKER_FLAGS_INIT} ${CMAKE_EXE_LINKER_FLAGS}" CACHE STRING "" FORCE)





# For time tracing using ClangBuildAnalyzer
#if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
#    list(APPEND CMAKE_CXX_FLAGS -ftime-trace)
#endif()

###############################
# Settings for shared builds
###############################

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already (but later on when installing)
# Note: Since DMRG++ is often run from the build folder we want to keep the build-folder RPATH in the executable.
#       Therefore it makes sense to keep this setting "FALSE" here but "TRUE" for dependencies that are
#       installed with in "cmake" mode with externalproject_add
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)


