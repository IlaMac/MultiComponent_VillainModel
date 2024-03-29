function(build_dependency dep_name install_dir extra_flags)
    set(build_dir    ${CMAKE_BINARY_DIR}/deps-build/${dep_name})

    if (GL_DEPS_IN_SUBDIR) # h5pp is run with append libsuffix so we don't need to append it again
        set(install_dir ${install_dir}/${dep_name})
        mark_as_advanced(install_dir)
    endif()
    include(cmake/GetNumThreads.cmake)
    get_num_threads(num_threads)

    if(NOT DEFINED ENV{CC} AND DEFINED CMAKE_C_COMPILER)
        set(ENV{CC} ${CMAKE_C_COMPILER})
    endif()
    if(NOT DEFINED ENV{CXX} AND DEFINED CMAKE_CXX_COMPILER)
        set(ENV{CXX} ${CMAKE_CXX_COMPILER})
    endif()
    if(NOT DEFINED ENV{FC} AND DEFINED CMAKE_Fortran_COMPILER)
        set(ENV{FC} ${CMAKE_Fortran_COMPILER})
    endif()

    execute_process( COMMAND  ${CMAKE_COMMAND} -E remove ${build_dir}/CMakeCache.txt)
    execute_process( COMMAND  ${CMAKE_COMMAND} -E make_directory ${build_dir})
    execute_process(
            COMMAND  ${CMAKE_COMMAND}
            # CMake flags
            -DCMAKE_POLICY_DEFAULT_CMP0074=NEW
            -DCMAKE_EXE_LINKER_FLAGS_INIT=${CMAKE_EXE_LINKER_FLAGS}
            -DCMAKE_SHARED_LINKER_FLAGS_INIT=${CMAKE_SHARED_LINKER_FLAGS}
            -DCMAKE_STATIC_LINKER_FLAGS_INIT=${CMAKE_STATIC_LINKER_FLAGS}
            -DCMAKE_MODULE_LINKER_FLAGS_INIT=${CMAKE_MODULE_LINKER_FLAGS}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
            -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=${CMAKE_POSITION_INDEPENDENT_CODE}
            -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
            -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
            -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
            -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
            -DCMAKE_CXX_FLAGS_INIT:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_CXX_FLAGS_RELEASE_INIT:STRING=${CMAKE_CXX_FLAGS_RELEASE}
            -DCMAKE_CXX_FLAGS_DEBUG_INIT:STRING=${CMAKE_CXX_FLAGS_DEBUG}
            -DCMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT:STRING=${CMAKE_CXX_FLAGS_RELWITHDEBINFO}
            -DCMAKE_CXX_FLAGS_MINSIZEREL_INIT:STRING=${CMAKE_CXX_FLAGS_MINSIZEREL}
            -DCMAKE_INSTALL_MESSAGE=NEVER #Avoid unnecessary output to console
            -DCMAKE_GENERATOR=${CMAKE_GENERATOR}
            -DCMAKE_GENERATOR_PLATFORM=${CMAKE_GENERATOR_PLATFORM}
            -DCMAKE_INSTALL_PREFIX:PATH=${install_dir}
            ${extra_flags}
            ${PROJECT_SOURCE_DIR}/cmake/external_${dep_name}
            WORKING_DIRECTORY ${build_dir}
            RESULT_VARIABLE config_result
    )
    if(config_result)
        message(STATUS "Got non-zero exit code while configuring ${dep_name}")
        message(STATUS  "build_dir         : ${build_dir}")
        message(STATUS  "install_dir       : ${install_dir}")
        message(STATUS  "extra_flags       : ${extra_flags}")
        message(STATUS  "config_result     : ${config_result}")
        message(FATAL_ERROR "Failed to configure ${dep_name}")
    endif()




    set(ENV{CMAKE_BUILD_PARALLEL_LEVEL} ${num_threads})
    execute_process(COMMAND  ${CMAKE_COMMAND} --build . --parallel ${num_threads}
            WORKING_DIRECTORY "${build_dir}"
            RESULT_VARIABLE build_result
            )

    if(build_result)
        message(STATUS "Got non-zero exit code while building ${dep_name}")
        message(STATUS  "build_dir         : ${build_dir}")
        message(STATUS  "install_dir       : ${install_dir}")
        message(STATUS  "extra_flags       : ${extra_flags}")
        message(STATUS  "build_result      : ${build_result}")
        message(FATAL_ERROR "Failed to build ${dep_name}")
    endif()

    # Copy the install manifest if it exists
    file(GLOB_RECURSE INSTALL_MANIFEST "${build_dir}/*/install_manifest*.txt")
    foreach(manifest ${INSTALL_MANIFEST})
        get_filename_component(manifest_filename ${manifest} NAME_WE)
        message(STATUS "Copying install manifest: ${manifest}")
        configure_file(${manifest} ${CMAKE_CURRENT_BINARY_DIR}/${manifest_filename}_${dep_name}.txt)
    endforeach()

endfunction()