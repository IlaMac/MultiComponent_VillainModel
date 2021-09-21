
if(GL_PACKAGE_MANAGER MATCHES "conan")

    unset(CONAN_BUILD_INFO)
    unset(CONAN_BUILD_INFO CACHE)
    find_file(CONAN_BUILD_INFO
            conanbuildinfo.cmake
            HINTS ${CMAKE_BINARY_DIR} ${CMAKE_CURRENT_LIST_DIR}
            NO_DEFAULT_PATH)

    if(CONAN_BUILD_INFO)
        ##################################################################
        ### Use pre-existing conanbuildinfo.cmake                      ###
        ### This avoids having to run conan again                      ###
        ##################################################################
        message(STATUS "Detected Conan build info: ${CONAN_BUILD_INFO}")
        include(${CONAN_BUILD_INFO})
        conan_basic_setup(TARGETS KEEP_RPATHS NO_OUTPUT_DIRS)
    else()

        ##################################################################
        ### Install conan-modules/conanfile.txt dependencies          ###
        ##################################################################

        find_program (
                CONAN_COMMAND
                conan
                HINTS ${CONAN_PREFIX} $ENV{CONAN_PREFIX} ${CONDA_PREFIX} $ENV{CONDA_PREFIX}
                PATHS $ENV{HOME}/anaconda3 $ENV{HOME}/miniconda3 $ENV{HOME}/.conda
                PATH_SUFFIXES bin envs/dmrg/bin
        )
        message(STATUS "Found conan: ${CONAN_COMMAND}")

        # Download cmake-conan automatically, you can also just copy the conan.cmake file
        if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
            message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
            file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.16.1/conan.cmake"
                    "${CMAKE_BINARY_DIR}/conan.cmake")
        endif()

        if(BUILD_SHARED_LIBS)
            list(APPEND GL_CONAN_OPTIONS OPTIONS "*:shared=True")
        else()
            list(APPEND GL_CONAN_OPTIONS OPTIONS "*:shared=False")
        endif()

        set(GL_CONANFILE_TXT conanfile.txt)
        if(GL_ENABLE_OPENGL)
            set(GL_CONANFILE_TXT conanfile_opengl.txt)
        endif()


        include(${CMAKE_BINARY_DIR}/conan.cmake)

        conan_cmake_run(
                CONANFILE ${GL_CONANFILE_TXT}
                CONAN_COMMAND ${CONAN_COMMAND}
                BUILD_TYPE ${CMAKE_BUILD_TYPE}
                BASIC_SETUP CMAKE_TARGETS
                SETTINGS compiler.cppstd=17
                SETTINGS compiler.libcxx=libstdc++11
                PROFILE_AUTO ALL
                ${GL_CONAN_OPTIONS}
                KEEP_RPATHS
                NO_OUTPUT_DIRS
                BUILD missing
        )
    endif()


    # Make aliases or link
    add_library(h5pp::h5pp          ALIAS CONAN_PKG::h5pp)
    if(GL_ENABLE_OPENGL)
        target_link_libraries(gl-opengl PUBLIC CONAN_PKG::glad)
        target_link_libraries(gl-opengl PUBLIC CONAN_PKG::glm)
        target_link_libraries(gl-opengl PUBLIC CONAN_PKG::glfw)
        target_link_libraries(gl-opengl PUBLIC CONAN_PKG::opengl)
        target_link_libraries(gl-opengl PUBLIC CONAN_PKG::xorg)
    endif()

endif()
