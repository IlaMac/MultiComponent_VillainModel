function(detect_os OS)
    # it could be cross compilation
    message(STATUS "CMake-Conan: cmake_system_name=${CMAKE_SYSTEM_NAME}")
    if(CMAKE_SYSTEM_NAME AND NOT CMAKE_SYSTEM_NAME STREQUAL "Generic")
        if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
            set(${OS} Macos PARENT_SCOPE)
        elseif(${CMAKE_SYSTEM_NAME} STREQUAL "QNX")
            set(${OS} Neutrino PARENT_SCOPE)
        else()
            set(${OS} ${CMAKE_SYSTEM_NAME} PARENT_SCOPE)
        endif()
    endif()
endfunction()


function(detect_cxx_standard CXX_STANDARD)
    set(${CXX_STANDARD} ${CMAKE_CXX_STANDARD} PARENT_SCOPE)
    if (CMAKE_CXX_EXTENSIONS)
        set(${CXX_STANDARD} "gnu${CMAKE_CXX_STANDARD}" PARENT_SCOPE)
    endif()
endfunction()


function(detect_compiler COMPILER COMPILER_VERSION)
    if(DEFINED CMAKE_CXX_COMPILER_ID)
        set(_COMPILER ${CMAKE_CXX_COMPILER_ID})
        set(_COMPILER_VERSION ${CMAKE_CXX_COMPILER_VERSION})
    else()
        if(NOT DEFINED CMAKE_C_COMPILER_ID)
            message(FATAL_ERROR "C or C++ compiler not defined")
        endif()
        set(_COMPILER ${CMAKE_C_COMPILER_ID})
        set(_COMPILER_VERSION ${CMAKE_C_COMPILER_VERSION})
    endif()

    message(STATUS "CMake-Conan: CMake compiler=${_COMPILER}")
    message(STATUS "CMake-Conan: CMake compiler version=${_COMPILER_VERSION}")

    if(_COMPILER MATCHES MSVC)
        set(_COMPILER "msvc")
        string(SUBSTRING ${MSVC_VERSION} 0 3 _COMPILER_VERSION)
    elseif(_COMPILER MATCHES AppleClang)
        set(_COMPILER "apple-clang")
        string(REPLACE "." ";" VERSION_LIST ${CMAKE_CXX_COMPILER_VERSION})
        list(GET VERSION_LIST 0 _COMPILER_VERSION)
    elseif(_COMPILER MATCHES Clang)
        set(_COMPILER "clang")
        string(REPLACE "." ";" VERSION_LIST ${CMAKE_CXX_COMPILER_VERSION})
        list(GET VERSION_LIST 0 _COMPILER_VERSION)
    elseif(_COMPILER MATCHES GNU)
        set(_COMPILER "gcc")
        string(REPLACE "." ";" VERSION_LIST ${CMAKE_CXX_COMPILER_VERSION})
        list(GET VERSION_LIST 0 _COMPILER_VERSION)
    endif()

    message(STATUS "CMake-Conan: [settings] compiler=${_COMPILER}")
    message(STATUS "CMake-Conan: [settings] compiler.version=${_COMPILER_VERSION}")

    set(${COMPILER} ${_COMPILER} PARENT_SCOPE)
    set(${COMPILER_VERSION} ${_COMPILER_VERSION} PARENT_SCOPE)
endfunction()

function(detect_build_type BUILD_TYPE)
    if(NOT CMAKE_CONFIGURATION_TYPES)
        # Only set when we know we are in a single-configuration generator
        # Note: we may want to fail early if `CMAKE_BUILD_TYPE` is not defined
        set(${BUILD_TYPE} ${CMAKE_BUILD_TYPE} PARENT_SCOPE)
    endif()
endfunction()


function(detect_host_profile output_file)
    detect_os(MYOS)
    detect_compiler(MYCOMPILER MYCOMPILER_VERSION)
    detect_cxx_standard(MYCXX_STANDARD)
    detect_build_type(MYBUILD_TYPE)

    set(PROFILE "")
    string(APPEND PROFILE "include(default)\n")
    string(APPEND PROFILE "[settings]\n")
    if(MYOS)
        string(APPEND PROFILE os=${MYOS} "\n")
    endif()
    if(MYCOMPILER)
        string(APPEND PROFILE compiler=${MYCOMPILER} "\n")
    endif()
    if(MYCOMPILER_VERSION)
        string(APPEND PROFILE compiler.version=${MYCOMPILER_VERSION} "\n")
    endif()
    if(MYCXX_STANDARD)
        string(APPEND PROFILE compiler.cppstd=${MYCXX_STANDARD} "\n")
    endif()
    if(MYBUILD_TYPE)
        string(APPEND PROFILE "build_type=${MYBUILD_TYPE}\n")
    endif()

    if(NOT DEFINED output_file)
        set(_FN "${CMAKE_BINARY_DIR}/profile")
    else()
        set(_FN ${output_file})
    endif()
    string(TOUPPER "${CMAKE_BUILD_TYPE}" CONAN_BUILD_TYPE)
    set(CONAN_CXXFLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${CONAN_BUILD_TYPE}}")
    set(CONAN_CFLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${CONAN_BUILD_TYPE}}")
    set(CONAN_LDFLAGS "${CMAKE_EXE_LINKER_FLAGS}")
    string(APPEND PROFILE "[buildenv]\n")
    if(CONAN_CFLAGS)
        string(APPEND PROFILE "CFLAGS+=${CONAN_CFLAGS}\n")
    endif()
    if(CONAN_CXXFLAGS)
        string(APPEND PROFILE "CXXFLAGS+=${CONAN_CXXFLAGS}\n")
    endif()
    if(CONAN_LDFLAGS)
        string(APPEND PROFILE "LDFLAGS+=${CONAN_LDFLAGS}\n")
    endif()
    if(CMAKE_C_COMPILER)
        string(APPEND PROFILE "CC=${CMAKE_C_COMPILER}\n")
    endif()
    if(CMAKE_CXX_COMPILER)
        string(APPEND PROFILE "CXX=${CMAKE_CXX_COMPILER}\n")
    endif()
    if(CMAKE_Fortran_COMPILER)
        string(APPEND PROFILE "FC=${CMAKE_Fortran_COMPILER}\n")
    endif()
    string(APPEND PROFILE "[conf]\n")
    string(APPEND PROFILE "tools.cmake.cmaketoolchain:generator=${CMAKE_GENERATOR}\n")

    message(STATUS "CMake-Conan: Creating profile ${_FN}")
    file(WRITE ${_FN} ${PROFILE})
    message(STATUS "CMake-Conan: Profile: \n${PROFILE}")
endfunction()


function(conan_profile_detect_default)
    message(STATUS "CMake-Conan: Checking if a default profile exists")
    execute_process(COMMAND ${CONAN_COMMAND} profile path default
            RESULT_VARIABLE return_code
            OUTPUT_VARIABLE conan_stdout
            ERROR_VARIABLE conan_stderr
            ECHO_ERROR_VARIABLE    # show the text output regardless
            ECHO_OUTPUT_VARIABLE
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    if(NOT ${return_code} EQUAL "0")
        message(STATUS "CMake-Conan: The default profile doesn't exist, detecting it.")
        execute_process(COMMAND ${CONAN_COMMAND} profile detect
                RESULT_VARIABLE return_code
                OUTPUT_VARIABLE conan_stdout
                ERROR_VARIABLE conan_stderr
                ECHO_ERROR_VARIABLE    # show the text output regardless
                ECHO_OUTPUT_VARIABLE
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    endif()
endfunction()


function(conan_install)
    cmake_parse_arguments(ARGS CONAN_ARGS ${ARGN})
    set(CONAN_OUTPUT_FOLDER ${CMAKE_BINARY_DIR}/conan)
    # Invoke "conan install" with the provided arguments
    set(CONAN_ARGS ${CONAN_ARGS} -of=${CONAN_OUTPUT_FOLDER})
    message(STATUS "CMake-Conan: conan install ${CMAKE_SOURCE_DIR} ${CONAN_ARGS} ${ARGN}")
    execute_process(COMMAND ${CONAN_COMMAND} install ${CMAKE_SOURCE_DIR} ${CONAN_ARGS} ${ARGN} --format=json
            RESULT_VARIABLE return_code
            OUTPUT_VARIABLE conan_stdout
            ERROR_VARIABLE conan_stderr
            ECHO_ERROR_VARIABLE    # show the text output regardless
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
    if(NOT "${return_code}" STREQUAL "0")
        message(FATAL_ERROR "Conan install failed='${return_code}'")
    else()
        # the files are generated in a folder that depends on the layout used, if
        # one is specified, but we don't know a priori where this is.
        # TODO: this can be made more robust if Conan can provide this in the json output
        string(JSON CONAN_GENERATORS_FOLDER GET ${conan_stdout} graph nodes 0 generators_folder)
        # message("conan stdout: ${conan_stdout}")
        message(STATUS "CMake-Conan: CONAN_GENERATORS_FOLDER=${CONAN_GENERATORS_FOLDER}")
        set_property(GLOBAL PROPERTY CONAN_GENERATORS_FOLDER "${CONAN_GENERATORS_FOLDER}")
        # reconfigure on conanfile changes
        string(JSON CONANFILE GET ${conan_stdout} graph nodes 0 label)
        message(STATUS "CMake-Conan: CONANFILE=${CMAKE_SOURCE_DIR}/${CONANFILE}")
        set_property(DIRECTORY ${CMAKE_SOURCE_DIR} APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/${CONANFILE}")
        # success
        set_property(GLOBAL PROPERTY CONAN_INSTALL_SUCCESS TRUE)
    endif()
endfunction()


macro(conan_provide_dependency package_name)
    get_property(CONAN_INSTALL_SUCCESS GLOBAL PROPERTY CONAN_INSTALL_SUCCESS)
    if(NOT CONAN_INSTALL_SUCCESS)
        find_program(CONAN_COMMAND "conan" REQUIRED)
        message(STATUS "CMake-Conan: first find_package() found. Installing dependencies with Conan")
        conan_profile_detect_default()
        detect_host_profile(${CMAKE_BINARY_DIR}/conan_host_profile)
        if(NOT CMAKE_CONFIGURATION_TYPES)
            message(STATUS "CMake-Conan: Installing single configuration ${CMAKE_BUILD_TYPE}")
            conan_install(-pr ${CMAKE_BINARY_DIR}/conan_host_profile --build=missing -g CMakeDeps)
        else()
            message(STATUS "CMake-Conan: Installing both Debug and Release")
            conan_install(-pr ${CMAKE_BINARY_DIR}/conan_host_profile -s build_type=Release --build=missing -g CMakeDeps)
            conan_install(-pr ${CMAKE_BINARY_DIR}/conan_host_profile -s build_type=Debug --build=missing -g CMakeDeps)
        endif()
        get_property(CONAN_INSTALL_SUCCESS GLOBAL PROPERTY CONAN_INSTALL_SUCCESS)
        if (CONAN_INSTALL_SUCCESS)
            get_property(CONAN_GENERATORS_FOLDER GLOBAL PROPERTY CONAN_GENERATORS_FOLDER)
            list(PREPEND CMAKE_PREFIX_PATH "${CONAN_GENERATORS_FOLDER}")
        endif()
    else()
        message(STATUS "CMake-Conan: find_package(${ARGV1}) found, 'conan install' already ran")
    endif()

    find_package(${ARGN} BYPASS_PROVIDER)
endmacro()
