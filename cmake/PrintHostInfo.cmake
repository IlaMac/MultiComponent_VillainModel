
if(GL_PRINT_INFO)

    # Print host properties
    cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
    cmake_host_system_information(RESULT _proc_type QUERY PROCESSOR_DESCRIPTION)
    cmake_host_system_information(RESULT _os_name QUERY OS_NAME)
    cmake_host_system_information(RESULT _os_release QUERY OS_RELEASE)
    cmake_host_system_information(RESULT _os_version QUERY OS_VERSION)
    cmake_host_system_information(RESULT _os_platform QUERY OS_PLATFORM)
    message(STATUS "| GL HOST INFO")
    message(STATUS "|---------------------------")
    message(STATUS "| CMake Version ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}")
    message(STATUS "| ${_host_name}")
    message(STATUS "| ${_os_name} ${_os_platform} ${_os_release}")
    message(STATUS "| ${_proc_type}")
    message(STATUS "| ${_os_version}")
    message(STATUS "|---------------------------")

endif ()
