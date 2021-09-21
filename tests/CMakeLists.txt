cmake_minimum_required(VERSION 3.12)
project(tests CXX)
add_custom_target(gl-all-tests)


if(GL_ENABLE_H5PP)
    list(APPEND TESTS ${PROJECT_SOURCE_DIR}/h5pp.cpp)
endif()
if(GL_ENABLE_EIGEN3)
    list(APPEND TESTS ${PROJECT_SOURCE_DIR}/eigen.cpp)
endif()

if(GL_ENABLE_OPENMP)
    list(APPEND TESTS ${PROJECT_SOURCE_DIR}/openmp.cpp)
endif()

if(TESTS)
    foreach (test ${TESTS})
        get_filename_component(test_src ${test} NAME)
        get_filename_component(test_nwe ${test} NAME_WE)
        add_executable(gl-${test_nwe} ${test_src} ${test_inc})
        target_link_libraries(gl-${test_nwe} PRIVATE gl-deps gl-flags)
        add_test(NAME gl-${test_nwe} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} COMMAND gl-${test_nwe})
        if (MSVC)
            set_target_properties(gl-${test_nwe} PROPERTIES LINK_FLAGS "/ignore:4099")
        endif()
        add_dependencies(gl-all-tests gl-${test_nwe})
    endforeach ()
    add_custom_command(
            TARGET gl-all-tests
            POST_BUILD
            COMMENT "Running Tests"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            DEPENDS gl-flags gl-deps
            COMMAND ${CMAKE_CTEST_COMMAND} -C $<CONFIG> --output-on-failure)
endif()