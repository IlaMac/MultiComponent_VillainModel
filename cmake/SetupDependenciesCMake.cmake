if (GL_PACKAGE_MANAGER MATCHES "cmake")
    include(cmake/InstallPackage.cmake)

    list(APPEND h5pp_ARGS -DEigen3_ROOT:PATH=${GL_DEPS_INSTALL_DIR})
    list(APPEND h5pp_ARGS -DH5PP_PACKAGE_MANAGER:STRING=${GL_PACKAGE_MANAGER})
    list(APPEND h5pp_ARGS -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE})

    # h5pp for writing to file binary in format
    install_package(h5pp VERSION 1.11.0 CMAKE_ARGS ${h5pp_ARGS})
endif ()

