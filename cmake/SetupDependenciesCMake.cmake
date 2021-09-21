
if (GL_PACKAGE_MANAGER MATCHES "cmake")
    include(cmake/InstallPackage.cmake)
    if(GL_PREFIX_ADD_PKGNAME)
        set(INSTALL_PREFIX_PKGNAME INSTALL_PREFIX_PKGNAME)
    endif()
    if(GL_ENABLE_H5PP)
        # h5pp for writing to file binary in format
        install_package(h5pp VERSION 1.9.0
                CMAKE_ARGS
                -DEigen3_ROOT:PATH=${GL_DEPS_INSTALL_DIR}
                -DH5PP_PACKAGE_MANAGER:STRING=cmake
                -DH5PP_ENABLE_EIGEN3:BOOL=${GL_ENABLE_EIGEN3}
                -DH5PP_ENABLE_SPDLOG:BOOL=${GL_ENABLE_SPDLOG}
                -DH5PP_ENABLE_MPI:BOOL=${GL_ENABLE_MPI}
                -DCMAKE_VERBOSE_MAKEFILE=${CMAKE_VERBOSE_MAKEFILE}
                ${INSTALL_PREFIX_PKGNAME})
    endif()

    if(GL_ENABLE_OPENGL)
        find_package(OpenGL COMPONENTS EGL GLX OpenGL REQUIRED)
        find_package(X11 REQUIRED)

        install_package(glm      VERSION         ${INSTALL_PREFIX_PKGNAME})
        install_package(Freetype                 ${INSTALL_PREFIX_PKGNAME})
        install_package(glad     VERSION 0.1.34  ${INSTALL_PREFIX_PKGNAME})
        install_package(glfw3    VERSION 3.3.4   TARGET_NAME glfw ${INSTALL_PREFIX_PKGNAME})

        target_link_libraries(gl-opengl PUBLIC glm::glm)
        target_link_libraries(gl-opengl PUBLIC Freetype::Freetype)
        target_link_libraries(gl-opengl PUBLIC glad::glad)
        target_link_libraries(gl-opengl PUBLIC glfw)

        if(OpenGL_FOUND)
            target_link_libraries(gl-opengl INTERFACE OpenGL::EGL OpenGL::GL OpenGL::GLX OpenGL::OpenGL)
        endif()

        if(X11_FOUND)
            target_link_libraries(gl-opengl PUBLIC
                    X11::Xext
                    X11::X11
                    X11::Xi
                    X11::Xrandr
                    X11::Xinerama
                    X11::Xcursor
                    rt m)
            if(X11_Xxf86vm_FOUND) # This one seems a bit problematic
                target_link_libraries(gl-opengl PUBLIC X11::Xxf86vm)
            endif()
        endif()

        target_compile_definitions(gl-opengl PUBLIC GLM_ENABLE_EXPERIMENTAL)
        target_compile_definitions(gl-opengl PUBLIC GUI_ENABLED)

    endif()

endif ()
