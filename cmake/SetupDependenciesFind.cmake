

##################################################################
### Preempt Threads::Threads                                   ###
### It's looked for in dependencies, so we make it right       ###
### before it's done wrong, i.e. with pthread instead of       ###
### -lpthread which causes link errors downstream with         ###
###    -Wl,--whole-archive.... -Wl,--no-whole-archive          ###
##################################################################
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
target_link_libraries(Threads::Threads INTERFACE rt dl)

if(GL_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS CXX REQUIRED)
endif()

if(GL_ENABLE_MPI)
    find_package(MPI REQUIRED)
endif()

if(GL_ENABLE_OPENGL)
    if(NOT TARGET gl-opengl)
        add_library(gl-opengl INTERFACE)
    endif()

    find_package(OpenGL COMPONENTS EGL GLX OpenGL REQUIRED)
    find_package(X11 REQUIRED)
endif()


if(GL_PACKAGE_MANAGER MATCHES "find")

    if(GL_PACKAGE_MANAGER STREQUAL "find")
        set(REQUIRED REQUIRED)
    endif()

    if(GL_ENABLE_H5PP)
        find_package(h5pp 1.9.0 ${REQUIRED})
    endif()

    if(GL_ENABLE_OPENGL)

        set(GLM_ROOT        ${PROJECT_SOURCE_DIR}/include)
        set(GLFW3_ROOT      ${PROJECT_SOURCE_DIR}/include)
        set(Freetype_ROOT   ${PROJECT_SOURCE_DIR}/include)

        find_package(GLM ${REQUIRED})
        find_package(GLFW3 ${REQUIRED})
        find_package(Freetype ${REQUIRED})

        add_library(GLAD OBJECT)
        target_sources(GLAD PUBLIC ${PROJECT_SOURCE_DIR}/thirdparty/glad.cpp)
        set(GLAD_FOUND)

        if(GLFW3_FOUND)
            target_link_libraries(gl-opengl PUBLIC ${GLFW3_LIBRARY} )
            target_include_directories(gl-opengl PUBLIC ${GLFW3_INCLUDE_DIR} )
        endif()

        if(GLM_FOUND)
            target_include_directories(gl-opengl PUBLIC ${GLM_INCLUDE_DIRS})
            target_compile_definitions(gl-opengl PUBLIC GLM_ENABLE_EXPERIMENTAL)
        endif()

        if(GLAD_FOUND)
            target_compile_definitions(gl-opengl PUBLIC $<TARGET_OBJECTS:GLAD>)
        endif()
        if(Freetype_FOUND)
            target_link_libraries(gl-opengl PUBLIC Freetype::Freetype)
        endif()
        if(OpenGL_FOUND)
            target_link_libraries(gl-opengl PUBLIC OpenGL::EGL OpenGL::GL OpenGL::GLX OpenGL::OpenGL)
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

        target_compile_definitions(gl-opengl PUBLIC GUI_ENABLED)
    endif()

endif()


