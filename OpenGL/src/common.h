////
//// remove this imho
////


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <random>
#include <chrono>


#include <glad/glad.h>    // extension loader library (replacing the OpenGL headers of GLFW which might not detect the newest ones)
#include <GLFW/glfw3.h>   // includes OpenGL header
#pragma GCC diagnostic push                   // disables warning temporarily
#pragma GCC diagnostic ignored "-Wpedantic"   //
#define GLM_FORCE_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/perpendicular.hpp>
#pragma GCC diagnostic pop                    // enables warning once again


#include "etc/cout.h"