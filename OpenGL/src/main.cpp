// #include "common.h"
// #include "Window/Window.h"

// using namespace std;


// ////
// //// error callback function
// ////
// static void error_callback (
//   int error,
//   const char* description
// ) {
//   fprintf(stderr, "Error: %s\n", description);
// }


// int main (
//   int argc,
//   char *argv[]
// ) {
//   ////
//   //// set error callback function
//   ////
//   glfwSetErrorCallback(error_callback);


//   ////
//   //// ensure successful GLFW initialization
//   ////
//   if ( ! glfwInit()) {
//     cout << "GLFW initialization failed" << endl;
//     exit(EXIT_FAILURE);
//   }


//   ////
//   //// require a minimum OpenGL version
//   ////
//   glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
//   glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
//   glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


//   ////
//   //// forward compatibility for OSX
//   ////
//   #ifdef __APPLE__
//     glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
//   #endif


//   Window window{1000, 800, "Mr. Window"};
//   window.show();

//   ////
//   //// terminate GLFW
//   ////
//   glfwTerminate();
//   exit(EXIT_SUCCESS);
// }