#include "Window.h"

using namespace std;


////
//// error callback function
////
static void error_callback (
  int error,
  const char* description
) {
  fprintf(stderr, "Error: %s\n", description);
}

Window::Window (
  const int      width,
  const int      height,
  const string & title
) {

  ////
  //// set error callback function
  ////
  glfwSetErrorCallback(error_callback);


  ////
  //// ensure successful GLFW initialization
  ////
  if ( ! glfwInit()) {
    cout << "GLFW initialization failed" << endl;
    exit(EXIT_FAILURE);
  }


  ////
  //// require a minimum OpenGL version
  ////
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


  ////
  //// forward compatibility for OSX
  ////
  #ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
  #endif


  ////
  ////
  ////
  this->glfwWindow = glfwCreateWindow(width, height, title.data(), NULL, NULL);


  ////
  //// update width and height
  ////
  glfwGetFramebufferSize(this->glfwWindow, &this->width, &this->height);


  // callback when window is eventually closed
  glfwSetWindowCloseCallback(this->glfwWindow, Window::GLFWwindowclosefun);

  if ( ! this->glfwWindow) {
    cout << "Failed to create GLFW window" << endl;
    glfwTerminate();
    exit(EXIT_FAILURE);
  }


  ////
  //// initialize context etc.
  ////
  this->init();
}



Window::~Window () {
  ////
  //// terminate window
  ////
  glfwDestroyWindow(this->glfwWindow);

  ////
  //// terminate GLFW context
  ////
  glfwTerminate();


  ////
  //// delete VOAs and shader programs
  ////
  delete this->sphere;
  delete this->cylinder;
  delete this->rectangle;
  delete this->line;
  delete this->text;
  delete this->line_2D;
  delete this->rectangle_2D;
  delete this->arrowHead_2D;
}

unsigned long Window::getFrameNumber () const { return this->frameNum; }

float Window::getFramerate () const { return this->framerate; }

float Window::getWidth () const { return this->width; }
float Window::getHeight () const { return this->height; }