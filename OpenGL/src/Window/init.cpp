#include "Window.h"

using namespace std;


void Window::init () {

  ////
  //// set callbacks
  ////
  glfwSetWindowUserPointer(this->glfwWindow, this);
  glfwSetKeyCallback(this->glfwWindow, onKey);
  glfwSetCursorPosCallback(this->glfwWindow, onMouseMove);
  glfwSetMouseButtonCallback(this->glfwWindow, onMouseButton);
  glfwSetScrollCallback(this->glfwWindow, onScroll);


  ////
  //// need to be assigned before using GLAD
  ////
  glfwMakeContextCurrent(this->glfwWindow);


  ////
  //// load modern OpenGL extensions using GLAD
  ////
  if ( ! gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    cout << "Failed to initialize GLAD" << endl;
    glfwTerminate();
    exit(EXIT_FAILURE);
  }


  ////
  //// output OpenGL version
  ////
  if (false) {
    printf("OpenGL v%d.%d\n", GLVersion.major, GLVersion.minor);
    printf("Vendor:   %s\n",    glGetString(GL_VENDOR));
    printf("Renderer: %s\n",  glGetString(GL_RENDERER));
    printf("Version:  %s\n",   glGetString(GL_VERSION));
    printf("GLSL:     %s\n",      glGetString(GL_SHADING_LANGUAGE_VERSION));
  }


  // the number of screen updates to wait from the time glfwSwapBuffers was called before swapping the buffers and returning
  // (GLFW by default use double buffering)
  glfwSwapInterval(1);

  // set the background color (RGBA)
  glClearColor(0.f, 0.0f, 0.f, 1.0f);



  // // Set OpenGL options
  // glEnable(GL_CULL_FACE);
  // glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);



  ////
  //// initialize VOAs and shader programs
  ////
  this->sphere       = new Sphere{20, 20};
  this->cylinder     = new Cylinder{20, 1};
  this->rectangle    = new Rectangle{};
  this->line         = new Line{};
  this->text         = new Text{};
  this->line_2D      = new Line_2D{};
  this->rectangle_2D = new Rectangle_2D{};
  this->arrowHead_2D = new ArrowHead_2D{};
}