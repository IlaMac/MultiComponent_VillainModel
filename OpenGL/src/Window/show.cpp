#include "Window.h"

using namespace std;


void Window::show (
  bool (* updateAndRender) (void *, const double &),
  void * self
) {

  ////
  //// start timer
  ////
  const auto t_init = chrono::high_resolution_clock::now();
  auto t_prev_frame = t_init;
  while ( ! glfwWindowShouldClose(this->glfwWindow)) {
    /* TEMP */ glfwMakeContextCurrent(this->glfwWindow);

    ////
    //// get number seconds since init
    ////
    const auto now = chrono::high_resolution_clock::now();
    const auto t = chrono::duration_cast<chrono::duration<double> >(now - t_init).count();

    ////
    //// get number of seconds since last frame
    ////
    const auto inv_frame_rate = chrono::duration_cast<chrono::duration<double> >(now - t_prev_frame).count();
    this->framerate = 1 / inv_frame_rate;

    ////
    //// get size of window and retrieve frame buffer size
    ////
    glfwGetFramebufferSize(this->glfwWindow, &this->width, &this->height);
    glViewport(0, 0, this->width, this->height);


    ////
    //// overwrite all pixels in the buffer using this color
    ////
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    ////
    //// enable depth test
    ////
    glEnable(GL_DEPTH_TEST);


    ////
    //// process event queue with appropriate callbacks
    ////
    glfwPollEvents();


    ////
    //// update world
    ////

    // projection
    const glm::mat4 proj = glm::perspective(camera.getVerticalFOV(),
                                            (float) this->width / (float) height,
                                            0.1f,
                                            1000.0f);

    // view
    const glm::mat4 view = glm::lookAt(camera.pos,
                                       camera.pos + camera.front,   // target = position + direction
                                       camera.up);                  // up vector

    // combine
    this->VP = proj * view;


    ////
    //// update and render the world from wherever
    ////
    const auto ok = (*updateAndRender)(self, t);

    // if we should continue drawing the frame
    if ( ! ok) continue;


    ////
    //// update camera position and direction
    ////
    camera.update();


    ////
    //// swap buffers
    ////
    glfwSwapBuffers(this->glfwWindow);


    ////
    //// update
    ////
    t_prev_frame = now;
    this->frameNum++;
  }

}
