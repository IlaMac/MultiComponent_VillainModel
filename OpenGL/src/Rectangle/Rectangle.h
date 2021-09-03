#ifndef RECTANGLE_H
#define RECTANGLE_H

#include "../common.h"
#include "../util/util.h"

class Rectangle {
  private:

    GLuint VAO, shaderProgram_fill, shaderProgram_frame;

  public:

  private:
    void init ();


  public:
    Rectangle ();

    void draw (const glm::vec3 &,
               const glm::mat4 &,
               const glm::vec3 &,
               const glm::vec3 &,
               const glm::vec3 &,
               const bool = false);
};

#endif