#ifndef RECTANGLE_2D_H
#define RECTANGLE_2D_H

#include "../common.h"
#include "../util/util.h"

class Rectangle_2D {
  private:

    GLuint VAO, shaderProgram;

  public:

  private:
    void init ();


  public:
    Rectangle_2D ();

    void draw (const glm::mat4 &,
               const glm::vec2 &,
               const float,
               const float,
               const float,
               const glm::vec3 &);
};


#endif