#ifndef ARROW_HEAD_2D_H
#define ARROW_HEAD_2D_H

#include "../common.h"
#include "../util/util.h"

class ArrowHead_2D {
  private:

    GLuint VAO, shaderProgram;

  public:

  private:
    void init ();


  public:
    ArrowHead_2D ();

    void draw (const glm::mat4 &,
               const glm::vec2 &,
               const float,
               const float,
               const float,
               const glm::vec3 &);
};


#endif