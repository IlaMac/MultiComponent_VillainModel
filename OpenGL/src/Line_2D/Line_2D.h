#ifndef LINE_2D_H
#define LINE_2D_H

#include "../common.h"
#include "../util/util.h"

class Line_2D {
  private:

    GLuint VAO, shaderProgram;

  public:

  private:
    void init ();


  public:
    Line_2D ();

    void draw (const glm::mat4 &,
               const glm::vec2 &,
               const glm::vec2 &,
               const glm::vec3 &,
               const float = 1);
};


#endif