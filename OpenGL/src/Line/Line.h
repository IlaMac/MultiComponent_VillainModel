#ifndef LINE_H
#define LINE_H

#include "../common.h"
#include "../util/util.h"

class Line {
  private:

    GLuint VAO, shaderProgram;

  public:

  private:
    void init ();


  public:
    Line ();

    void draw (const glm::mat4 &,
               const glm::vec3 &,
               const glm::vec3 &,
               const glm::vec3 &,
               const float = 1);
};


#endif