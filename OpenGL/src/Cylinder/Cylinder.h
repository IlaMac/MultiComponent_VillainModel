#ifndef CYLINDER_H
#define CYLINDER_H

#include "../common.h"
#include "../util/util.h"

class Cylinder {
  private:
    const unsigned numPhis,
                   numStacks;


    GLuint VAO, shaderProgram_fill, shaderProgram_frame;

  public:

  private:
    void init ();


  public:
    Cylinder (const unsigned = 5,
              const unsigned = 1);

    void draw (const glm::vec3 &,
               const glm::mat4 &,
               const glm::vec3 &,
               const glm::vec3 &,
               const float,
               const glm::vec3 &,
               const bool = false);

};


#endif