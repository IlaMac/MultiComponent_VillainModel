#ifndef SPHERE_H
#define SPHERE_H

#include "../common.h"
#include "../util/util.h"

class Sphere {
  private:
    const unsigned numThetas,
                   numPhis;


    GLuint VAO, shaderProgram_fill, shaderProgram_frame;

  public:

  private:
    void init ();


  public:
    Sphere (const unsigned = 20,
            const unsigned = 40);

    void draw (const glm::vec3 &,
               const glm::mat4 &,
               const glm::vec3 &,
               const float,
               const glm::vec3 &,
               const bool = false);

};


#endif