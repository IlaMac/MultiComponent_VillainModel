#ifndef CONE_H
#define CONE_H

#include "../common.h"
#include "../util/util.h"

class Cone {
  private:
    const unsigned numPhis;


    GLuint VAO, shaderProgram_fill, shaderProgram_frame;

  public:

  private:
    void init ();


  public:
    Cone (const unsigned = 20);

    void draw (const glm::vec3 &,
               const glm::mat4 &,
               const glm::vec3 &,
               const glm::vec3 &,
               const float,
               const glm::vec3 &,
               const bool = false);

};


#endif