#include "ArrowHead_2D.h"

using namespace std;


void ArrowHead_2D::draw (
  const glm::mat4 & PV,
  const glm::vec2 & p0,
  const float       length,
  const float       rotAngle,
  const float       shearAngle,
  const glm::vec3 & color
) {
  // bind VAO
  glBindVertexArray(this->VAO);

  // select shader program
  glUseProgram(this->shaderProgram);

  // orient in proper direction
  const glm::mat4 rotationMat = glm::rotate(glm::mat4{1}, (float) (rotAngle - 0.75*M_PI), glm::vec3(0.f, 0.f, 1.f));

  // shear and scale
  const float a = 0.25*M_PI - 0.5*shearAngle;
  const float w = length / cos(0.25*M_PI - a);
  float shearAndScale_data[16] = {
     w*cos(a), w*sin(a), 0, 0,
     w*sin(a), w*cos(a), 0, 0,
     0,        0,        1, 0,
     0,        0,        0, 1
  };
  glm::mat4 shearAndScale;
  memcpy(glm::value_ptr(shearAndScale), shearAndScale_data, sizeof(shearAndScale_data));

  // position line
  const auto model = glm::translate(glm::mat4{1}, glm::vec3{p0, 0})
                   * rotationMat
                   * shearAndScale;

  // compute mvp
  glm::mat4 mvp = PV * model;

  // send mvp
  int mvpLoc = glGetUniformLocation(this->shaderProgram, "mvp");
  glUniformMatrix4fv(mvpLoc, 1, GL_FALSE, glm::value_ptr(mvp));

  // send color
  int colorLoc = glGetUniformLocation(this->shaderProgram, "color");
  glUniform3fv(colorLoc, 1, glm::value_ptr(color));

  // fill mode
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // disable blending
  glDisable(GL_BLEND);

  // draw
  glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);


  // unbind VAO
  glBindVertexArray(0);
}