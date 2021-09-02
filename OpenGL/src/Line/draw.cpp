#include "Line.h"

using namespace std;


void Line::draw (
  const glm::mat4 & PV,
  const glm::vec3 & p1,
  const glm::vec3 & p2,
  const glm::vec3 & color,
  const float       width
) {
  // bind VAO
  glBindVertexArray(this->VAO);

  // select shader program
  glUseProgram(this->shaderProgram);

  // orient in proper direction
  const auto dot = glm::dot({0, 0, 1}, glm::normalize(p2 - p1));
  glm::mat4 rotationMat(1);
  if (1 - 1E-5 > dot) {
    const auto cross = glm::cross({0, 0, 1}, glm::normalize(p2 - p1));
    rotationMat = glm::rotate(rotationMat, acos(dot), cross);
  }

  // position line
  const auto model = glm::translate(glm::mat4{1}, p1)
                   * rotationMat
                   * glm::scale(glm::mat4{1}, glm::vec3{1, 1, glm::distance(p1, p2)});

  // compute mvp
  const auto mvp = PV * model;

  // send mvp
  int mvpLoc = glGetUniformLocation(this->shaderProgram, "mvp");
  glUniformMatrix4fv(mvpLoc, 1, GL_FALSE, glm::value_ptr(mvp));

  // send color
  int colorLoc = glGetUniformLocation(this->shaderProgram, "color");
  glUniform3fv(colorLoc, 1, glm::value_ptr(color));

  glDepthMask(GL_TRUE);

  // disable blending
  glDisable(GL_BLEND);

  // draw
  glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, 0);


  // unbind VAO
  glBindVertexArray(0);
}