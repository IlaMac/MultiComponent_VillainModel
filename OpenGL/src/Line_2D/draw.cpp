#include "Line_2D.h"

using namespace std;


void Line_2D::draw (
  const glm::mat4 & PV,
  const glm::vec2 & p1,
  const glm::vec2 & p2,
  const glm::vec3 & color,
  const float       width
) {
  // bind VAO
  glBindVertexArray(this->VAO);

  // select shader program
  glUseProgram(this->shaderProgram);

  // orient in proper direction
  const float alpha = atan2(p2[1] - p1[1], p2[0] - p1[0]);
  const glm::mat4 rotationMat = glm::rotate(glm::mat4{1}, alpha, glm::vec3(0.f, 0.f, 1.f));

  // position line
  const auto model = glm::translate(glm::mat4{1}, glm::vec3{p1, 0})
                   * rotationMat
                   * glm::scale(glm::mat4{1}, glm::vec3{glm::distance(p1, p2), 1, 1});

  // compute mvp
  glm::mat4 mvp = PV * model;

  // send mvp
  int mvpLoc = glGetUniformLocation(this->shaderProgram, "mvp");
  glUniformMatrix4fv(mvpLoc, 1, GL_FALSE, glm::value_ptr(mvp));

  // send color
  int colorLoc = glGetUniformLocation(this->shaderProgram, "color");
  glUniform3fv(colorLoc, 1, glm::value_ptr(color));

  // disable blending
  glDisable(GL_BLEND);

  // line width
  glLineWidth(width);

  // draw
  glDrawElements(GL_LINES, 2, GL_UNSIGNED_INT, 0);


  // unbind VAO
  glBindVertexArray(0);
}