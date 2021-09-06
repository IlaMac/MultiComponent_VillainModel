#include "Cone.h"

using namespace std;


void Cone::draw (
  const glm::vec3 & camPos,
  const glm::mat4 & PV,
  const glm::vec3 & p_base,
  const glm::vec3 & p_top,
  const float       radius,
  const glm::vec3 & color,
  const bool        useWireFrames
) {
  // bind VAO
  glBindVertexArray(this->VAO);

  // orient in proper direction
  const auto dot = glm::dot({0, 0, 1}, glm::normalize(p_top - p_base));
  glm::mat4 rotationMat(1);
  if (1 - 1E-5 > dot) {
    const auto cross = glm::cross({0, 0, 1}, glm::normalize(p_top - p_base));
    rotationMat = glm::rotate(rotationMat, acos(dot), cross);
  }

  // position cylinder
  const auto model = glm::translate(glm::mat4{1}, p_base)
                   * rotationMat
                   * glm::scale(glm::mat4{1}, glm::vec3{radius, radius, glm::distance(p_base, p_top)});


  if (useWireFrames) {
    // select shader program
    glUseProgram(this->shaderProgram_frame);

    // send mvp
    const auto mvp = PV * model;
    int mvpLoc = glGetUniformLocation(this->shaderProgram_frame, "mvp");
    glUniformMatrix4fv(mvpLoc, 1, GL_FALSE, glm::value_ptr(mvp));

    // send color
    int colorLoc = glGetUniformLocation(this->shaderProgram_frame, "color");
    glUniform3fv(colorLoc, 1, glm::value_ptr(color));


    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // disable blending
    glDisable(GL_BLEND);
  } else {
    // select shader program
    glUseProgram(this->shaderProgram_fill);

    // send normal matrix
    const auto normalMatrix = glm::mat3(glm::transpose(glm::inverse(model)));
    int normalMatrixLoc = glGetUniformLocation(this->shaderProgram_fill, "normalMatrix");
    glUniformMatrix3fv(normalMatrixLoc, 1, GL_FALSE, glm::value_ptr(normalMatrix));

    // send model
    int modelLoc = glGetUniformLocation(this->shaderProgram_fill, "model");
    glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

    // send projection * view
    int pvLoc = glGetUniformLocation(this->shaderProgram_fill, "pv");
    glUniformMatrix4fv(pvLoc, 1, GL_FALSE, glm::value_ptr(PV));

    // send  object color
    int objectColorLoc = glGetUniformLocation(this->shaderProgram_fill, "objectColor");
    glUniform3fv(objectColorLoc, 1, glm::value_ptr(color));

    // send light color
    int lightColorLoc = glGetUniformLocation(this->shaderProgram_fill, "lightColor");
    glUniform3fv(lightColorLoc, 1, glm::value_ptr(glm::vec3{0.7, 0.7, 0.7}));

    // send light position
    int lightPosLoc = glGetUniformLocation(this->shaderProgram_fill, "lightPos");
    glUniform3fv(lightPosLoc, 1, glm::value_ptr(camPos));


    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  // draw
  glDrawElements(GL_TRIANGLES, 3 * this->numPhis, GL_UNSIGNED_INT, 0);
  // draw
  glDrawElements(GL_TRIANGLE_FAN, this->numPhis + 2, GL_UNSIGNED_INT, (void*)((3 * this->numPhis) * sizeof(GLuint)));


  // unbind VAO
  glBindVertexArray(0);
}