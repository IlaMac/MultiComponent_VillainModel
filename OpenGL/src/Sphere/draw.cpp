#include "Sphere.h"

using namespace std;


void Sphere::draw (
  const glm::vec3 & camPos,
  const glm::mat4 & PV,
  const glm::vec3 & pos,
  const float       radius,
  const glm::vec3 & color,
  const bool        useWireFrames
) {
  // bind VAO
  glBindVertexArray(this->VAO);


  // position sphere
  const auto model = glm::translate(glm::mat4{1}, pos)
                   * glm::scale(glm::mat4{1}, glm::vec3{radius, radius, radius});


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


  // depth
  glDepthMask(GL_TRUE);


  // draw
  glDrawElements(GL_TRIANGLES, 6 * (this->numThetas - 1) * this->numPhis, GL_UNSIGNED_INT, 0);


  // unbind VAO
  glBindVertexArray(0);
}