#include "Text.h"

using namespace std;


void Text::draw (
  const int         width,
  const int         height,
  const string &    text,
  float &           x,
  float &           y,
  const float       scale,
  const glm::vec3 & color
) {
  // select shader program
  glUseProgram(this->shaderProgram);

  // send text color
  int textColorLoc = glGetUniformLocation(this->shaderProgram, "textColor");
  glUniform3f(textColorLoc, color.x, color.y, color.z);

  // send projection matrix
  glm::mat4 projection = glm::ortho(0.0f, static_cast<GLfloat>(width), 0.0f, static_cast<GLfloat>(height), 0.f, 1.f);
  int projectionLoc = glGetUniformLocation(this->shaderProgram, "projection");
  glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));

  // for some reason this must be active, otherwise there will only be triangles in certain cases
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // enable blending
  glEnable(GL_BLEND);


  glActiveTexture(GL_TEXTURE0);
  glBindVertexArray(this->VAO);

  // Iterate through all characters
  std::string::const_iterator c;
  for (c = text.begin(); c != text.end(); c++) {
    auto ch = this->characters[*c];

    GLfloat xpos = x + ch.Bearing.x * scale;
    GLfloat ypos = y - (ch.Size.y - ch.Bearing.y) * scale;

    GLfloat w = ch.Size.x * scale;
    GLfloat h = ch.Size.y * scale;

    // Update VBO for each character
    GLfloat vertices[6][4] = {
      { xpos,     ypos + h,   0.0, 0.0 },
      { xpos,     ypos,       0.0, 1.0 },
      { xpos + w, ypos,       1.0, 1.0 },

      { xpos,     ypos + h,   0.0, 0.0 },
      { xpos + w, ypos,       1.0, 1.0 },
      { xpos + w, ypos + h,   1.0, 0.0 }
    };

    // Render glyph texture over quad
    glBindTexture(GL_TEXTURE_2D, ch.TextureID);

    // Update content of VBO memory
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // Be sure to use glBufferSubData and not glBufferData

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    // Render quad
    glDrawArrays(GL_TRIANGLES, 0, 6);
    // Now advance cursors for next glyph (note that advance is number of 1/64 pixels)
    x += (ch.Advance >> 6) * scale; // Bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
  }

  glBindVertexArray(0);
  glBindTexture(GL_TEXTURE_2D, 0);
}