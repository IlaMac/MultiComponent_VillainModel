#ifndef TEXT_H
#define TEXT_H

#include <string>

#include "../common.h"
#include "../util/util.h"

// FreeType
#include <ft2build.h>
#include FT_FREETYPE_H

class Text {
  private:

    GLuint VAO, VBO, shaderProgram;

    /// Holds all state information relevant to a character as loaded using FreeType
    struct Character {
        GLuint TextureID;     // ID handle of the glyph texture
        glm::ivec2 Size;      // Size of glyph
        glm::ivec2 Bearing;   // Offset from baseline to left/top of glyph
        GLuint Advance;       // Horizontal offset to advance to next glyph
    };

    std::map<GLchar, Character> characters;

  public:

  private:
    void init ();


  public:
    Text ();

    void draw (const int,
               const int,
               const std::string &,
               float &,
               float &,
               const float,
               const glm::vec3 &);

};


#endif