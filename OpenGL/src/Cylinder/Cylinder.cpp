#include "Cylinder.h"

using namespace std;


Cylinder::Cylinder (
  const unsigned numPhis,
  const unsigned numStacks
) :
  numPhis{numPhis},
  numStacks{numStacks}
{
  this->init();
}



void Cylinder::init () {
  ////
  //// create vertices
  ////
  vector<float> vertices(2 * 3 * (2 + (this->numStacks + 3) * this->numPhis));

  const float dPhi = 2 * M_PI / this->numPhis,
              dZ = 1 / (float) this->numStacks;

  // add first and last
  vertices[0] = 0;
  vertices[1] = 0;
  vertices[2] = 0;
  vertices[3] = 0;
  vertices[4] = 0;
  vertices[5] = -1;
  vertices[vertices.size() - 2 * 3 + 0] = 0;
  vertices[vertices.size() - 2 * 3 + 1] = 0;
  vertices[vertices.size() - 2 * 3 + 2] = 1;
  vertices[vertices.size() - 2 * 3 + 3] = 0;
  vertices[vertices.size() - 2 * 3 + 4] = 0;
  vertices[vertices.size() - 2 * 3 + 5] = 1;


  // the ringlike ones
  for (unsigned i = 0; i < this->numStacks + 3; i++) {

    // take into account that the two first and last rings have same z coordinate
    const unsigned _i = min(max(0, ((int) i) - 1), (int) this->numStacks);

    for (unsigned j = 0; j < this->numPhis; j++) {
      unsigned index = 2 * 3 * (i * this->numPhis + j + 1);

      // position
      vertices[index + 0] = cos(j * dPhi);
      vertices[index + 1] = sin(j * dPhi);
      vertices[index + 2] = _i * dZ;

      // normal vectors
      if (i == 0) {
        vertices[index + 3] = 0;
        vertices[index + 4] = 0;
        vertices[index + 5] = -1;
      } else if (i == this->numStacks + 2) {
        vertices[index + 3] = 0;
        vertices[index + 4] = 0;
        vertices[index + 5] = 1;
      } else {
        vertices[index + 3] = cos(j * dPhi);
        vertices[index + 4] = sin(j * dPhi);
        vertices[index + 5] = 0;
      }
    }
  }


  ////
  //// create indices
  ////
  vector<unsigned> indices(2 * 3 * this->numPhis + 6 * this->numStacks * this->numPhis);

  // bottom
  for (unsigned j = 0; j < this->numPhis; j++) {
    const unsigned index = 3 * j;

    indices[index + 0] = 0;
    indices[index + 1] = j + 1;
    indices[index + 2] = (j + 1 == this->numPhis) ? 1 : j + 2;
  }

  // top
  for (unsigned j = 0; j < this->numPhis; j++) {
    const unsigned index = 3 * this->numPhis + 6 * this->numStacks * this->numPhis + 3 * j;
    const unsigned o = 1 + (this->numStacks + 2) * this->numPhis;

    indices[index + 0] = vertices.size() / 6 - 1;
    indices[index + 1] = o + j;
    indices[index + 2] = o + ((j == this->numPhis - 1) ? 0 : j + 1);
  }

  // sides
  for (unsigned i = 0; i < this->numStacks; i++) {
    for (unsigned j = 0; j < this->numPhis; j++) {
      const unsigned index = 3 * this->numPhis + 6 * i * this->numPhis + 6 * j;
      const unsigned o1 = 1 + (i + 1) * this->numPhis;
      const unsigned o2 = 1 + (i + 2) * this->numPhis;

      indices[index + 0] = o1 + j;
      indices[index + 1] = o2 + j;
      indices[index + 2] = o2 + ((j == this->numPhis - 1) ? 0 : j + 1);

      indices[index + 3] = o2 + ((j == this->numPhis - 1) ? 0 : j + 1);
      indices[index + 4] = o1 + ((j == this->numPhis - 1) ? 0 : j + 1);
      indices[index + 5] = o1 + j;
    }
  }

  ////
  //// load data and instructions
  ////

  // generate VBO and EBO
  GLuint VBO;
  GLuint EBO;
  glGenBuffers(1, &VBO);
  glGenBuffers(1, &EBO);

  // generate and bind the triangle VAO
  glGenVertexArrays(1, &this->VAO);

  // bind and fill VOA with instructions
  glBindVertexArray(this->VAO);

  // load data
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), &vertices[0], GL_STATIC_DRAW);
  // load  indices
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * indices.size(), &indices[0], GL_STATIC_DRAW);
  // data interpretation
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
  glEnableVertexAttribArray(1);

  // unbind VOA
  glBindVertexArray(0);

  // obtain current directory of this cpp file
  size_t indexOfSlash = string{__FILE__}.find_last_of("/");
  string currentDir = string{__FILE__}.substr(0, indexOfSlash + 1);

  // load shaders
  this->shaderProgram_fill = agp::util::loadShaders((currentDir + "../shaders/vertex_diffuse.glsl").c_str(), (currentDir + "../shaders/fragment_diffuse.glsl").c_str());
  this->shaderProgram_frame = agp::util::loadShaders((currentDir + "../shaders/vertex_ambient.glsl").c_str(), (currentDir + "../shaders/fragment.glsl").c_str());
}