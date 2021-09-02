#include "Sphere.h"

using namespace std;


Sphere::Sphere (
  const unsigned numThetas,
  const unsigned numPhis
) :
  numThetas{numThetas},
  numPhis{numPhis}
{
  this->init();
}



void Sphere::init () {
  ////
  //// create vertices and indices
  ////
  vector<unsigned> indices(6 * (this->numThetas - 1) * this->numPhis);
  vector<float> vertices(6 * this->numThetas * this->numPhis);

  const float dTheta = M_PI / (this->numThetas - 1),
              dPhi   = 2 * M_PI / this->numPhis;

  for (unsigned i = 0; i < this->numThetas; i++) {
    for (unsigned j = 0; j < this->numPhis; j++) {
      unsigned index = 6*(this->numPhis*i + j);
      // position
      vertices[index + 0] = sin(i * dTheta) * cos(j * dPhi);
      vertices[index + 2] = sin(i * dTheta) * sin(j * dPhi);
      vertices[index + 1] = cos(i * dTheta);

      // normal
      vertices[index + 3] = sin(i * dTheta) * cos(j * dPhi);
      vertices[index + 5] = sin(i * dTheta) * sin(j * dPhi);
      vertices[index + 4] = cos(i * dTheta);

      if (i + 1 < this->numThetas) {
        unsigned index = this->numPhis*i + j;

        // first triangle
        indices[6*index + 0] = index;
        indices[6*index + 1] = index + this->numPhis;
        if (j + 1 < this->numPhis) {
          indices[6*index + 2] = index + this->numPhis + 1;
        } else {
          indices[6*index + 2] = index + 1;
        }

        // second triangle
        indices[6*index + 3] = index;
        if (j + 1 < this->numPhis) {
          indices[6*index + 4] = index + this->numPhis + 1;
          indices[6*index + 5] = index + 1;
        } else {
          indices[6*index + 4] = index + 1;
          indices[6*index + 5] = index + 1- this->numPhis;
        }
      }
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