#include "Rectangle_2D.h"

using namespace std;


Rectangle_2D::Rectangle_2D () {
  this->init();
}



void Rectangle_2D::init () {
  ////
  //// create vertices and indices
  ////
  vector<float> vertices{0, 0, 1, 0, 1, 1, 0, 1};
  vector<unsigned> indices{0, 1, 2, 0, 2, 3};


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
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);

  // unbind VOA
  glBindVertexArray(0);

  // obtain current directory of this cpp file
  size_t indexOfSlash = string{__FILE__}.find_last_of("/");
  string currentDir = string{__FILE__}.substr(0, indexOfSlash + 1);

  // load shaders
  this->shaderProgram = agp::util::loadShaders((currentDir + "../shaders/vertex_2D.glsl").c_str(), (currentDir + "../shaders/fragment.glsl").c_str());
}