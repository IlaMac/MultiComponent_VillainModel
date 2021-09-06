#include "Rectangle.h"

using namespace std;


Rectangle::Rectangle () {
  this->init();
}



void Rectangle::init () {
  ////
  //// create vertices and indices
  ////
  vector<float> vertices{0,0,0, 0,0,-1,
                         1,0,0, 0,0,-1,
                         1,1,0, 0,0,-1,
                         0,1,0, 0,0,-1,

                         0,0,0, 0,-1,0,
                         1,0,0, 0,-1,0,
                         1,0,1, 0,-1,0,
                         0,0,1, 0,-1,0,

                         0,0,0, -1,0,0,
                         0,1,0, -1,0,0,
                         0,1,1, -1,0,0,
                         0,0,1, -1,0,0,

                         0,0,1, 0,0,1,
                         1,0,1, 0,0,1,
                         1,1,1, 0,0,1,
                         0,1,1, 0,0,1,

                         0,1,0, 0,1,0,
                         1,1,0, 0,1,0,
                         1,1,1, 0,1,0,
                         0,1,1, 0,1,0,

                         1,0,0, 1,0,0,
                         1,1,0, 1,0,0,
                         1,1,1, 1,0,0,
                         1,0,1, 1,0,0};
  vector<unsigned> indices{0,1,2,    0,2,3,
                           4,5,6,    4,6,7,
                           8,9,10,   8,10,11,
                           12,13,15, 13,14,15,
                           16,17,19, 17,18,19,
                           20,21,23, 21,22,23};


  // wire
  vector<float> vertices_wire{0,0,0,
                              1,0,0,
                              1,1,0,
                              0,1,0,

                              0,0,1,
                              1,0,1,
                              1,1,1,
                              0,1,1
                            };
  vector<unsigned> indices_wire{0,1, 1,2, 2,3, 3,0,
                                4,5, 5,6, 6,7, 7,4,
                                0,4, 1,5, 2,6, 3,7};


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




  // generate VBO and EBO
  GLuint VBO_wire;
  GLuint EBO_wire;
  glGenBuffers(1, &VBO_wire);
  glGenBuffers(1, &EBO_wire);

  // generate and bind the triangle VAO
  glGenVertexArrays(1, &this->VAO_wire);

  // bind and fill VOA with instructions
  glBindVertexArray(this->VAO_wire);

  // load data
  glBindBuffer(GL_ARRAY_BUFFER, VBO_wire);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertices_wire.size(), &vertices_wire[0], GL_STATIC_DRAW);
  // load  indices
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO_wire);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * indices_wire.size(), &indices_wire[0], GL_STATIC_DRAW);
  // data interpretation
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
  glEnableVertexAttribArray(0);
  // glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)(3 * sizeof(float)));
  // glEnableVertexAttribArray(1);

  // unbind VOA
  glBindVertexArray(0);







  // obtain current directory of this cpp file
  size_t indexOfSlash = string{__FILE__}.find_last_of("/");
  string currentDir = string{__FILE__}.substr(0, indexOfSlash + 1);

  // load shaders
  this->shaderProgram_fill  = agp::util::loadShaders((currentDir + "../shaders/vertex_diffuse.glsl").c_str(), (currentDir + "../shaders/fragment_diffuse.glsl").c_str());
  this->shaderProgram_frame = agp::util::loadShaders((currentDir + "../shaders/vertex_ambient.glsl").c_str(), (currentDir + "../shaders/fragment.glsl").c_str());
}