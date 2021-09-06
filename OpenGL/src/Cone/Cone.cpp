#include "Cone.h"

using namespace std;


Cone::Cone (
  const unsigned numPhis
) :
  numPhis{numPhis}
{
  this->init();
}



void Cone::init () {
  ////
  //// create vertices
  ////
  vector<float> vertices(2 * 3 * 3 * this->numPhis + 2 * 3 * (this->numPhis + 2));
  const float dPhi = 2 * M_PI / this->numPhis;

  // top
  for (unsigned i = 0; i < this->numPhis; i++) {

    // coordinates
    const glm::vec3 p0 = {cos(i * dPhi), sin(i * dPhi), 0};
    const glm::vec3 p1 = {cos((i + 1) * dPhi), sin((i + 1) * dPhi), 0};
    const glm::vec3 p2 = {0, 0, 1};

    // compute normal
    const glm::vec3 n0 = glm::vec3{cos((i + 0.0) * dPhi), sin((i + 0.0) * dPhi), 1} / (float) sqrt(2);
    const glm::vec3 n1 = glm::vec3{cos((i + 1.0) * dPhi), sin((i + 1.0) * dPhi), 1} / (float) sqrt(2);
    const glm::vec3 n2 = glm::vec3{cos((i + 0.5) * dPhi), sin((i + 0.5) * dPhi), 1} / (float) sqrt(2);

    // bottom
    vertices[i*18 + 0] = p0.x;
    vertices[i*18 + 1] = p0.y;
    vertices[i*18 + 2] = p0.z;

    vertices[i*18 + 3] = n0.x;
    vertices[i*18 + 4] = n0.y;
    vertices[i*18 + 5] = n0.z;

    // bottom
    vertices[i*18 + 6] = p1.x;
    vertices[i*18 + 7] = p1.y;
    vertices[i*18 + 8] = p1.z;

    vertices[i*18 + 9]  = n1.x;
    vertices[i*18 + 10] = n1.y;
    vertices[i*18 + 11] = n1.z;

    // top
    vertices[i*18 + 12] = p2.x;
    vertices[i*18 + 13] = p2.y;
    vertices[i*18 + 14] = p2.z;

    vertices[i*18 + 15] = n2.x;
    vertices[i*18 + 16] = n2.y;
    vertices[i*18 + 17] = n2.z;
  }


  // base
  const unsigned I = 2 * 3 * 3 * this->numPhis;
  for (unsigned i = 0; i < this->numPhis + 2; i++) {

    glm::vec3 p = {cos(i * dPhi), sin(i * dPhi), 0};
    if (i == 0) p = {0, 0, 0};

    vertices[I + i*6 + 0] = p.x;
    vertices[I + i*6 + 1] = p.y;
    vertices[I + i*6 + 2] = p.z;

    vertices[I + i*6 + 3] = 0;
    vertices[I + i*6 + 4] = 0;
    vertices[I + i*6 + 5] = -1;
  }


  ////
  //// create indices
  ////
  vector<unsigned> indices(3 * this->numPhis + this->numPhis + 2);
  for (unsigned i = 0; i < indices.size(); i++) indices[i] = i;



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
  string currentDir   = string{__FILE__}.substr(0, indexOfSlash + 1);

  // load shaders
  this->shaderProgram_fill  = agp::util::loadShaders((currentDir + "../shaders/vertex_diffuse.glsl").c_str(), (currentDir + "../shaders/fragment_diffuse.glsl").c_str());
  this->shaderProgram_frame = agp::util::loadShaders((currentDir + "../shaders/vertex_ambient.glsl").c_str(), (currentDir + "../shaders/fragment.glsl").c_str());
}