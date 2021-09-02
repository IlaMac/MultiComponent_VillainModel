#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

out vec3 FragPos;
out vec3 Normal;

uniform mat4 model;
uniform mat4 pv;
uniform mat3 normalMatrix;

void main () {

  gl_Position = pv * model * vec4(aPos, 1.0);

  FragPos = vec3(model * vec4(aPos, 1.0));
  Normal = normalMatrix * aNormal;

}