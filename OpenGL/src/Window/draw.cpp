#include "Window.h"

using namespace std;


void Window::drawSphere (
  const glm::vec3 & pos,
  const float       radius,
  const glm::vec3 & color,
  const bool        useWireFrames
) {
  this->sphere->draw(this->camera.pos, this->VP, pos, radius, color, useWireFrames);
}


void Window::drawCylinder (
  const glm::vec3 & p1,
  const glm::vec3 & p2,
  const float       radius,
  const glm::vec3 & color,
  const bool        useWireFrames
) {
  this->cylinder->draw(this->camera.pos, this->VP, p1, p2, radius, color, useWireFrames);
}


void Window::drawRectangle (
  const glm::vec3 & pos,
  const glm::vec3 & dim,
  const glm::vec3 & color,
  const bool        useWireFrames
) {
  this->rectangle->draw(this->camera.pos, this->VP, pos, dim, color, useWireFrames);
}


void Window::drawLine (
  const glm::vec3 & p1,
  const glm::vec3 & p2,
  const glm::vec3 & color,
  const float       width
) {
  this->line->draw(this->VP, p1, p2, color, width);
}


void Window::drawLine_2D (
  const glm::vec2 & p1,
  const glm::vec2 & p2,
  const glm::vec3 & color,
  const float       width
) {
  ////
  //// construct orthogonal projection
  ////
  glm::mat4 projection = glm::ortho(0., (double) this->width, 0., (double) this->height, 0., 1.);

  this->line_2D->draw(projection, p1, p2, color, width);
}


void Window::drawRectangle_2D (
  const glm::vec2 & p0,
  const float       width,
  const float       height,
  const glm::vec3 & color,
  const float       angle
) {
  ////
  //// construct orthogonal projection
  ////
  glm::mat4 projection = glm::ortho(0., (double) this->width, 0., (double) this->height, 0., 1.);

  this->rectangle_2D->draw(projection, p0, width, height, angle, color);
}


void Window::drawArrowHead_2D (
  const glm::vec2 & p0,
  const float       length,
  const glm::vec3 & color,
  const float       rotAngle,
  const float       shearAngle
) {
  ////
  //// construct orthogonal projection
  ////
  glm::mat4 projection = glm::ortho(0., (double) this->width, 0., (double) this->height, 0., 1.);

  this->arrowHead_2D->draw(projection, p0, length, rotAngle, shearAngle, color);
}


void Window::drawArrow_2D (
  const glm::vec2 & beg,
  const glm::vec2 & end,
  const glm::vec3 & color,
  const float       tailWidth,
  const float       headLength
) {
  const float rotAngle = atan2(beg.x - end.x, end.y - beg.y);
  const float length   = glm::distance(beg, end);

  // make line tiny by longer to ensure overlap
  // this->drawLine_2D(beg, beg + (length - headLength)*glm::normalize(end - beg), glm::vec3{1, 1, 1}, tailWidth);

  const float _rotAngle = rotAngle + M_PI/2;
  this->drawRectangle_2D(beg + 0.5f*tailWidth*glm::vec2{sin(_rotAngle), -cos(_rotAngle)},
                         length - headLength,
                         tailWidth,
                         color,
                         _rotAngle);

  this->drawArrowHead_2D(end, headLength, color, rotAngle, M_PI/2);
}


void Window::drawText(
  const string &  text,
  float           x,
  float           y,
  const float     scale,
  const glm::vec3 & color
) {
  this->text->draw(this->width, this->height, text, x, y, scale, color);
}
void Window::_drawText(
  const string &  text,
  float &         x,
  float &         y,
  const float     scale,
  const glm::vec3 & color
) {
  this->text->draw(this->width, this->height, text, x, y, scale, color);
}