#include "Camera.h"

using namespace std;


Camera::Camera (
  const int * width,
  const int * height,
  const glm::vec3 & initPos,
  const glm::vec3 & initFront,
  const glm::vec3 & initUp
) :
  width{width},
  height{height},
  initPos{initPos},
  initFront{initFront},
  initUp{initUp},
  pos{initPos},
  front{initFront},
  up{initUp}
{ }


float Camera::getVerticalFOV () { return this->verticalFOV; }
float Camera::getHorizontalFOV () {
  // http://paulbourke.net/miscellaneous/lens/
  return 2 * atan(*(this->width) * tan(0.5 * this->verticalFOV) / *(this->height));
}


void Camera::activateForward      () { this->forward    = true;  }
void Camera::deactivateForward    () { this->forward    = false; }
void Camera::activateBackward     () { this->backward   = true;  }
void Camera::deactivateBackward   () { this->backward   = false; }
void Camera::activateLeftwards    () { this->leftwards  = true;  }
void Camera::deactivateLeftwards  () { this->leftwards  = false; }
void Camera::activateRightwards   () { this->rightwards = true;  }
void Camera::deactivateRightwards () { this->rightwards = false; }
void Camera::activateUpward       () { this->upward     = true;  }
void Camera::deactivateUpward     () { this->upward     = false; }
void Camera::activateDownward     () { this->downward   = true;  }
void Camera::deactivateDownward   () { this->downward   = false; }
void Camera::activateBoost        () { this->boost      = true;  }
void Camera::deactivateBoost      () { this->boost      = false; }

void Camera::update () {


  if (
    this->forward != this->backward ||
    this->rightwards != this->leftwards ||
    this->upward != this->downward
  ) {

    auto right = glm::cross(this->front, this->up);

    glm::vec3 dir = this->front * (float) (this->forward    - this->backward)
                  +       right * (float) (this->rightwards - this->leftwards)
                  + this->up    * (float) (this->upward     - this->downward);

    this->pos += glm::normalize(dir) * (this->boost ? this->fast : this->slow);
  }
}

void Camera::reset () {
  this->pos   = this->initPos;
  this->front = this->initFront;
  this->up    = this->initUp;
}


void Camera::activateOrientation (
  const double & x,
  const double & y
  // const int width,
  // const int height
) {
  this->orientation = true;
  this->x_pre = x;
  this->y_pre = y;
  // this->width  = width;
  // *(this->height) = height;
  this->front_pre = this->front;
  this->up_pre    = this->up;
}


void Camera::deactivateOrientation () {
  this->orientation = false;
}


void Camera::orient (
  const double & x,
  const double & y
) {
  if (this->orientation) {
    // pan
    float phi = 2 * atan( (x - x_pre) / (float) *(this->height) * tan(0.5 * this->verticalFOV));
    this->front = glm::normalize(glm::rotate(this->front_pre, phi, this->up_pre));

    // tilt
    float theta = 2 * atan( (y - y_pre) / (float) *(this->height) * tan(0.5 * this->verticalFOV));
    auto right_pre = glm::cross(this->front_pre, this->up_pre);
    this->front = glm::normalize(glm::rotate(this->front, theta, right_pre));
    this->up = glm::normalize(glm::rotate(this->up_pre, theta, right_pre));
  }
}


void Camera::roll (
  const double & scrollDir,
  const double & x,
  const double & y,
  const int width,
  const int height
) {
  if ( ! this->orientation) {
    // disabled during orientation

    // angle to cursor from each axis
    float phi   = 2 * atan( (0.5 * width - x)   / (float) height * tan(0.5 * this->verticalFOV));
    float theta = 2 * atan( ( 0.5 * height - y) / (float) height * tan(0.5 * this->verticalFOV));

    // compute right axis
    glm::vec3 right = glm::cross(this->front, this->up);

    // plane tangents
    glm::vec3 tang_front_up    = glm::rotate(this->front, theta, right);
    glm::vec3 tang_front_right = glm::rotate(this->front, phi, this->up);

    // normal vectors
    glm::vec3 norm_front_up    = glm::cross(right, tang_front_up);
    glm::vec3 norm_front_right = glm::cross(tang_front_right, up);

    // direction
    glm::vec3 dir = glm::cross(norm_front_up, norm_front_right);


    float alpha = glm::radians(scrollDir * 2);
    this->up = glm::normalize(glm::rotate(this->up, alpha, dir));
    this->front = glm::normalize(glm::rotate(this->front, alpha, dir));
  }
}




void Camera::setInitialPosition (const glm::vec3 & vec) { this->initPos   = vec; }
void Camera::setInitialFront    (const glm::vec3 & vec) { this->initFront = vec; }
void Camera::setInitialUp       (const glm::vec3 & vec) { this->initUp    = vec; }