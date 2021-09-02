#ifndef CAMERA_H
#define CAMERA_H


#include "../common.h"


class Camera {
  private:

    // field of view (vertical)
    const float verticalFOV = glm::radians(45.0f);

    const float slow = 0.1,
                fast = 1;

    const int * width;
    const int * height;


    ////
    //// position quantities
    ////
    bool forward    = false,
         backward   = false,
         leftwards  = false,
         rightwards = false,
         upward     = false,
         downward   = false,
         boost      = false;


    ////
    //// orientation quantities
    ////
    bool orientation = false;
    double x_pre, y_pre;
    // int width, height;
    glm::vec3 front_pre, up_pre;


    ////
    //// initial position and camera orientation
    ////
    glm::vec3 initPos, initFront, initUp;


  public:
    glm::vec3 pos, front, up;

    float getHorizontalFOV ();
    float getVerticalFOV ();

  private:

  public:
    Camera (const int *,
            const int *,
            const glm::vec3 &,
            const glm::vec3 &,
            const glm::vec3 &);

    void update ();

    void activateForward ();
    void deactivateForward ();
    void activateBackward ();
    void deactivateBackward ();
    void activateLeftwards ();
    void deactivateLeftwards ();
    void activateRightwards ();
    void deactivateRightwards ();
    void activateUpward ();
    void deactivateUpward ();
    void activateDownward ();
    void deactivateDownward ();
    void activateBoost ();
    void deactivateBoost ();

    void activateOrientation (const double &, const double &); //, const int, const int);
    void deactivateOrientation ();
    void orient (const double &, const double &);

    void roll (const double &, const double &, const double &, const int, const int);

    void reset ();

    void setInitialPosition (const glm::vec3 &);
    void setInitialFront    (const glm::vec3 &);
    void setInitialUp       (const glm::vec3 &);
};


#endif