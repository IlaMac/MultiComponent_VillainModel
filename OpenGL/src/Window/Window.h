#ifndef WINDOW_H
#define WINDOW_H

#include <vector>
#include <utility>

#include "../common.h"
#include "../Camera/Camera.h"
#include "../Sphere/Sphere.h"
#include "../Cylinder/Cylinder.h"
#include "../Line/Line.h"
#include "../Text/Text.h"
#include "../Line_2D/Line_2D.h"
#include "../Rectangle_2D/Rectangle_2D.h"
#include "../ArrowHead_2D/ArrowHead_2D.h"

class Window {
  private:
    GLFWwindow * glfwWindow;

    // projection * view matrix
    glm::mat4 VP;

    // current width and height
    int width, height;


    ////
    //// callback functions
    ////
    std::vector<std::pair<void (*) (void *, const int , const int, const int), void *> > onKeyPressCallbacks;
    std::vector<std::pair<void (*) (void *, const int , const int, const int), void *> > onKeyRepeatCallbacks;

    unsigned long frameNum = 0;

    float framerate;


  public:
    Camera camera{&width, &height, {0, 0, -5}, {0, 0, 1}, {0, 1, 0}};

    Sphere *       sphere;
    Cylinder *     cylinder;
    Line *         line;
    Text *         text;
    Line_2D *      line_2D;
    Rectangle_2D * rectangle_2D;
    ArrowHead_2D * arrowHead_2D;

  private:

    void init ();

    ////
    //// user defined events
    ////
    virtual void onKeyPress (const int, const int, const int);
    virtual void onKeyRelease (const int, const int, const int);
    virtual void onKeyRepeat (const int, const int, const int);
    virtual void onMouseMove (const double &, const double &);
    virtual void onMousePress (const double &, const double &, const int, const int);
    virtual void onMouseRelease (const double &, const double &, const int, const int);
    virtual void onScroll (const double &, const double &);


    ////
    //// glfw events
    ////
    static void onKey (GLFWwindow *, const int, const int, const int, const int);
    static void onMouseMove (GLFWwindow *, const double, const double);
    static void onMouseButton (GLFWwindow *, const int, const int, const int);
    static void onScroll (GLFWwindow *, const double, const double);


  public:
    Window (const int,
            const int,
            const std::string &);

    virtual ~Window ();


    ////
    //// main loop
    ////
    void show (bool (*) (void *, const double &),
               void *);

    unsigned long getFrameNumber () const;

    float getFramerate () const;

    float getWidth () const;
    float getHeight () const;


    ////
    //// additional event listeners
    ////
    void listenOnKeyPress (void (*) (void *, const int , const int, const int),
                           void *);
    void listenOnKeyRepeat (void (*) (void *, const int , const int, const int),
                            void *);


    ////
    //// 3D primitives
    ////
    void drawSphere (const glm::vec3 &,
                     const float,
                     const glm::vec3 &,
                     const bool = false);

    void drawCylinder (const glm::vec3 &,
                       const glm::vec3 &,
                       const float,
                       const glm::vec3 &,
                       const bool = false);

    void drawLine (const glm::vec3 &,
                   const glm::vec3 &,
                   const glm::vec3 &,
                   const float = 1);

    void drawText (const std::string &,
                   float,
                   float,
                   const float,
                   const glm::vec3 &);
    void _drawText (const std::string &,
                   float &,
                   float &,
                   const float,
                   const glm::vec3 &);


    void drawLine_2D (const glm::vec2 &,
                      const glm::vec2 &,
                      const glm::vec3 &,
                      const float = 1);

    void drawRectangle_2D (const glm::vec2 &,
                           const float,
                           const float,
                           const glm::vec3 &,
                           const float = 0);

    void drawArrowHead_2D (const glm::vec2 &,
                           const float,
                           const glm::vec3 &,
                           const float = 0,
                           const float = M_PI/2);

    void drawArrow_2D (const glm::vec2 &,
                       const glm::vec2 &,
                       const glm::vec3 &,
                       const float = 1,
                       const float = 10);


    ////
    ////
    ////
    static void GLFWwindowclosefun (GLFWwindow * window) { std::cout << "window closed!" << std::endl; };

};


#endif