#pragma once
#include "GLFW/glfw3.h"
#include "tetgen.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

class Vertex {
public:
    double x, y, z, id;

    Vertex(double x, double y, double z, double id) : x(x), y(y), z(z), id(id) {}
};

extern Eigen::Quaternionf rotation;
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
std::string createEdgeId(Vertex* vertex1, Vertex* vertex2);
void drawEdge(Vertex* vertex1, Vertex* vertex2, float r, float g, float b);
