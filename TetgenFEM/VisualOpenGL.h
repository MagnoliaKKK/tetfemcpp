#pragma once
#include "GLFW/glfw3.h"
#include "tetgen.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
extern Eigen::Quaternionf rotation;
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
std::string createEdgeId(int vertexIndex1, int vertexIndex2);
void drawEdge(tetgenio& out, int vertexIndex1, int vertexIndex2, float r, float g, float b);