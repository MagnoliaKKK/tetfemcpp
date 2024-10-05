#pragma once
#ifndef MOUSECALLBACKS_HPP
#define MOUSECALLBACKS_HPP

#include "GLFW/glfw3.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

// Function declarations for mouse button and cursor position callbacks
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
void cursorPosCallback(GLFWwindow* window, double xpos, double ypos);

// Function declaration for framebuffer size callback
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// Function declaration for scroll callback
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

// Global variables to store zoom factor, transformation matrix, and aspect ratio
extern float zoomFactor;
extern Eigen::Matrix4f transformationMatrix;
extern float aspectRatio;

#endif // MOUSECALLBACKS_HPP
