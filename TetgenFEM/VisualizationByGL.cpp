#include <VisualizationByGL.h>

// Global variables definition

void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
    static bool mousePressed = false;  // Static variable to keep track of mouse state
    static double lastX = 0.0, lastY = 0.0; // Static variables to store the last cursor position

    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            mousePressed = true;
            glfwGetCursorPos(window, &lastX, &lastY);
        }
        else if (action == GLFW_RELEASE) {
            mousePressed = false;
        }
    }
}

Eigen::Quaternionf rotation = Eigen::Quaternionf::Identity();

void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
    static bool mousePressed = false;  // Static variable to keep track of mouse state
    static double lastX = 0.0, lastY = 0.0; // Static variables to store the last cursor position
   
   
   
    if (mousePressed) {
        float dx = (xpos - lastX) * 0.005f;
        float dy = (ypos - lastY) * 0.005f;
        Eigen::AngleAxisf aaX(dy, Eigen::Vector3f::UnitX());
        Eigen::AngleAxisf aaY(dx, Eigen::Vector3f::UnitY());
        // Assuming rotation is another global variable
        rotation = aaY * aaX * rotation;
        lastX = xpos;
        lastY = ypos;
    }
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
    aspectRatio = static_cast<float>(width) / static_cast<float>(height ? height : 1);  // Avoid division by zero
}
float zoomFactor = 1.0f;
Eigen::Matrix4f transformationMatrix = Eigen::Matrix4f::Identity();
float aspectRatio = 1.0f;
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    

    zoomFactor *= (1.0f + 0.1f * yoffset);  // Adjust zoom factor

    // Update projection matrix
    float left = -zoomFactor * aspectRatio;
    float right = zoomFactor * aspectRatio;
    float bottom = -zoomFactor;
    float top = zoomFactor;
    float nearVal = -1.0f;
    float farVal = 1.0f;

    Eigen::Matrix4f projectionMatrix;
    projectionMatrix <<
        2.0f / (right - left), 0.0f, 0.0f, -(right + left) / (right - left),
        0.0f, 2.0f / (top - bottom), 0.0f, -(top + bottom) / (top - bottom),
        0.0f, 0.0f, -2.0f / (farVal - nearVal), -(farVal + nearVal) / (farVal - nearVal),
        0.0f, 0.0f, 0.0f, 1.0f;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glMultMatrixf(projectionMatrix.data());
    glMatrixMode(GL_MODELVIEW);
}
