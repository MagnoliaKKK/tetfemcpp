#include "VisualOpenGL.h"


double lastX, lastY;
bool mousePressed = false;
// Global variables to hold rotation state
Eigen::Quaternionf rotation = Eigen::Quaternionf::Identity();
// Callback to handle mouse button events
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
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


// Callback to handle mouse motion events
void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
	if (mousePressed) {
		float dx = (xpos - lastX) * 0.005f;
		float dy = (ypos - lastY) * 0.005f;
		Eigen::AngleAxisf aaX(dy, Eigen::Vector3f::UnitX());
		Eigen::AngleAxisf aaY(dx, Eigen::Vector3f::UnitY());
		rotation = aaY * aaX * rotation;
		lastX = xpos;
		lastY = ypos;
	}
}

float aspectRatio = 1.0f;
float zoomFactor = 1.0f;

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
	aspectRatio = static_cast<float>(width) / static_cast<float>(height ? height : 1);  // Avoid division by zero
}
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

std::string createEdgeId(Vertex* vertex1, Vertex* vertex2) {
	std::ostringstream idStream;
	if (vertex1 < vertex2) {
		idStream << vertex1 << "_" << vertex2;
	}
	else {
		idStream << vertex2 << "_" << vertex1;
	}
	return idStream.str();
}



// Function to draw a line between two vertices with a specified color
void drawEdge(Vertex* vertex1, Vertex* vertex2, float r, float g, float b) {
	glColor3f(r, g, b);
	glVertex3f(
		vertex1->x,
		vertex1->y,
		vertex1->z
	);
	glVertex3f(
		vertex2->x,
		vertex2->y,
		vertex2->z
	);
}

void drawAxis(float length) {
	glPushMatrix();  // 保存当前的模型视图矩諄E
	glTranslatef(-length * 3, -length * 3, 0);  // 将坐眮E嵩阋贫酱翱诘挠蚁陆?

	glBegin(GL_LINES);
	// X axis in red
	glColor3f(1.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(length, 0.0f, 0.0f);
	// Y axis in green
	glColor3f(0.0f, 1.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, length, 0.0f);
	// Z axis in blue
	glColor3f(0.0f, 0.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(0.0f, 0.0f, length);
	glEnd();

	glPopMatrix();  // 恢复之前保存的模型视图矩諄E
}
