#include "VisualOpenGL.h"
#include <Windows.h>


double lastX, lastY;
bool mousePressed = false;
// Global variables to hold rotation state
Eigen::Quaternionf rotation = Eigen::Quaternionf::Identity();


#define MAX_CHAR    128
GLuint TextFont;

//Ӣ�ġ�����
void XPrintString(const char* s)
{

	glListBase(TextFont);
	glCallLists(strlen(s), GL_UNSIGNED_BYTE, s);
}




//�������֣���֧�ֺ��֡�unicode
void initFontData()
{
	TextFont = glGenLists(MAX_CHAR);
	wglUseFontBitmaps(wglGetCurrentDC(), 0, MAX_CHAR, TextFont);
}



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
	glPushMatrix();  // ���浱ǰ��ģ����ͼ��ՁE
	glTranslatef(-length * 3, -length * 3, 0);  // ������E�ԭ���ƶ������ڵ����½?
	

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

	glPopMatrix();
}
void drawAxis1(float length, const Eigen::Matrix3f& rotationMatrix) {
	glPushMatrix();  // ���浱ǰ��ģ����ͼ��ՁE
	glTranslatef(-length * 3, -length * 3, 0);  // ������E�ԭ���ƶ������ڵ����½?
	Eigen::Matrix4f matrix4f = Eigen::Matrix4f::Identity(); // ����һ��4x4��λ����
	matrix4f.block<3, 3>(0, 0) = rotationMatrix;
	glMultMatrixf(matrix4f.data());//��ת������

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

	glPopMatrix();  // �ָ�֮ǰ�����ģ����ͼ��ՁE
}
float getRotationAngleZ(const Eigen::Matrix3f& rotationmatrix) {
	// ʹ�� atan2 ����������Ƕ�
	float angle = std::atan2(rotationmatrix(1, 0), rotationmatrix(0, 0));

	// ���Ƕȴӻ���ת��Ϊ�����������Ҫ�Ļ�
	angle = angle * (180.0 / 3.1415926535f);

	return angle; // ���ػ���ֵ
}


// Function to convert HSV to RGB
void hsvToRgb(float h, float s, float v, float& r, float& g, float& b) {
	int i = int(h / 60.0f) % 6;
	float f = (h / 60.0f) - i;
	float p = v * (1.0f - s);
	float q = v * (1.0f - f * s);
	float t = v * (1.0f - (1.0f - f) * s);

	switch (i) {
	case 0: r = v, g = t, b = p; break;
	case 1: r = q, g = v, b = p; break;
	case 2: r = p, g = v, b = t; break;
	case 3: r = p, g = q, b = v; break;
	case 4: r = t, g = p, b = v; break;
	case 5: r = v, g = p, b = q; break;
	}
}