#include <iostream>
#include <vector>
#include <cstring>  // for std::strcpy
#include "tetgen.h"  // Include the TetGen header file
#include <fstream>
#include "GLFW/glfw3.h"
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <unordered_set>

// Global variables to hold rotation state
Eigen::Quaternionf rotation = Eigen::Quaternionf::Identity();
double lastX, lastY;
bool mousePressed = false;

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
float zoomFactor = 1.0f;
Eigen::Matrix4f transformationMatrix = Eigen::Matrix4f::Identity();
float aspectRatio = 1.0f;

void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
	glViewport(0, 0, width, height);
	aspectRatio = static_cast<float>(width) / static_cast<float>(height ? height : 1);  // Avoid division by zero
}

// Global variables to store zoom factor and transformation matrix


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

// Function to read STL file
void readSTL(const std::string& filename, tetgenio& in) {
	std::ifstream file(filename, std::ios::binary);
	if (!file) {
		std::cerr << "Error: Could not open file: " << filename << std::endl;
		return;
	}

	char header[80];
	file.read(header, 80);  // Read header

	uint32_t numTriangles;
	file.read(reinterpret_cast<char*>(&numTriangles), 4);  // Read number of triangles

	in.numberofpoints = numTriangles * 3;  // Each triangle has 3 vertices
	in.pointlist = new REAL[in.numberofpoints * 3];
	in.numberoffacets = numTriangles;
	in.facetlist = new tetgenio::facet[in.numberoffacets];
	in.facetmarkerlist = new int[in.numberoffacets];

	for (uint32_t i = 0; i < numTriangles; ++i) {
		float normal[3];
		file.read(reinterpret_cast<char*>(normal), 12);  // Read normal vector

		tetgenio::facet& f = in.facetlist[i];
		f.numberofpolygons = 1;
		f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
		f.numberofholes = 0;
		f.holelist = NULL;
		tetgenio::polygon& p = f.polygonlist[0];
		p.numberofvertices = 3;
		p.vertexlist = new int[p.numberofvertices];

		for (int j = 0; j < 3; ++j) {  // Read vertices
			float vertex[3];
			file.read(reinterpret_cast<char*>(vertex), 12);
			int index = i * 3 + j;
			in.pointlist[index * 3] = vertex[0];
			in.pointlist[index * 3 + 1] = vertex[1];
			in.pointlist[index * 3 + 2] = vertex[2];
			p.vertexlist[j] = index + 1;  // Indices in TetGen start from 1
		}

		uint16_t attributeByteCount;
		file.read(reinterpret_cast<char*>(&attributeByteCount), 2);  // Read attribute byte count

		in.facetmarkerlist[i] = 0;  // No facet markers
	}

	file.close();
}

// Function to create a unique identifier for an edge
std::string createEdgeId(int vertexIndex1, int vertexIndex2) {
	return vertexIndex1 < vertexIndex2 ?
		std::to_string(vertexIndex1) + "_" + std::to_string(vertexIndex2) :
		std::to_string(vertexIndex2) + "_" + std::to_string(vertexIndex1);
}

// Function to draw a line between two vertices with a specified color
void drawEdge(tetgenio& out, int vertexIndex1, int vertexIndex2, float r, float g, float b) {
	glColor3f(r, g, b);
	glVertex3f(
		out.pointlist[vertexIndex1 * 3],
		out.pointlist[vertexIndex1 * 3 + 1],
		out.pointlist[vertexIndex1 * 3 + 2]
	);
	glVertex3f(
		out.pointlist[vertexIndex2 * 3],
		out.pointlist[vertexIndex2 * 3 + 1],
		out.pointlist[vertexIndex2 * 3 + 2]
	);
}





int main() {

	// Initialize the GLFW library
	if (!glfwInit()) {
		return -1;
	}

	// Create a windowed mode window and its OpenGL context
	GLFWwindow* window = glfwCreateWindow(1080, 1080, "Tetrahedral Mesh Visualization", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}

	// Make the window's context current
	glfwMakeContextCurrent(window);
	// Set scroll callback
	glfwSetScrollCallback(window, scroll_callback);
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	glfwSetMouseButtonCallback(window, mouseButtonCallback);
	glfwSetCursorPosCallback(window, cursorPosCallback);

	tetgenio in, out;
	in.firstnumber = 1;  // All indices start from 1
	readSTL("C:/Users/XYX/Documents/VSProgramming/bunny.stl", in);

	// Configure TetGen behavior
	tetgenbehavior behavior;
	char args[] = "pq1.414a0.1";
	behavior.parse_commandline(args);

	// Call TetGen to tetrahedralize the geometry
	tetrahedralize(&behavior, &in, &out);


	// Inside your main function and before the render loop:
	std::unordered_set<std::string> surfaceEdges;
	for (int i = 0; i < out.numberoftrifaces; ++i) {
		for (int j = 0; j < 3; ++j) {
			int vertexIndex1 = out.trifacelist[i * 3 + j] - 1;
			int vertexIndex2 = out.trifacelist[i * 3 + (j + 1) % 3] - 1;
			surfaceEdges.insert(createEdgeId(vertexIndex1, vertexIndex2));
		}
	}

	while (!glfwWindowShouldClose(window)) {

		// Render here
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		


		// Enable wireframe mode
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		Eigen::Matrix4f mat = Eigen::Matrix4f::Identity();
		mat.block<3, 3>(0, 0) = rotation.toRotationMatrix();
		glMultMatrixf(mat.data());

		// Draw vertices
		glPointSize(5.0f);  // Set point size
		glPushMatrix();
		glMultMatrixf(transformationMatrix.data());
		glBegin(GL_POINTS);
		for (int i = 0; i < out.numberofpoints; ++i) {
			glVertex3f(
				out.pointlist[i * 3],
				out.pointlist[i * 3 + 1],
				out.pointlist[i * 3 + 2]
			);
		}
		glEnd();

		// Draw edges
		glBegin(GL_LINES);
		for (int i = 0; i < out.numberoftetrahedra; ++i) {
			int vertexIndices[4];
			for (int j = 0; j < 4; ++j) {
				vertexIndices[j] = out.tetrahedronlist[i * 4 + j] - 1;  // Indices in TetGen start from 1
			}
			// Draw all 6 edges of the tetrahedron
			for (int j = 0; j < 4; ++j) {
				for (int k = j + 1; k < 4; ++k) {
					bool isSurfaceEdge = surfaceEdges.count(createEdgeId(vertexIndices[j], vertexIndices[k])) > 0;
					drawEdge(out, vertexIndices[j], vertexIndices[k], isSurfaceEdge ? 1.0f : 1.0f, isSurfaceEdge ? 1.0f : 0.0f, isSurfaceEdge ? 1.0f : 0.0f);
				}
			}
		}
		glEnd();

		glPopMatrix();

		// Swap front and back buffers
		glfwSwapBuffers(window);

		// Poll for and process events
		glfwPollEvents();

	}

	glfwTerminate();

	return 0;
}
