#include <iostream>
#include <vector>
#include <cstring>  // for std::strcpy
#include "tetgen.h"  // Include the TetGen header file
#include <fstream>
#include "GLFW/glfw3.h"
#include <cmath>
#include <unordered_set>
#include "VisualOpenGL.h"
#include "ReadSTL.h"


// Global variables to store zoom factor and transformation matrix

Eigen::Matrix4f transformationMatrix = Eigen::Matrix4f::Identity();


// Function to create a unique identifier for an edge


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
	readSTL("C:/tetFEMCpp/bunny.stl", in);

	// Configure TetGen behavior
	tetgenbehavior behavior;
	char args[] = "pq1.1/15a0.0005"; // pq1.414a0.1 minratio 1/ mindihedral -q maxvolume -a switches='pq1.1/15a0.003'
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
