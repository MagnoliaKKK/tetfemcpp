#include <iostream>
#include <vector>
#include <cstring>  // for std::strcpy
#include "tetgen.h"  // Include the TetGen header file
#include <fstream>
#include "GLFW/glfw3.h"
#include <cmath>
#include "VisualOpenGL.h"
#include "ReadSTL.h"
#include "GroupDivision.h"

// Global variables to store zoom factor and transformation matrix
Eigen::Matrix4f transformationMatrix = Eigen::Matrix4f::Identity();
double youngs = 10000000;
double poisson = 0.49;
double density = 1000;

int main() {


	tetgenio in, out;
	in.firstnumber = 1;  // All indices start from 1
	readSTL("stls/cube.stl", in);

	// Configure TetGen behavior
	tetgenbehavior behavior;
	//char args[] = "pq1.414a0.1";
	char args[] = "pq1.1/15a0.1"; // pq1.414a0.1 minratio 1/ mindihedral -q maxvolume -a switches='pq1.1/15a0.003' "pq1.1/15a0.0005"
	behavior.parse_commandline(args);

	// Call TetGen to tetrahedralize the geometry
	tetrahedralize(&behavior, &in, &out);

	int groupNum = 3; //ｷﾖｼｸﾗ・ﾄｿﾇｰﾏﾈﾗ鋐・ﾗ・Objectﾀ犲ﾍﾑﾕﾉｫｶｼﾐｴﾋﾀﾁﾋ
	Object object;
	divideIntoGroups(out, object, groupNum);

	// Accessing and printing the groups and their tetrahedra
	for (int i = 0; i < groupNum; ++i) {  // Loop over the groups
		Group& group = object.getGroup(i);
		std::cout << "Group " << i << " has " << group.tetrahedra.size() << " tetrahedra." << std::endl;

		for (size_t j = 0; j < group.tetrahedra.size(); ++j) {  // Loop over the tetrahedra in each group
			Tetrahedron* tet = group.tetrahedra[j];
			//std::cout << "  Tetrahedron " << j << " vertices:" << std::endl;

			for (int k = 0; k < 4; ++k) {  // Loop over the vertices in each tetrahedron
				Vertex* vertex = tet->vertices[k];
				//std::cout << "    Vertex " << k << ": (" << vertex->x << ", " << vertex->y << ", " << vertex->z << ")" << std::endl;
			}
		}
	}


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

	Eigen::MatrixXd Ke;
	//Ke = object.groups[1].tetrahedra[0]->createElementK(youngs, poisson, );
	
	Eigen::Matrix4f mat;
	while (!glfwWindowShouldClose(window)) {

		/*object.groups[0].tetrahedra[1]->vertices[2]->x += 0.01;
		object.groups[1].getUniqueVertices()[3]->y += 0.001;*/

		//object.getGroup(1).tetrahedra[0]->vertices[0]->x += 0.01;
		// Render here
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		
		// ｻ贍ﾆﾗ・・
		drawAxis(0.3f);

		// Enable wireframe mode
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		mat = Eigen::Matrix4f::Identity();
		mat.block<3, 3>(0, 0) = rotation.toRotationMatrix(); 
		glMultMatrixf(mat.data());

		// Draw vertices
		// ﾉ靹ﾃｵ羞ﾄｴ｡ﾎｪ10ﾏﾘ
		glPointSize(5.0f);
		// ｻ贍ﾆｰﾗﾉｫｵ・
		glColor3f(1.0f, 1.0f, 1.0f);
		glBegin(GL_POINTS);
		for (int groupIdx = 0; groupIdx < 3; ++groupIdx) {
			Group& group = object.getGroup(groupIdx);
			std::vector<Vertex*> uniqueVertices = group.getUniqueVertices();
			for (Vertex* vertex : uniqueVertices) {
				glVertex3f(vertex->x, vertex->y, vertex->z);
			}
		}
		glEnd();


		// Draw edges
		glBegin(GL_LINES);




		//ｷﾖﾗ鰕ﾔﾊｾ
		for (int groupIdx = 0; groupIdx < 3; ++groupIdx) {
			Group& group = object.getGroup(groupIdx);
			for (Tetrahedron* tet : group.tetrahedra) {
				for (int edgeIdx = 0; edgeIdx < 6; ++edgeIdx) {  // Loop through each edge in the tetrahedron
					Edge* edge = tet->edges[edgeIdx];
					Vertex* vertex1 = edge->vertices[0];
					Vertex* vertex2 = edge->vertices[1];
					bool isSurfaceEdge = edge->isBoundary;

					float red = 0.0f, green = 0.0f, blue = 0.0f;

					// Assign color based on groupIdx
					if (groupIdx == 0) {
						red = 1.0f;  // Red for group 0
					}
					else if (groupIdx == 1) {
						green = 1.0f;  // Green for group 1
					}
					else if (groupIdx == 2) {
						blue = 1.0f;  // Blue for group 2
					}

					// If it's a boundary edge, you may want to adjust the color or keep as is
					// For example, make the color brighter if it's a boundary edge
					if (isSurfaceEdge) {
						red = std::min(1.0f, red + 0.5f);
						green = std::min(1.0f, green + 0.5f);
						blue = std::min(1.0f, blue + 0.5f);
					}

					drawEdge(vertex1, vertex2, red, green, blue);
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
