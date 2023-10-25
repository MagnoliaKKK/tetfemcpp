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


class Tetrahedron {
public:
	Vertex* vertices[4];

	Tetrahedron(Vertex* v1, Vertex* v2, Vertex* v3, Vertex* v4) {
		vertices[0] = v1;
		vertices[1] = v2;
		vertices[2] = v3;
		vertices[3] = v4;
	}
};

class Group {
public:
	std::vector<Tetrahedron*> tetrahedra;

	void addTetrahedron(Tetrahedron* tet) {
		tetrahedra.push_back(tet);
	}
};

class Object {
public:
	Group groups[3];

	Group& getGroup(int index) {
		return groups[index];
	}
};

void divideIntoGroups(tetgenio& out, Object& object) {
	double minX = out.pointlist[0];
	double maxX = out.pointlist[0];

	// Find min and max X coordinates
	for (int i = 0; i < out.numberofpoints; ++i) {
		double x = out.pointlist[i * 3];
		if (x < minX) minX = x;
		if (x > maxX) maxX = x;
	}

	double range = maxX - minX;
	double groupRange = range / 3;

	// Create vertices
	std::vector<Vertex*> vertices;
	for (int i = 0; i < out.numberofpoints; ++i) {
		double x = out.pointlist[i * 3];
		double y = out.pointlist[i * 3 + 1];
		double z = out.pointlist[i * 3 + 2];
		vertices.push_back(new Vertex(x, y, z, i));
	}

	// Create tetrahedra and assign to groups
	for (int i = 0; i < out.numberoftetrahedra; ++i) {
		Vertex* v1 = vertices[out.tetrahedronlist[i * 4] - 1];
		Vertex* v2 = vertices[out.tetrahedronlist[i * 4 + 1] - 1];
		Vertex* v3 = vertices[out.tetrahedronlist[i * 4 + 2] - 1];
		Vertex* v4 = vertices[out.tetrahedronlist[i * 4 + 3] - 1];

		Tetrahedron* tet = new Tetrahedron(v1, v2, v3, v4);

		// Determine group based on average X coordinate
		double avgX = (v1->x + v2->x + v3->x + v4->x) / 4;
		int groupIndex = (avgX - minX) / groupRange;
		if (groupIndex > 2) groupIndex = 2;  // Ensure groupIndex is within bounds

		object.getGroup(groupIndex).addTetrahedron(tet);
	}
}

// Global variables to store zoom factor and transformation matrix
Eigen::Matrix4f transformationMatrix = Eigen::Matrix4f::Identity();

int main() {


	tetgenio in, out;
	in.firstnumber = 1;  // All indices start from 1
	readSTL("D:/docs/tetfemcpp/cube.stl", in);

	// Configure TetGen behavior
	tetgenbehavior behavior;
	//char args[] = "pq1.414a0.1";
	char args[] = "pq1.1/15a0.0005"; // pq1.414a0.1 minratio 1/ mindihedral -q maxvolume -a switches='pq1.1/15a0.003'
	behavior.parse_commandline(args);

	// Call TetGen to tetrahedralize the geometry
	tetrahedralize(&behavior, &in, &out);

	Object object;
	divideIntoGroups(out, object);

	// Accessing and printing the groups and their tetrahedra
	for (int i = 0; i < 3; ++i) {  // Loop over the groups
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

	std::unordered_set<std::string> surfaceEdges;

	// Iterate through all groups and tetrahedra to collect surface edges
	for (int groupIdx = 0; groupIdx < 3; ++groupIdx) {
		Group& group = object.getGroup(groupIdx);
		for (Tetrahedron* tet : group.tetrahedra) {
			for (int vertexIdx1 = 0; vertexIdx1 < 4; ++vertexIdx1) {
				for (int vertexIdx2 = vertexIdx1 + 1; vertexIdx2 < 4; ++vertexIdx2) {
					Vertex* vertex1 = tet->vertices[vertexIdx1];
					Vertex* vertex2 = tet->vertices[vertexIdx2];
					surfaceEdges.insert(createEdgeId(vertex1, vertex2));
				}
			}
		}
	}

	while (!glfwWindowShouldClose(window)) {
		object.getGroup(0).tetrahedra[0]->vertices[0]->x += 0.01;
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
		glBegin(GL_POINTS);
		for (int groupIdx = 0; groupIdx < 3; ++groupIdx) {
			Group& group = object.getGroup(groupIdx);
			for (Tetrahedron* tet : group.tetrahedra) {
				for (int vertexIdx = 0; vertexIdx < 4; ++vertexIdx) {
					Vertex* vertex = tet->vertices[vertexIdx];
					glVertex3f(vertex->x, vertex->y, vertex->z);
				}
			}
		}
		glEnd();

		// Draw edges
		glBegin(GL_LINES);
		for (int groupIdx = 0; groupIdx < 3; ++groupIdx) {
			Group& group = object.getGroup(groupIdx);
			for (Tetrahedron* tet : group.tetrahedra) {
				for (int vertexIdx1 = 0; vertexIdx1 < 4; ++vertexIdx1) {
					for (int vertexIdx2 = vertexIdx1 + 1; vertexIdx2 < 4; ++vertexIdx2) {
						Vertex* vertex1 = tet->vertices[vertexIdx1];
						Vertex* vertex2 = tet->vertices[vertexIdx2];
						bool isSurfaceEdge = surfaceEdges.count(createEdgeId(vertex1, vertex2)) > 0;
						drawEdge(vertex1, vertex2, isSurfaceEdge ? 1.0f : 1.0f, isSurfaceEdge ? 1.0f : 0.0f, isSurfaceEdge ? 1.0f : 0.0f);
					}
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
