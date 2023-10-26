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
#include <unordered_map>


class Edge {
public:
	Vertex* vertices[2];
	bool isBoundary;

	Edge(Vertex* v1, Vertex* v2) : isBoundary(false) {
		vertices[0] = v1;
		vertices[1] = v2;
	}
};

class Tetrahedron {
public:
	Vertex* vertices[4];
	Edge* edges[6];  // Each tetrahedron has six edges

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
	std::unordered_map<int, Vertex*> verticesMap;  // Map to store unique vertices by index 哈希映射 同时也保证没有重复的点

	void addTetrahedron(Tetrahedron* tet) {
		tetrahedra.push_back(tet);
		for (int i = 0; i < 4; ++i) {
			verticesMap[tet->vertices[i]->index] = tet->vertices[i]; 
		}
	}

	std::vector<Vertex*> getUniqueVertices() {
		std::vector<Vertex*> uniqueVertices;
		for (auto& pair : verticesMap) {
			uniqueVertices.push_back(pair.second);
		}
		return uniqueVertices;
	}
	//用法 auto aaa = object.getGroup(0).getUniqueVertices();

};

class Object {
public:
	Group groups[3];

	Group& getGroup(int index) {
		return groups[index];
	}
};

std::unordered_set<std::string> boundaryEdgesSet;  // Set to store boundary edges

void findBoundaryEdges(tetgenio& out) {
	int indexOffset = out.firstnumber;  // Get the index offset (0 or 1)
	for (int i = 0; i < out.numberoftrifaces; ++i) {
		for (int j = 0; j < 3; ++j) {
			int vertexIndex1 = out.trifacelist[i * 3 + j] - indexOffset;
			int vertexIndex2 = out.trifacelist[i * 3 + ((j + 1) % 3)] - indexOffset;
			std::string edgeKey = vertexIndex1 < vertexIndex2 ?
				std::to_string(vertexIndex1) + "-" + std::to_string(vertexIndex2) :
				std::to_string(vertexIndex2) + "-" + std::to_string(vertexIndex1);
			boundaryEdgesSet.insert(edgeKey);
		}
	}
}

void divideIntoGroups(tetgenio& out, Object& object, int numGroups) {

	findBoundaryEdges(out);  // Populate the boundaryEdgesSet

	double minX = out.pointlist[0];
	double maxX = out.pointlist[0];

	// Find min and max X coordinates
	for (int i = 0; i < out.numberofpoints; ++i) {
		double x = out.pointlist[i * 3];
		if (x < minX) minX = x;
		if (x > maxX) maxX = x;
	}

	double range = maxX - minX;
	double groupRange = range / numGroups;

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
		if (groupIndex >= numGroups) groupIndex = numGroups - 1;  // Update here

		object.getGroup(groupIndex).addTetrahedron(tet);

		// Set up edges for each tetrahedron
		static int edgeIndices[6][2] = { {0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3} };
		for (int j = 0; j < 6; ++j) {
			Vertex* vertex1 = tet->vertices[edgeIndices[j][0]];
			Vertex* vertex2 = tet->vertices[edgeIndices[j][1]];
			Edge* edge = new Edge(vertex1, vertex2);

			std::string edgeKey = vertex1->index < vertex2->index ?
				std::to_string(vertex1->index) + "-" + std::to_string(vertex2->index) :
				std::to_string(vertex2->index) + "-" + std::to_string(vertex1->index);

			edge->isBoundary = boundaryEdgesSet.count(edgeKey) > 0;
			tet->edges[j] = edge;
		}
	}
}

// Global variables to store zoom factor and transformation matrix
Eigen::Matrix4f transformationMatrix = Eigen::Matrix4f::Identity();

int main() {


	tetgenio in, out;
	in.firstnumber = 1;  // All indices start from 1
	readSTL("stls/cubeLong.stl", in);

	// Configure TetGen behavior
	tetgenbehavior behavior;
	//char args[] = "pq1.414a0.1";
	char args[] = "pq1.1/15a0.0005"; // pq1.414a0.1 minratio 1/ mindihedral -q maxvolume -a switches='pq1.1/15a0.003'
	behavior.parse_commandline(args);

	// Call TetGen to tetrahedralize the geometry
	tetrahedralize(&behavior, &in, &out);

	int groupNum = 3; //分几组 目前先最多3组 Object类和颜色都写死了
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

	
	
	Eigen::Matrix4f mat;
	while (!glfwWindowShouldClose(window)) {
		//object.getGroup(1).tetrahedra[0]->vertices[0]->y += 0.01;
		// Render here
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		


		// Enable wireframe mode
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		mat = Eigen::Matrix4f::Identity();
		mat.block<3, 3>(0, 0) = rotation.toRotationMatrix(); 
		glMultMatrixf(mat.data());

		// Draw vertices
		// 设置点的大小为10像素
		glPointSize(5.0f);
		// 绘制白色点
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




		//分组显示
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
