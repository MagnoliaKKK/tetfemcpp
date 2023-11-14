#include <iostream>
#include <vector>
#include <cstring>  // for std::strcpy
#include "tetgen.h"  // Include the TetGen header file
#include <fstream>
#include "GL/glew.h" 
#include "GLFW/glfw3.h"

#include <cmath>

#include "VisualOpenGL.h"
#include "ReadSTL.h"
#include "GroupDivision.h"


//C:/Users/xu_yu/Desktop/tmp/arial.ttf

// Global variables to store zoom factor and transformation matrix
Eigen::Matrix4f transformationMatrix = Eigen::Matrix4f::Identity();
double youngs = 10000000;
double poisson = 0.49;
double density = 1000;




int main() {



	tetgenio in, out;
	in.firstnumber = 1;  // All indices start from 1
	readSTL("stls/cubeLong.stl", in);

	// Configure TetGen behavior
	tetgenbehavior behavior;
	//char args[] = "pq1.414a0.1";
	char args[] = "pq1.1/15a0.1"; // pq1.414a0.1 minratio 1/ mindihedral -q maxvolume -a switches='pq1.1/15a0.003' "pq1.1/15a0.0005"
	behavior.parse_commandline(args);

	// Call TetGen to tetrahedralize the geometry
	tetrahedralize(&behavior, &in, &out);

	int groupNum = 2; //Object�����ɫ��д����
	Object object;
	divideIntoGroups(out, object, groupNum); //convert tetgen to our data structure
	//�������for��Object Group Vertex��������
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
	initFontData();
	while (!glfwWindowShouldClose(window)) {
		
		//double aa = object.groups[0].tetrahedra[0]->calMassTetra(density);
		//object.groups[0].calMassMatrix(density);
		//object.groups[0].calMassGroup();
		//object.groups[0].calDampingMatrix();
		//object.groups[0].calCenterofMass();
		//Eigen::Vector3d groupCOM = object.groups[0].centerofMass;
		//Ke = object.groups[0].tetrahedra[0]->createElementK(youngs, poisson, groupCOM);
		//object.groups[0].calGroupK(youngs, poisson);
		//object.groups[0].setVertexMassesFromMassMatrix();
		//object.groups[0].calMassDistributionMatrix();
		//object.groups[0].calPrimeVec();
		//object.groups[0].calInitCOM();
		//object.groups[0].calLocalPos();
		//object.groups[0].calRotationMatrix();
		//object.groups[0].calLHS();
		//object.groups[0].calRHS();
		
		//auto aaaaa = object.findCommonVertices(object.getGroup(0), object.getGroup(1));
		//object.groups[0].calDeltaX();

		/*object.groups[0].tetrahedra[1]->vertices[2]->x += 0.01;
		object.groups[1].getUniqueVertices()[3]->y += 0.001;*/

		//object.getGroup(1).tetrahedra[0]->vertices[0]->x += 0.01;
		// Render here
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		// ��������E�E
		drawAxis(0.3f);

		// Enable wireframe mode
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		mat = Eigen::Matrix4f::Identity();
		mat.block<3, 3>(0, 0) = rotation.toRotationMatrix();
		glMultMatrixf(mat.data());
		
		
		// Draw vertices
		// ���õ�Ĵ�СΪ10����
		glPointSize(5.0f);
		// ���ư�ɫ��
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


		for (int groupIdx = 0; groupIdx < 3; ++groupIdx) { //д��ı�ţ�����
			Group& group = object.getGroup(groupIdx);

			//�����ظ��İ汾
			//std::vector<Vertex*> uniqueVertices = group.getUniqueVertices();
			//for (size_t i = 0; i < uniqueVertices.size(); ++i) {
			//	Vertex* vertex = uniqueVertices[i];
			//	char buffer[5]; // �����㹻��Ļ�����
			//	sprintf_s(buffer, "%d", vertex->index); // ��intת��Ϊchar*
			//	glColor3f(1, 0.0f, 0.0f);
			//	glRasterPos3f(vertex->x, vertex->y, vertex->z);
			//	XPrintString(buffer);
			//}


			//���ظ��İ汾

			for (const auto& pair : group.verticesMap) {
				int key = pair.first;          // ����������
				Vertex* value = pair.second;   // ֵ��Vertex*��

				char buffer[5]; // �����㹻��Ļ�����
				sprintf_s(buffer, "%d", value->index); // ��intת��Ϊchar*
				glColor3f(1, 0.0f, 0.0f);
				glRasterPos3f(value->x, value->y, value->z);
				XPrintString(buffer);
			}

		}

		
		// Draw edges
		glBegin(GL_LINES);




		//������ʾ
		for (int groupIdx = 0; groupIdx < 3; ++groupIdx) { //����Ҳ��д����
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
