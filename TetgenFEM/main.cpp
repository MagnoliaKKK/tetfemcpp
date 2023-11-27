#include <iostream>
#include <vector>
#include <cstring>  // for std::strcpy
#include "tetgen.h"  // Include the TetGen header file
#include <fstream>
#include "GL/glew.h" 
#include "GLFW/glfw3.h"

#include <cmath>
#include <random>
#include "VisualOpenGL.h"
#include "ReadSTL.h"
#include "GroupDivision.h"


//C:/Users/xu_yu/Desktop/tmp/arial.ttf

// Global variables to store zoom factor and transformation matrix
Eigen::Matrix4f transformationMatrix = Eigen::Matrix4f::Identity();
double youngs = 100000000;
double poisson = 0.49;
double density = 1000;




int main() {



	tetgenio in, out;
	in.firstnumber = 1;  // All indices start from 1
	readSTL("stls/cubeLong.stl", in);

	// Configure TetGen behavior
	tetgenbehavior behavior;
	//char args[] = "pq1.414a0.1";
	char args[] = "pq1.1/15a0.005"; // pq1.414a0.1 minratio 1/ mindihedral -q maxvolume -a switches='pq1.1/15a0.003' "pq1.1/15a0.0005"
	behavior.parse_commandline(args);

	// Call TetGen to tetrahedralize the geometry
	tetrahedralize(&behavior, &in, &out);

	int groupNum = 2; //Object类和颜色都写死了 不能超出class Object {里的组数
	Object object;
	divideIntoGroups(out, object, groupNum); //convert tetgen to our data structure
	object.updateIndices(); // 每个点分配一个独立index，重复的改新的index
	object.assignLocalIndicesToAllGroups(); //分配Local index
	object.generateUniqueVertices();//产生UniqueVertices

	
	// Accessing and printing the groups and their tetrahedra
	for (int i = 0; i < groupNum; ++i) {  // Loop over the groups
		Group& group = object.getGroup(i);
		std::cout << "Group " << i << " has " << group.tetrahedra.size() << " tetrahedra." << std::endl;
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
		object.commonPoints = object.findCommonVertices(object.groups[0], object.groups[1]);
		//object.commonPoints1 = object.findCommonVertices(object.groups[1], object.groups[2]);
		//固定点设置
		for (Group& g : object.groups) {
			// 遍历Group中的每个Vertex
			for (const auto& vertexPair : g.verticesMap) {
				// 对每个顶点调用setFixedIfBelowThreshold方法
				Vertex* vertex = vertexPair.second;
				
				vertex->setFixedIfBelowThreshold();
			}
			
		}


		
		//double aa = object.groups[0].tetrahedra[0]->calMassTetra(density);
		object.groups[0].calMassMatrix(density);
		object.groups[0].calMassGroup();
		object.groups[0].calDampingMatrix();
		object.groups[0].calCenterofMass();
		Eigen::Vector3d groupCOM = object.groups[0].centerofMass;
		object.groups[0].calGroupK(youngs, poisson);
		object.groups[0].setVertexMassesFromMassMatrix();
		object.groups[0].calMassDistributionMatrix();
		object.groups[0].calPrimeVec();
		object.groups[0].calInitCOM();
		object.groups[0].calLocalPos();
		object.groups[0].calRotationMatrix();
		object.groups[0].calLHS();
		//object.groups[0].updateVertexPositions();
		//object.groups[0].calRHS();
		//.groups[0].calDeltaX();
		//object.groups[0].calFbind(object.commonPoints.first, object.commonPoints.second, 1000);

		object.groups[1].calMassMatrix(density);
		object.groups[1].calMassGroup();
		object.groups[1].calDampingMatrix();
		object.groups[1].calCenterofMass();
		Eigen::Vector3d groupCOM1 = object.groups[1].centerofMass;
		object.groups[1].calGroupK(youngs, poisson);
		object.groups[1].setVertexMassesFromMassMatrix();
		object.groups[1].calMassDistributionMatrix();
		object.groups[1].calPrimeVec();
		object.groups[1].calInitCOM();
		object.groups[1].calLocalPos();
		object.groups[1].calRotationMatrix();
		object.groups[1].calLHS();
		//object.groups[1].calRHS();
		//object.groups[1].calDeltaX();
		//object.groups[1].updateVertexPositions();

		/*object.groups[2].calMassMatrix(density);
		object.groups[2].calMassGroup();
		object.groups[2].calDampingMatrix();
		object.groups[2].calCenterofMass();
		Eigen::Vector3d groupCOM2 = object.groups[2].centerofMass;
		object.groups[2].calGroupK(youngs, poisson);
		object.groups[2].setVertexMassesFromMassMatrix();
		object.groups[2].calMassDistributionMatrix();
		object.groups[2].calPrimeVec();
		object.groups[2].calInitCOM();
		object.groups[2].calLocalPos();
		object.groups[2].calRotationMatrix();
		object.groups[2].calLHS();*/
		//object.groups[0].calRHS();
		//object.groups[2].calDeltaX();
		//object.groups[2].updateVertexPositions();

		object.PBDLOOP(2);

		
		
		// Render here
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		// 绘制坐眮E丒
		drawAxis(0.3f);

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


		for (int groupIdx = 0; groupIdx < 3; ++groupIdx) { //写点的标号，画字
			Group& group = object.getGroup(groupIdx);

			//画不重复的版本
			//std::vector<Vertex*> uniqueVertices = group.getUniqueVertices();
			//for (size_t i = 0; i < uniqueVertices.size(); ++i) {
			//	Vertex* vertex = uniqueVertices[i];
			//	char buffer[5]; // 分配足够大的缓冲区
			//	sprintf_s(buffer, "%d", vertex->index); // 将int转换为char*
			//	glColor3f(1, 0.0f, 0.0f);
			//	glRasterPos3f(vertex->x, vertex->y, vertex->z);
			//	XPrintString(buffer);
			//}


			//画重复的版本


			for (Tetrahedron* tetra : group.tetrahedra) { // 遍历组中的每个四面体
				for (int i = 0; i < 4; ++i) { // 遍历四面体的每个顶点
					Vertex* vertex = tetra->vertices[i];
					char buffer[5]; // 分配足够大的缓冲区
					sprintf_s(buffer, "%d", vertex->index); // 将int转换为char*
					//if (groupIdx == 0) {
					//	glColor3f(1, 0.0f, 0.0f);
					//	glRasterPos3f(vertex->x, vertex->y, vertex->z);
					//	XPrintString(buffer);
					//}
						
					//if(groupIdx == 1)
					//{
					//	glColor3f(0, 1, 0.0f);
					//	glRasterPos3f(vertex->x, vertex->y, vertex->z);
					//	XPrintString(buffer);
					//}
					//	

					//std::default_random_engine generator(vertex->index);//随机数发生器，用于字符偏移防重叠
					//std::uniform_real_distribution<float> distribution(0, 0.05);
					//float random_number = distribution(generator);
					glColor3f(1, 0.0f, 0.0f);
					glRasterPos3f(vertex->x + 0, vertex->y + 0, vertex->z + 0);
					XPrintString(buffer);
					
				}
			}

		}

		
		// Draw edges
		glBegin(GL_LINES);




		//分组显示
		for (int groupIdx = 0; groupIdx < 3; ++groupIdx) { //这里也是写死了
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
