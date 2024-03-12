
//#define EIGEN_USE_MKL_ALL
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
float youngs = 100000;
float poisson = 0.49;
float density = 1000;
int groupNum, groupNumX =1, groupNumY = 1, groupNumZ = 1;//Object�����ɫ��д���� ���ܳ���class Object {�������
int wKey = 0;


int main() {



	tetgenio in, out;
	in.firstnumber = 1;  // All indices start from 1
	readSTL("stls/cubeLong.stl", in);

	// Configure TetGen behavior
	tetgenbehavior behavior;
	//char args[] = "pq1.414a0.1";
	char args[] = "pq1.15a0.0005";  // pq1.414a0.1 minratio 1/ mindihedral -q maxvolume -a switches='pq1.1/15a0.003' "pq1.1/15a0.0005"
	behavior.parse_commandline(args);

	// Call TetGen to tetrahedralize the geometry
	tetrahedralize(&behavior, &in, &out);

	
	Object object;
	groupNum = groupNumX * groupNumY * groupNumZ;
	object.groupNum = groupNum;
	object.groupNumX = groupNumX;
	object.groupNumY = groupNumY;
	object.groupNumZ = groupNumZ;
	divideIntoGroups(out, object, groupNumX, groupNumY, groupNumZ); //convert tetgen to our data structure
	object.updateIndices(); // ÿ�������һ������index���ظ��ĸ��µ�index
	object.assignLocalIndicesToAllGroups(); //����Local index
	object.generateUniqueVertices();//����UniqueVertices
	
	object.updateAdjacentGroupIndices(groupNumX, groupNumY, groupNumZ);
	for (int i = 0; i < groupNum; ++i) {
		// ��ÿ������� storeAdjacentGroupsCommonVertices ����
		object.storeAdjacentGroupsCommonVertices(i);
	}
	
	// Accessing and printing the groups and their tetrahedra
//#pragma omp parallel for
	for (int i = 0; i < groupNum; ++i) {  // Loop over the groups
		Group& group = object.getGroup(i);
		std::cout << "Group " << i << " has " << group.tetrahedra.size() << " tetrahedra." << std::endl;
		group.LHS_I = Eigen::MatrixXf::Identity(3 * group.verticesMap.size(), 3 * group.verticesMap.size()); //��ʡʱ��С����
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
	//object.findCommonVertices();
	//object.commonPoints = object.findCommonVertices1(object.groups[0], object.groups[1]);
	//object.commonPoints1 = object.findCommonVertices1(object.groups[1], object.groups[2]);
	//object.commonPoints2 = object.findCommonVertices1(object.groups[2], object.groups[3]);
	//object.commonPoints3 = object.findCommonVertices1(object.groups[3], object.groups[4]);
	//std::pair<std::vector<Vertex*>, std::vector<Vertex*>> commonVertices2 = object.findCommonVertices1(object.groups[0], object.groups[1]);
	for (Group& g : object.groups) {
		 //����Group�е�ÿ��Vertex
		for (const auto& vertexPair : g.verticesMap) {
			// ��ÿ���������setFixedIfBelowThreshold����
			Vertex* vertex = vertexPair.second;

			vertex->setFixedIfBelowThreshold();
		}

	}
	
	
	
	/////////��ͷ���̶���
	//float maxY = -std::numeric_limits<float>::infinity(); // ��ʼ��Ϊ��Сֵ
	//Vertex* vertexWithMaxY = nullptr;
	//// ���� y ֵ���Ķ���
	//for (Group& g : object.groups) {
	//	for (const auto& vertexPair : g.verticesMap) {
	//		Vertex* vertex = vertexPair.second;
	//		if (vertex->y > maxY) {
	//			maxY = vertex->y;
	//			vertexWithMaxY = vertex;
	//		}
	//	}
	//}
	//// �� y ֵ���Ķ�����Ϊ�̶���
	//if (vertexWithMaxY != nullptr) {
	//	vertexWithMaxY->isFixed = true; // ������һ������ setFixed �����ö���Ĺ̶�״̬
	//}
	/////////
	
#pragma omp parallel for
	for (int i = 0; i < object.groupNum; ++i) {
		object.groups[i].calMassMatrix(density);
		object.groups[i].calDampingMatrix();
		object.groups[i].calCenterofMass();
		object.groups[i].calInitCOM();
		object.groups[i].calLocalPos(); // �����ʼλ�����ʼ���ĵĲ�ֵ
		object.groups[i].calGroupK(youngs, poisson);
		
		object.groups[i].setVertexMassesFromMassMatrix();
		object.groups[i].calMassGroup();
		object.groups[i].calMassDistributionMatrix();
		//object.groups[i].inverseTerm = (object.groups[i].massMatrix + object.groups[i].dampingMatrix * 0.01f).inverse(); //˳˳����������
		//object.groups[i].inverseTermSparse = object.groups[i].inverseTerm.sparseView();
		object.groups[i].calLHS();
	}

	//for calculate frame rate
	double lastTime = glfwGetTime();
	int nbFrames = 0;
	glfwSwapInterval(0); //��ֱͬ��

	
	while (!glfwWindowShouldClose(window)) {
		
		//object.commonPoints1 = object.findCommonVertices(object.groups[1], object.groups[2]);
		//�̶�������
	
	
		if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
			wKey = 1;
			
		}
		else if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
			wKey = 2;

		}
		else if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
			wKey = 3;

		}
		else if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
			wKey = 4;

		}
		else
			wKey = 0;
		//std::cout << wKey << std::endl;// �� W ������ʱ���߼�
		//double aa = object.groups[0].tetrahedra[0]->calMassTetra(density);
		
		#pragma omp parallel for
		for (int i = 0; i < groupNum; i++) {
			//object.groups[i].calGroupKFEM(youngs, poisson);
			object.groups[i].calPrimeVec(wKey);
			/*object.groups[i].calLHSFEM();
			object.groups[i].calRHSFEM();
			object.groups[i].calDeltaXFEM();
			object.groups[i].calculateCurrentPositionsFEM();
			object.groups[i].updateVelocityFEM();
			object.groups[i].updatePositionFEM();*/
			
			
			object.groups[i].calRotationMatrix();
		
		}
	
		object.PBDLOOP(2);


		
		
		// Render here
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		
		//drawAxis1(0.3f, object.groups[0].rotate_matrix);
		// ��������E�E
		drawAxis(0.3f);
		//std::cout << getRotationAngleZ(object.groups[0].rotate_matrix) << std::endl;;
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
		for (int groupIdx = 0; groupIdx < groupNum; ++groupIdx) {
			Group& group = object.getGroup(groupIdx);
			std::vector<Vertex*> uniqueVertices = group.getUniqueVertices();
			for (Vertex* vertex : uniqueVertices) {
				glVertex3f(vertex->x, vertex->y, vertex->z);


			}
		}
		glEnd();


		for (int groupIdx = 0; groupIdx < groupNum; ++groupIdx) { //д��ı�ţ�����
			//Group& group = object.getGroup(groupIdx);

			////�����ظ��İ汾
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


			//for (Tetrahedron* tetra : group.tetrahedra) { // �������е�ÿ��������
			//	for (int i = 0; i < 4; ++i) { // �����������ÿ������
			//		Vertex* vertex = tetra->vertices[i];
			//		char buffer[5]; // �����㹻��Ļ�����
			//		sprintf_s(buffer, "%d", vertex->index); // ��intת��Ϊchar*
			//		//if (groupIdx == 0) {
			//		//	glColor3f(1, 0.0f, 0.0f);
			//		//	glRasterPos3f(vertex->x, vertex->y, vertex->z);
			//		//	XPrintString(buffer);
			//		//}
			//			
			//		//if(groupIdx == 1)
			//		//{
			//		//	glColor3f(0, 1, 0.0f);
			//		//	glRasterPos3f(vertex->x, vertex->y, vertex->z);
			//		//	XPrintString(buffer);
			//		//}
			//		//	

			//		//std::default_random_engine generator(vertex->index);//������������������ַ�ƫ�Ʒ��ص�
			//		//std::uniform_real_distribution<float> distribution(0, 0.05);
			//		//float random_number = distribution(generator);
			//		glColor3f(1, 0.0f, 0.0f);
			//		glRasterPos3f(vertex->x + 0, vertex->y + 0, vertex->z + 0);
			//		XPrintString(buffer);
			//		
			//	}
			//}

		}
		//for (int groupIdx = 0; groupIdx < groupNum; ++groupIdx) {
		//	Group& group = object.getGroup(groupIdx);

		//	// ʹ�������Ѿ�����õ�����
		//	Eigen::Vector3f& center = group.centerofMass;

		//	// ��������ת��Ϊ�ַ���
		//	char groupNumber[10];
		//	sprintf_s(groupNumber, "%d", groupIdx);

		//	// Ϊ����������ɫ
		//	glColor3f(1.0f, 1.0f, 0.0f); // ��ɫ�����ı�

		//	// �������ŵ�λ�ò�����
		//	glRasterPos3f(center[0], center[1], center[2]);
		//	XPrintString(groupNumber);

		//	// ... [���Ʊ�Ե�Ͷ�����������] ...
		//}
		
		// Draw edges
		glBegin(GL_LINES);


		

		//������ʾ
		for (int groupIdx = 0; groupIdx < groupNum; ++groupIdx) {
			Group& group = object.getGroup(groupIdx);
			for (Tetrahedron* tet : group.tetrahedra) {
				for (int edgeIdx = 0; edgeIdx < 6; ++edgeIdx) {  // Loop through each edge in the tetrahedron
					Edge* edge = tet->edges[edgeIdx];
					Vertex* vertex1 = edge->vertices[0];
					Vertex* vertex2 = edge->vertices[1];
					bool isSurfaceEdge = edge->isBoundary;

					// Use HSV to RGB conversion to create a unique color for each group
					float hue = (360.0f * groupIdx) / groupNum;  // Distribute hues evenly across the spectrum
					float saturation = 1.0f;  // Full saturation
					float value = 1.0f;      // Full brightness

					// Convert HSV to RGB
					float red, green, blue;
					hsvToRgb(hue, saturation, value, red, green, blue);

					// If it's a boundary edge, you may want to adjust the color or keep as is
					// For example, make the color brighter if it's a boundary edge
					if (isSurfaceEdge) {
						red = std::min(1.0f, red + 0.3f);
						green = std::min(1.0f, green + 0.3f);
						blue = std::min(1.0f, blue + 0.3f);
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

		//calculate frame rate
		double currentTime = glfwGetTime();
		nbFrames++;
		if (currentTime - lastTime >= 1.0) { // �����ȥ������1��
			printf("%d frames/sec\n", nbFrames);
			nbFrames = 0;
			lastTime += 1.0;
		}

	}

	glfwTerminate();

	return 0;
}
