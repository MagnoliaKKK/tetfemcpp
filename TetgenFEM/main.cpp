
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
float poisson = 0.45;
float density = 1000;
int groupNum, groupNumX =4, groupNumY = 1, groupNumZ =1;//Objectﾀ犲ﾍﾑﾕﾉｫｶｼﾐｴﾋﾀﾁﾋ ｲｻﾄﾜｳｬｳlass Object {ﾀ・ﾄﾗ鯡?
int wKey = 0;


void saveOBJ(const std::string& filename, std::vector<Group>& groups) {
	std::ofstream objFile(filename);
	if (!objFile.is_open()) {
		std::cerr << "Failed to open file for writing.\n";
		return;
	}

	std::unordered_map<Vertex*, int> vertexIndexMap;
	int currentIndex = 1;

	// 遍历组，找出所有边界边并记录其顶点
	for (const auto& group : groups) {
		for (const auto* tet : group.tetrahedra) {
			for (const auto* edge : tet->edges) {
				if (edge->isBoundary) {
					for (Vertex* vertex : edge->vertices) {
						if (vertexIndexMap.find(vertex) == vertexIndexMap.end()) {
							vertexIndexMap[vertex] = currentIndex++;
							objFile << "v " << vertex->x << " " << vertex->y << " " << vertex->z << "\n";
						}
					}
				}
			}
		}
	}

	// 再次遍历，这次是为了构建面
	for (const auto& group : groups) {
		for (const auto* tet : group.tetrahedra) {
			for (const auto* edge : tet->edges) {
				if (edge->isBoundary) {
					objFile << "f";
					for (Vertex* vertex : edge->vertices) {
						objFile << " " << vertexIndexMap[vertex];
					}
					objFile << "\n";
				}
			}
		}
	}

	objFile.close();
	std::cout << "OBJ file saved: " << filename << "\n";
}



int main() {



	tetgenio in, out;
	in.firstnumber = 1;  // All indices start from 1
	readSTL("stls/cubeLong.stl", in);
	//readOBJ("C:/Users/76739/Downloads/VegaFEM-v4.0/VegaFEM-v4.0/models/beam3/bunnyHDLow.obj", in);
	// Configure TetGen behavior
	tetgenbehavior behavior;
	//char args[] = "pq1.414a0.1";
	char args[] = "pq1.414a0.001";  // pq1.414a0.1 minratio 1/ mindihedral -q maxvolume -a switches='pq1.1/15a0.003' "pq1.1/15a0.0005 pq1.15a0.0001"
	behavior.parse_commandline(args);

	//char argsNode[] = "C:/Users/76739/Desktop/tetfemcpp/TetgenFEM/cubeThin";
	//char argsEle[] = "C:/Users/76739/Desktop/tetfemcpp/TetgenFEM//cubeThin";
	//if (!in.load_node(argsNode)) {
	//    std::cerr << "Error loading .node file!" << std::endl;
	//    return 1;
	//}

	////// Load the ele file
	//if (!in.load_tet(argsEle)) {
	//    std::cerr << "Error loading .ele file!" << std::endl;
	//    return 1;
	//}

	// Call TetGen to tetrahedralize the geometry
	tetrahedralize(&behavior, &in, &out);
	out.save_nodes("bunnyHDLow");
	out.save_elements("bunnyHDLow");

	//out = in;

	Object object;
	groupNum = groupNumX * groupNumY * groupNumZ;
	object.groupNum = groupNum;
	object.groupNumX = groupNumX;
	object.groupNumY = groupNumY;
	object.groupNumZ = groupNumZ;
	divideIntoGroups(out, object, groupNumX, groupNumY, groupNumZ); //convert tetgen to our data structure
	object.updateIndices(); // ﾃｿｸ羚ﾖﾅ萪ｻｸﾀﾁ｢index｣ｬﾖﾘｸｴｵﾄｸﾄﾐﾂｵﾄindex
	object.assignLocalIndicesToAllGroups(); //ｷﾖﾅ膈ocal index
	object.generateUniqueVertices();//ｲ揵￤niqueVertices
	
	object.updateAdjacentGroupIndices(groupNumX, groupNumY, groupNumZ);
	for (int i = 0; i < groupNum; ++i) {
		// ｶﾔﾃｿｸ魴ﾃ storeAdjacentGroupsCommonVertices ｺｯﾊ
		object.storeAdjacentGroupsCommonVertices(i);
	}
	
	// Accessing and printing the groups and their tetrahedra
//#pragma omp parallel for
	for (int i = 0; i < groupNum; ++i) {  // Loop over the groups
		Group& group = object.getGroup(i);
		std::cout << "Group " << i << " has " << group.tetrahedra.size() << " tetrahedra." << std::endl;
		group.LHS_I = Eigen::MatrixXf::Identity(3 * group.verticesMap.size(), 3 * group.verticesMap.size()); //ｽﾚﾊ｡ﾊｱｼ菻｡ﾄﾜﾊﾖ
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
		 //ｱ鯊ⅷroupﾖﾐｵﾄﾃｿｸertex
		for (const auto& vertexPair : g.verticesMap) {
			// ｶﾔﾃｿｸ･ｵ羞ﾃsetFixedIfBelowThresholdｷｽｷｨ
			Vertex* vertex = vertexPair.second;

			vertex->setFixedIfBelowThreshold();
		}

	}
	
	
	
	/////////ｾｾﾍｷｷ｢ｹﾌｶｨｷｨ
	//float maxY = -std::numeric_limits<float>::infinity(); // ｳｼｻｯﾎｪｼｫﾐ｡ﾖｵ
	//Vertex* vertexWithMaxY = nullptr;
	//// ｲ鰈ﾒ y ﾖｵﾗ鋗ﾄｶ･ｵ・
	//for (Group& g : object.groups) {
	//	for (const auto& vertexPair : g.verticesMap) {
	//		Vertex* vertex = vertexPair.second;
	//		if (vertex->y > maxY) {
	//			maxY = vertex->y;
	//			vertexWithMaxY = vertex;
	//		}
	//	}
	//}
	//// ｽｫ y ﾖｵﾗ鋗ﾄｶ･ｵ翹靜ｪｹﾌｶｨｵ・
	//if (vertexWithMaxY != nullptr) {
	//	vertexWithMaxY->isFixed = true; // ｼﾙﾉ勒ﾐﾒｻｸｽｷｨ setFixed ﾀｴﾉ靹ﾃｶ･ｵ羞ﾄｹﾌｶｨﾗｴﾌｬ
	//}
	/////////
	
#pragma omp parallel for
	for (int i = 0; i < object.groupNum; ++i) {
		object.groups[i].calMassMatrix(density);
		object.groups[i].calDampingMatrix();
		object.groups[i].calCenterofMass();
		object.groups[i].calInitCOM();//initial com
		object.groups[i].calLocalPos(); // initial local positions
		object.groups[i].calGroupK(youngs, poisson);		
		object.groups[i].setVertexMassesFromMassMatrix();//vertex mass
		object.groups[i].calMassGroup();
		object.groups[i].calMassDistributionMatrix();
		//object.groups[i].inverseTerm = (object.groups[i].massMatrix + object.groups[i].dampingMatrix * 0.01f).inverse(); //ﾋｳﾋｳｱ羃ﾑﾕ篋翆ﾋ
		//object.groups[i].inverseTermSparse = object.groups[i].inverseTerm.sparseView();
		object.groups[i].calLHS();
	}

	//for calculate frame rate
	double lastTime = glfwGetTime();
	int nbFrames = 0;
	glfwSwapInterval(0); //ｴｹﾖｱﾍｬｲｽ

	int frame = 1;
	while (!glfwWindowShouldClose(window)) {

		//object.commonPoints1 = object.findCommonVertices(object.groups[1], object.groups[2]);
		//ｹﾌｶｨｵ翹靹ﾃ


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
		//std::cout << wKey << std::endl;// ｵｱ W ｱｻｰｴﾏﾂﾊｱｵﾄﾂﾟｼｭ
		//double aa = object.groups[0].tetrahedra[0]->calMassTetra(density);
#pragma omp parallel for
		for (int i = 0; i < groupNum; i++) {
			//object.groups[i].calGroupKFEM(youngs, poisson);
			object.groups[i].calPrimeVec();

			//object.groups[i].calPrimeVec2(wKey);
			//object.groups[i].calPrimeVec(wKey);
			//object.groups[i].calPrimeVecS(wKey);
			//object.groups[i].calPrimeVecT(wKey);
			/*object.groups[i].calLHSFEM();
			object.groups[i].calRHSFEM();
			object.groups[i].calDeltaXFEM();
			object.groups[i].calculateCurrentPositionsFEM();
			object.groups[i].updateVelocityFEM();
			object.groups[i].updatePositionFEM();*/

			object.groups[i].calRotationMatrix();

		}
		/*for (int i = 0; i < groupNum; i++) {
			std::cout << "Group" << i << "Prime vector is" << std::endl << object.groups[i].primeVec;
		}*/


		object.PBDLOOP(8);

		if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
			for (int i = 0; i < groupNum; i++)
				std::cout << "Group Velocity-" << i << ":\n" << object.groups[i].groupVelocity << std::endl;

			system("pause");
		}


		// Render here
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//drawAxis1(0.3f, object.groups[0].rotate_matrix);
		// ｻ贍ﾆﾗ・・
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
		// ﾉ靹ﾃｵ羞ﾄｴ｡ﾎｪ10ﾏﾘ
		glPointSize(5.0f);
		// ｻ贍ﾆｰﾗﾉｫｵ・
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


		for (int groupIdx = 0; groupIdx < groupNum; ++groupIdx) { //ﾐｴｵ羞ﾄｱ・ﾅ｣ｬｻｭﾗ?
			Group& group = object.getGroup(groupIdx);

			////ｻｭｲｻﾖﾘｸｴｵﾄｰ豎ｾ
			//std::vector<Vertex*> uniqueVertices = group.getUniqueVertices();
			//for (size_t i = 0; i < uniqueVertices.size(); ++i) {
			//	Vertex* vertex = uniqueVertices[i];
			//	char buffer[5]; // ｷﾖﾅ葫羯ｻｴﾄｻｺｳ衂・
			//	sprintf_s(buffer, "%d", vertex->index); // ｽｫintﾗｪｻｻﾎｪchar*
			//	glColor3f(1, 0.0f, 0.0f);
			//	glRasterPos3f(vertex->x, vertex->y, vertex->z);
			//	XPrintString(buffer);
			//}


			//ｻｭﾖﾘｸｴｵﾄｰ豎ｾ


			//for (Tetrahedron* tetra : group.tetrahedra) { // ｱ鯊晙鰒ﾐｵﾄﾃｿｸﾄﾃ賣・
			//	for (int i = 0; i < 4; ++i) { // ｱ鯊撝ﾄﾃ賣蠏ﾄﾃｿｸ･ｵ・
			//		Vertex* vertex = tetra->vertices[i];
			//		char buffer[5]; // ｷﾖﾅ葫羯ｻｴﾄｻｺｳ衂・
			//		sprintf_s(buffer, "%d", vertex->index); // ｽｫintﾗｪｻｻﾎｪchar*
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

			//		//std::default_random_engine generator(vertex->index);//ﾋ貊摠ｷ｢ﾉ憘ｬﾓﾃﾓﾚﾗﾖｷ鈼ｫﾒﾆｷﾀﾖﾘｵ
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

		//	// ﾊｹﾓﾃﾗ鰒ﾐﾒﾑｾｭｼﾆﾋ羲ﾃｵﾄﾖﾊﾐﾄ
		//	Eigen::Vector3f& center = group.centerofMass;

		//	// ｽｫﾗ鰺ﾗｪｻｻﾎｪﾗﾖｷ逸ｮ
		//	char groupNumber[10];
		//	sprintf_s(groupNumber, "%d", groupIdx);

		//	// ﾎｪﾗ魍犲ﾅﾉ靹ﾃﾑﾕﾉｫ
		//	glColor3f(1.0f, 1.0f, 0.0f); // ｰﾗﾉｫﾓﾃﾓﾚﾎﾄｱｾ

		//	// ﾉ靹ﾃﾗ魍犲ﾅｵﾄﾎｻﾖﾃｲ｢ｻ贍ﾆ
		//	glRasterPos3f(center[0] + 1, center[1], center[2]);
		//	XPrintString(groupNumber);

		//	// ... [ｻ贍ﾆｱﾟﾔｵｺﾍｶ･ｵ羞ﾄﾆ萼犇惲・ ...
		//}

		// Draw edges
		glBegin(GL_LINES);




		//ｷﾖﾗ鰕ﾔﾊｾ
		for (int groupIdx = 0; groupIdx < groupNum; ++groupIdx) {
			float hhh;
			Group& group = object.getGroup(groupIdx);
			for (Tetrahedron* tet : group.tetrahedra) {
				for (int edgeIdx = 0; edgeIdx < 6; ++edgeIdx) {  // Loop through each edge in the tetrahedron
					Edge* edge = tet->edges[edgeIdx];
					Vertex* vertex1 = edge->vertices[0];
					Vertex* vertex2 = edge->vertices[1];
					bool isSurfaceEdge = edge->isBoundary;

					//Use HSV to RGB conversion to create a unique color for each group
					float hue = (360.0f * groupIdx) / groupNum;  // Distribute hues evenly across the spectrum
					hhh = hue;
					float saturation = 1.0f;  // Full saturation
					float value = 1.0f;      // Full brightness

					//Convert HSV to RGB
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
		//saveOBJ("43224.obj", object.groups);

		glPopMatrix();

		// Swap front and back buffers
		glfwSwapBuffers(window);

		// Poll for and process events
		glfwPollEvents();

		//calculate frame rate
		double currentTime = glfwGetTime();
		nbFrames++;
		if (currentTime - lastTime >= 1.0) { // ﾈ郢鄧ﾈ･ﾁﾋﾖﾁﾉﾙ1ﾃ・
			printf("%d frames/sec\n", nbFrames);
			nbFrames = 0;
			lastTime += 1.0;
		}
		printf("%d frame number\n", frame);
		frame++;

		

	}
	
	
	glfwTerminate();
	return 0;
}
