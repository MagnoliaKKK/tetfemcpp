#include "Tetrahedron.h"

Eigen::MatrixXf Tetrahedron::createElementKAni(float E1, float E2, float E3, float nu, const Eigen::Vector3f& groupCenterOfMass)
{
	float x1 = vertices[0]->x - groupCenterOfMass.x();
	float y1 = vertices[0]->y - groupCenterOfMass.y();
	float z1 = vertices[0]->z - groupCenterOfMass.z();
	float x2 = vertices[1]->x - groupCenterOfMass.x();
	float y2 = vertices[1]->y - groupCenterOfMass.y();
	float z2 = vertices[1]->z - groupCenterOfMass.z();
	float x3 = vertices[2]->x - groupCenterOfMass.x();
	float y3 = vertices[2]->y - groupCenterOfMass.y();
	float z3 = vertices[2]->z - groupCenterOfMass.z();
	float x4 = vertices[3]->x - groupCenterOfMass.x();
	float y4 = vertices[3]->y - groupCenterOfMass.y();
	float z4 = vertices[3]->z - groupCenterOfMass.z();

	// ?Œš‘Ì??ŽZ‹é? A
	Eigen::Matrix4d A;
	A << x1, y1, z1, 1,
		x2, y2, z2, 1,
		x3, y3, z3, 1,
		x4, y4, z4, 1;

	// ?ŽZŽl–Ê‘Ì“I‘Ì?
	float V = std::abs(A.determinant() / 6);

	// ’è? mbeta, mgamma, mdelta ‹é?
	Eigen::Matrix3f mbeta1, mbeta2, mbeta3, mbeta4, mgamma1, mgamma2, mgamma3, mgamma4, mdelta1, mdelta2, mdelta3, mdelta4;


	mbeta1 << 1, y2, z2, 1, y3, z3, 1, y4, z4;
	mbeta2 << 1, y1, z1, 1, y3, z3, 1, y4, z4;
	mbeta3 << 1, y1, z1, 1, y2, z2, 1, y4, z4;
	mbeta4 << 1, y1, z1, 1, y2, z2, 1, y3, z3;

	mgamma1 << 1, x2, z2, 1, x3, z3, 1, x4, z4;
	mgamma2 << 1, x1, z1, 1, x3, z3, 1, x4, z4;
	mgamma3 << 1, x1, z1, 1, x2, z2, 1, x4, z4;
	mgamma4 << 1, x1, z1, 1, x2, z2, 1, x3, z3;

	mdelta1 << 1, x2, y2, 1, x3, y3, 1, x4, y4;
	mdelta2 << 1, x1, y1, 1, x3, y3, 1, x4, y4;
	mdelta3 << 1, x1, y1, 1, x2, y2, 1, x4, y4;
	mdelta4 << 1, x1, y1, 1, x2, y2, 1, x3, y3;

	// ?ŽZ beta, gamma ˜a delta ?
	float beta1 = -mbeta1.determinant();
	float beta2 = mbeta2.determinant();
	float beta3 = -mbeta3.determinant();
	float beta4 = mbeta4.determinant();

	float gamma1 = mgamma1.determinant();
	float gamma2 = -mgamma2.determinant();
	float gamma3 = mgamma3.determinant();
	float gamma4 = -mgamma4.determinant();

	float delta1 = -mdelta1.determinant();
	float delta2 = mdelta2.determinant();
	float delta3 = -mdelta3.determinant();
	float delta4 = mdelta4.determinant();

	// ’è? B ‹é?
	Eigen::MatrixXf B(6, 12);

	B << beta1, 0, 0, beta2, 0, 0, beta3, 0, 0, beta4, 0, 0,
		0, gamma1, 0, 0, gamma2, 0, 0, gamma3, 0, 0, gamma4, 0,
		0, 0, delta1, 0, 0, delta2, 0, 0, delta3, 0, 0, delta4,
		gamma1, beta1, 0, gamma2, beta2, 0, gamma3, beta3, 0, gamma4, beta4, 0,
		0, delta1, gamma1, 0, delta2, gamma2, 0, delta3, gamma3, 0, delta4, gamma4,
		delta1, 0, beta1, delta2, 0, beta2, delta3, 0, beta3, delta4, 0, beta4;

	B /= (6 * V);

	// ’è?Þ—¿‘®«‹é? D
	// ”‘¼”ä
	Eigen::MatrixXf D = Eigen::MatrixXf::Zero(6, 6);

	Eigen::Matrix3f A1;
	A1 << E1 * (1 - nu), std::sqrt(E1 * E2)* nu, std::sqrt(E1 * E3)* nu,
		std::sqrt(E1 * E2)* nu, E2* (1 - nu), std::sqrt(E2 * E3)* nu,
		std::sqrt(E1 * E3)* nu, std::sqrt(E2 * E3)* nu, E3* (1 - nu);

	A1 /= (nu + 1) * (1 - 2 * nu);

	// Calculate the orthotropic matrix 'B' based on the uploaded image
	Eigen::Matrix3f B1;
	B1 << std::sqrt(E1 * E2) / 2, 0, 0,
		0, std::sqrt(E2 * E3) / 2, 0,
		0, 0, std::sqrt(E1 * E3) / 2;

	B1 /= (1 + nu);

	// Now construct the complete compliance matrix
	//Eigen::MatrixXf D = Eigen::MatrixXf::Zero(6, 6);
	D.topLeftCorner<3, 3>() = A1;
	D.bottomRightCorner<3, 3>() = B1;

	// ?ŽZ?“x‹é? k
	Eigen::MatrixXf k = V * (B.transpose() * D * B);

	elementK = k;
	return k;
}
Eigen::MatrixXf Tetrahedron::createElementK(float E, float nu, const Eigen::Vector3f& groupCenterOfMass) {
	// ’è??“_¿?
	float x1 = vertices[0]->x - groupCenterOfMass.x();
	float y1 = vertices[0]->y - groupCenterOfMass.y();
	float z1 = vertices[0]->z - groupCenterOfMass.z();
	float x2 = vertices[1]->x - groupCenterOfMass.x();
	float y2 = vertices[1]->y - groupCenterOfMass.y();
	float z2 = vertices[1]->z - groupCenterOfMass.z();
	float x3 = vertices[2]->x - groupCenterOfMass.x();
	float y3 = vertices[2]->y - groupCenterOfMass.y();
	float z3 = vertices[2]->z - groupCenterOfMass.z();
	float x4 = vertices[3]->x - groupCenterOfMass.x();
	float y4 = vertices[3]->y - groupCenterOfMass.y();
	float z4 = vertices[3]->z - groupCenterOfMass.z();

	// ?Œš‘Ì??ŽZ‹é? A
	Eigen::Matrix4d A;
	A << x1, y1, z1, 1,
		x2, y2, z2, 1,
		x3, y3, z3, 1,
		x4, y4, z4, 1;

	// ?ŽZŽl–Ê‘Ì“I‘Ì?
	float V = std::abs(A.determinant() / 6);

	// ’è? mbeta, mgamma, mdelta ‹é?
	Eigen::Matrix3f mbeta1, mbeta2, mbeta3, mbeta4, mgamma1, mgamma2, mgamma3, mgamma4, mdelta1, mdelta2, mdelta3, mdelta4;


	mbeta1 << 1, y2, z2, 1, y3, z3, 1, y4, z4;
	mbeta2 << 1, y1, z1, 1, y3, z3, 1, y4, z4;
	mbeta3 << 1, y1, z1, 1, y2, z2, 1, y4, z4;
	mbeta4 << 1, y1, z1, 1, y2, z2, 1, y3, z3;

	mgamma1 << 1, x2, z2, 1, x3, z3, 1, x4, z4;
	mgamma2 << 1, x1, z1, 1, x3, z3, 1, x4, z4;
	mgamma3 << 1, x1, z1, 1, x2, z2, 1, x4, z4;
	mgamma4 << 1, x1, z1, 1, x2, z2, 1, x3, z3;

	mdelta1 << 1, x2, y2, 1, x3, y3, 1, x4, y4;
	mdelta2 << 1, x1, y1, 1, x3, y3, 1, x4, y4;
	mdelta3 << 1, x1, y1, 1, x2, y2, 1, x4, y4;
	mdelta4 << 1, x1, y1, 1, x2, y2, 1, x3, y3;

	// ?ŽZ beta, gamma ˜a delta ?
	float beta1 = -mbeta1.determinant();
	float beta2 = mbeta2.determinant();
	float beta3 = -mbeta3.determinant();
	float beta4 = mbeta4.determinant();

	float gamma1 = mgamma1.determinant();
	float gamma2 = -mgamma2.determinant();
	float gamma3 = mgamma3.determinant();
	float gamma4 = -mgamma4.determinant();

	float delta1 = -mdelta1.determinant();
	float delta2 = mdelta2.determinant();
	float delta3 = -mdelta3.determinant();
	float delta4 = mdelta4.determinant();

	// ’è? B ‹é?
	Eigen::MatrixXf B(6, 12);

	B << beta1, 0, 0, beta2, 0, 0, beta3, 0, 0, beta4, 0, 0,
		0, gamma1, 0, 0, gamma2, 0, 0, gamma3, 0, 0, gamma4, 0,
		0, 0, delta1, 0, 0, delta2, 0, 0, delta3, 0, 0, delta4,
		gamma1, beta1, 0, gamma2, beta2, 0, gamma3, beta3, 0, gamma4, beta4, 0,
		0, delta1, gamma1, 0, delta2, gamma2, 0, delta3, gamma3, 0, delta4, gamma4,
		delta1, 0, beta1, delta2, 0, beta2, delta3, 0, beta3, delta4, 0, beta4;

	B /= (6 * V);

	// ’è?Þ—¿‘®«‹é? D
	// ”‘¼”ä
	Eigen::MatrixXf D = Eigen::MatrixXf::Zero(6, 6);

	D << 1 - nu, nu, nu, 0, 0, 0,
		nu, 1 - nu, nu, 0, 0, 0,
		nu, nu, 1 - nu, 0, 0, 0,
		0, 0, 0, (1 - 2 * nu) / 2, 0, 0,
		0, 0, 0, 0, (1 - 2 * nu) / 2, 0,
		0, 0, 0, 0, 0, (1 - 2 * nu) / 2;

	D *= (E / ((1 + nu) * (1 - 2 * nu)));

	// ?ŽZ?“x‹é? k
	Eigen::MatrixXf k = V * (B.transpose() * D * B);

	elementK = k;
	return k;
}

Eigen::MatrixXf Tetrahedron::createElementKFEM(float E, float nu) {
	// ’è??“_¿?
	float x1 = vertices[0]->x;
	float y1 = vertices[0]->y;
	float z1 = vertices[0]->z;
	float x2 = vertices[1]->x;
	float y2 = vertices[1]->y;
	float z2 = vertices[1]->z;
	float x3 = vertices[2]->x;
	float y3 = vertices[2]->y;
	float z3 = vertices[2]->z;
	float x4 = vertices[3]->x;
	float y4 = vertices[3]->y;
	float z4 = vertices[3]->z;

	// ?Œš‘Ì??ŽZ‹é? A
	Eigen::Matrix4d A;
	A << x1, y1, z1, 1,
		x2, y2, z2, 1,
		x3, y3, z3, 1,
		x4, y4, z4, 1;

	// ?ŽZŽl–Ê‘Ì“I‘Ì?
	float V = std::abs(A.determinant() / 6);

	// ’è? mbeta, mgamma, mdelta ‹é?
	Eigen::Matrix3f mbeta1, mbeta2, mbeta3, mbeta4, mgamma1, mgamma2, mgamma3, mgamma4, mdelta1, mdelta2, mdelta3, mdelta4;


	mbeta1 << 1, y2, z2, 1, y3, z3, 1, y4, z4;
	mbeta2 << 1, y1, z1, 1, y3, z3, 1, y4, z4;
	mbeta3 << 1, y1, z1, 1, y2, z2, 1, y4, z4;
	mbeta4 << 1, y1, z1, 1, y2, z2, 1, y3, z3;

	mgamma1 << 1, x2, z2, 1, x3, z3, 1, x4, z4;
	mgamma2 << 1, x1, z1, 1, x3, z3, 1, x4, z4;
	mgamma3 << 1, x1, z1, 1, x2, z2, 1, x4, z4;
	mgamma4 << 1, x1, z1, 1, x2, z2, 1, x3, z3;

	mdelta1 << 1, x2, y2, 1, x3, y3, 1, x4, y4;
	mdelta2 << 1, x1, y1, 1, x3, y3, 1, x4, y4;
	mdelta3 << 1, x1, y1, 1, x2, y2, 1, x4, y4;
	mdelta4 << 1, x1, y1, 1, x2, y2, 1, x3, y3;

	// ?ŽZ beta, gamma ˜a delta ?
	float beta1 = -mbeta1.determinant();
	float beta2 = mbeta2.determinant();
	float beta3 = -mbeta3.determinant();
	float beta4 = mbeta4.determinant();

	float gamma1 = mgamma1.determinant();
	float gamma2 = -mgamma2.determinant();
	float gamma3 = mgamma3.determinant();
	float gamma4 = -mgamma4.determinant();

	float delta1 = -mdelta1.determinant();
	float delta2 = mdelta2.determinant();
	float delta3 = -mdelta3.determinant();
	float delta4 = mdelta4.determinant();

	// ’è? B ‹é?
	Eigen::MatrixXf B(6, 12);

	B << beta1, 0, 0, beta2, 0, 0, beta3, 0, 0, beta4, 0, 0,
		0, gamma1, 0, 0, gamma2, 0, 0, gamma3, 0, 0, gamma4, 0,
		0, 0, delta1, 0, 0, delta2, 0, 0, delta3, 0, 0, delta4,
		gamma1, beta1, 0, gamma2, beta2, 0, gamma3, beta3, 0, gamma4, beta4, 0,
		0, delta1, gamma1, 0, delta2, gamma2, 0, delta3, gamma3, 0, delta4, gamma4,
		delta1, 0, beta1, delta2, 0, beta2, delta3, 0, beta3, delta4, 0, beta4;

	B /= (6 * V);

	// ’è?Þ—¿‘®«‹é? D
	// ”‘¼”ä
	Eigen::MatrixXf D = Eigen::MatrixXf::Zero(6, 6);

	D << 1 - nu, nu, nu, 0, 0, 0,
		nu, 1 - nu, nu, 0, 0, 0,
		nu, nu, 1 - nu, 0, 0, 0,
		0, 0, 0, (1 - 2 * nu) / 2, 0, 0,
		0, 0, 0, 0, (1 - 2 * nu) / 2, 0,
		0, 0, 0, 0, 0, (1 - 2 * nu) / 2;

	D *= (E / ((1 + nu) * (1 - 2 * nu)));

	// ?ŽZ?“x‹é? k
	Eigen::MatrixXf k = V * (B.transpose() * D * B);

	elementKFEM = k;
	return k;
}
float Tetrahedron::calMassTetra(float den) {

	//float volume;
	Eigen::Vector3f AB(vertices[1]->x - vertices[0]->x, vertices[1]->y - vertices[0]->y, vertices[1]->z - vertices[0]->z);
	Eigen::Vector3f AC(vertices[2]->x - vertices[0]->x, vertices[2]->y - vertices[0]->y, vertices[2]->z - vertices[0]->z);
	Eigen::Vector3f AD(vertices[3]->x - vertices[0]->x, vertices[3]->y - vertices[0]->y, vertices[3]->z - vertices[0]->z);

	// Calculate volume using the formula
	volumeTetra = (AB.cross(AC)).dot(AD) / 6.0f;
	volumeTetra = std::abs(volumeTetra);
	massTetra = volumeTetra * den;
	return massTetra;


}
float Tetrahedron::calVolumeTetra() {

	//float volume;
	Eigen::Vector3f AB(vertices[1]->x - vertices[0]->x, vertices[1]->y - vertices[0]->y, vertices[1]->z - vertices[0]->z);
	Eigen::Vector3f AC(vertices[2]->x - vertices[0]->x, vertices[2]->y - vertices[0]->y, vertices[2]->z - vertices[0]->z);
	Eigen::Vector3f AD(vertices[3]->x - vertices[0]->x, vertices[3]->y - vertices[0]->y, vertices[3]->z - vertices[0]->z);

	// Calculate volume using the formula
	volumeTetra = (AB.cross(AC)).dot(AD) / 6.0f;
	volumeTetra = std::abs(volumeTetra);
	return volumeTetra;
}