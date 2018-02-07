#include "Matrix10F10.hpp"
#include <math.h>

using namespace math;


const Matrix10F10::RotLookup_t Matrix10F10::RotLookup[] = {
	{  0,   0,   0 },
	{  0,   0,  45 },
	{  0,   0,  90 },
	{  0,   0, 135 },
	{  0,   0, 180 },
	{  0,   0, 225 },
	{  0,   0, 270 },
	{  0,   0, 315 },
	{180,   0,   0 },
	{180,   0,  45 },
	{180,   0,  90 },
	{180,   0, 135 },
	{  0, 180,   0 },
	{180,   0, 225 },
	{180,   0, 270 },
	{180,   0, 315 },
	{ 90,   0,   0 },
	{ 90,   0,  45 },
	{ 90,   0,  90 },
	{ 90,   0, 135 },
	{270,   0,   0 },
	{270,   0,  45 },
	{270,   0,  90 },
	{270,   0, 135 },
	{  0,  90,   0 },
	{  0, 270,   0 },
	{270,   0, 270 },
	{180, 270,   0 },
	{  0,  90, 180 },
	{ 90,  90,   0 },
	{ 90,  68, 293 },
	{270,  90,   0 },
	{  0,   9, 180 },
};

Matrix10F10::Matrix10F10(Vector10F m0, Vector10F m1, Vector10F m2, Vector10F m3, Vector10F m4, Vector10F m5, Vector10F m6, Vector10F m7, Vector10F m8, Vector10F m9) :
	data{m0, m1, m2, m3, m4, m5, m6, m7, m8, m9},
	nan{NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN}
{
};


Matrix10F10::Matrix10F10() :
	data{
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	},
    nan{NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN}
{
};


Matrix10F10::~Matrix10F10()
{
};


Vector10F& Matrix10F10::operator [] (uint32 i)
{
	if(i >= 10)
	{
		return nan;
	}
	else
	{
		return data[i];
	}
};


Vector10F Matrix10F10::operator [] (uint32 i) const
{
	if(i >= 10)
	{
		return nan;
	}
	else
	{
		return data[i];
	}
};


Matrix10F10 Matrix10F10::Transpose(void)
{
	Matrix10F10 res;

	res[0][0] = data[0][0];
	res[0][1] = data[1][0];
	res[0][2] = data[2][0];
	res[1][0] = data[0][1];
	res[1][1] = data[1][1];
	res[1][2] = data[2][1];
	res[2][0] = data[0][2];
	res[2][1] = data[1][2];
	res[2][2] = data[2][2];
	res[3][0] = data[0][3];
	res[3][1] = data[1][3];
	res[3][2] = data[2][3];
	res[4][0] = data[0][4];
	res[4][1] = data[1][4];
	res[4][2] = data[2][4];
	res[5][0] = data[0][5];
	res[5][1] = data[1][5];
	res[5][2] = data[2][5];
	res[6][0] = data[0][6];
	res[6][1] = data[1][6];
	res[6][2] = data[2][6];
	res[7][0] = data[0][7];
	res[7][1] = data[1][7];
	res[7][2] = data[2][7];
	res[8][0] = data[0][8];
	res[8][1] = data[1][8];
	res[8][2] = data[2][8];
	res[9][0] = data[0][9];
	res[9][1] = data[1][9];
	res[9][2] = data[2][9];

	return res;
}


Matrix10F10 Matrix10F10::Identity() {
    Matrix10F10 matOut(
    		{1.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    		{0.0, 1.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    		{0.0, 0.0, 1.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    		{0.0, 0.0, 0.0, 1.0f, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    		{0.0, 0.0, 0.0, 0.0, 1.0f, 0.0, 0.0, 0.0, 0.0, 0.0},
    		{0.0, 0.0, 0.0, 0.0, 0.0, 1.0f, 0.0, 0.0, 0.0, 0.0},
    		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0f, 0.0, 0.0, 0.0},
    		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0f, 0.0, 0.0},
    		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0f, 0.0},
    		{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0f});

    return matOut;
}


void Matrix10F10::Zero(void)
{
	data[0][0] = 0.0f;
	data[0][1] = 0.0f;
	data[0][2] = 0.0f;
	data[1][0] = 0.0f;
	data[1][1] = 0.0f;
	data[1][2] = 0.0f;
	data[2][0] = 0.0f;
	data[2][1] = 0.0f;
	data[2][2] = 0.0f;
	data[3][0] = 0.0f;
	data[3][1] = 0.0f;
	data[3][2] = 0.0f;
	data[4][0] = 0.0f;
	data[4][1] = 0.0f;
	data[4][2] = 0.0f;
	data[5][0] = 0.0f;
	data[5][1] = 0.0f;
	data[5][2] = 0.0f;
	data[6][0] = 0.0f;
	data[6][1] = 0.0f;
	data[6][2] = 0.0f;
	data[7][0] = 0.0f;
	data[7][1] = 0.0f;
	data[7][2] = 0.0f;
	data[8][0] = 0.0f;
	data[8][1] = 0.0f;
	data[8][2] = 0.0f;
	data[9][0] = 0.0f;
	data[9][1] = 0.0f;
	data[9][2] = 0.0f;
}


// overload * operator to provide a matrix product
Matrix10F10 Matrix10F10::operator*(const Matrix10F10 &matIn)
{
    Matrix10F10 matOut;

    matOut[0][0] = data[0][0]*matIn[0][0] + data[0][1]*matIn[1][0] + data[0][2]*matIn[2][0];
    matOut[0][1] = data[0][0]*matIn[0][1] + data[0][1]*matIn[1][1] + data[0][2]*matIn[2][1];
    matOut[0][2] = data[0][0]*matIn[0][2] + data[0][1]*matIn[1][2] + data[0][2]*matIn[2][2];

    matOut[1][0] = data[1][0]*matIn[0][0] + data[1][1]*matIn[1][0] + data[1][2]*matIn[2][0];
    matOut[1][1] = data[1][0]*matIn[0][1] + data[1][1]*matIn[1][1] + data[1][2]*matIn[2][1];
    matOut[1][2] = data[1][0]*matIn[0][2] + data[1][1]*matIn[1][2] + data[1][2]*matIn[2][2];

    matOut[2][0] = data[2][0]*matIn[0][0] + data[2][1]*matIn[1][0] + data[2][2]*matIn[2][0];
    matOut[2][1] = data[2][0]*matIn[0][1] + data[2][1]*matIn[1][1] + data[2][2]*matIn[2][1];
    matOut[2][2] = data[2][0]*matIn[0][2] + data[2][1]*matIn[1][2] + data[2][2]*matIn[2][2];

    matOut[3][0] = data[3][0]*matIn[0][0] + data[3][1]*matIn[1][0] + data[3][2]*matIn[2][0];
    matOut[3][1] = data[3][0]*matIn[0][1] + data[3][1]*matIn[1][1] + data[3][2]*matIn[2][1];
    matOut[3][2] = data[3][0]*matIn[0][2] + data[3][1]*matIn[1][2] + data[3][2]*matIn[2][2];

    matOut[4][0] = data[4][0]*matIn[0][0] + data[4][1]*matIn[1][0] + data[4][2]*matIn[2][0];
    matOut[4][1] = data[4][0]*matIn[0][1] + data[4][1]*matIn[1][1] + data[4][2]*matIn[2][1];
    matOut[4][2] = data[4][0]*matIn[0][2] + data[4][1]*matIn[1][2] + data[4][2]*matIn[2][2];

    matOut[5][0] = data[5][0]*matIn[0][0] + data[5][1]*matIn[1][0] + data[5][2]*matIn[2][0];
    matOut[5][1] = data[5][0]*matIn[0][1] + data[5][1]*matIn[1][1] + data[5][2]*matIn[2][1];
    matOut[5][2] = data[5][0]*matIn[0][2] + data[5][1]*matIn[1][2] + data[5][2]*matIn[2][2];

    matOut[6][0] = data[6][0]*matIn[0][0] + data[6][1]*matIn[1][0] + data[6][2]*matIn[2][0];
    matOut[6][1] = data[6][0]*matIn[0][1] + data[6][1]*matIn[1][1] + data[6][2]*matIn[2][1];
    matOut[6][2] = data[6][0]*matIn[0][2] + data[6][1]*matIn[1][2] + data[6][2]*matIn[2][2];

    matOut[7][0] = data[7][0]*matIn[0][0] + data[7][1]*matIn[1][0] + data[7][2]*matIn[2][0];
    matOut[7][1] = data[7][0]*matIn[0][1] + data[7][1]*matIn[1][1] + data[7][2]*matIn[2][1];
    matOut[7][2] = data[7][0]*matIn[0][2] + data[7][1]*matIn[1][2] + data[7][2]*matIn[2][2];

    matOut[8][0] = data[8][0]*matIn[0][0] + data[8][1]*matIn[1][0] + data[8][2]*matIn[2][0];
    matOut[8][1] = data[8][0]*matIn[0][1] + data[8][1]*matIn[1][1] + data[8][2]*matIn[2][1];
    matOut[8][2] = data[8][0]*matIn[0][2] + data[8][1]*matIn[1][2] + data[8][2]*matIn[2][2];

    matOut[9][0] = data[9][0]*matIn[0][0] + data[9][1]*matIn[1][0] + data[9][2]*matIn[2][0];
    matOut[9][1] = data[9][0]*matIn[0][1] + data[9][1]*matIn[1][1] + data[9][2]*matIn[2][1];
    matOut[9][2] = data[9][0]*matIn[0][2] + data[9][1]*matIn[1][2] + data[9][2]*matIn[2][2];

    return matOut;
}


// overload * operator to provide a matrix vector product
Vector10F Matrix10F10::operator*(const Vector10F &vecIn)
{
    Vector10F vecOut;
    vecOut[0] = data[0][0]*vecIn[0] + data[0][1]*vecIn[1] + data[0][2]*vecIn[2];
    vecOut[1] = data[1][0]*vecIn[0] + data[1][1]*vecIn[1] + data[1][2]*vecIn[2];
    vecOut[2] = data[2][0]*vecIn[0] + data[2][1]*vecIn[1] + data[2][2]*vecIn[2];
    vecOut[3] = data[3][0]*vecIn[0] + data[3][1]*vecIn[1] + data[3][2]*vecIn[2];
    vecOut[4] = data[4][0]*vecIn[0] + data[4][1]*vecIn[1] + data[4][2]*vecIn[2];
    vecOut[5] = data[5][0]*vecIn[0] + data[5][1]*vecIn[1] + data[5][2]*vecIn[2];
    vecOut[6] = data[6][0]*vecIn[0] + data[6][1]*vecIn[1] + data[6][2]*vecIn[2];
    vecOut[7] = data[7][0]*vecIn[0] + data[7][1]*vecIn[1] + data[7][2]*vecIn[2];
    vecOut[8] = data[8][0]*vecIn[0] + data[8][1]*vecIn[1] + data[8][2]*vecIn[2];
    vecOut[9] = data[9][0]*vecIn[0] + data[9][1]*vecIn[1] + data[9][2]*vecIn[2];

    return vecOut;
}


Matrix10F10 Matrix10F10::operator*(const float &scalar)
{
    Matrix10F10 matOut;

    matOut[0][0] = scalar*data[0][0];
    matOut[1][0] = scalar*data[1][0];
    matOut[2][0] = scalar*data[2][0];
    matOut[3][0] = scalar*data[3][0];
    matOut[4][0] = scalar*data[4][0];
    matOut[5][0] = scalar*data[5][0];
    matOut[6][0] = scalar*data[6][0];
    matOut[7][0] = scalar*data[7][0];
    matOut[8][0] = scalar*data[8][0];
    matOut[9][0] = scalar*data[9][0];

    matOut[0][1] = scalar*data[0][1];
    matOut[1][1] = scalar*data[1][1];
    matOut[2][1] = scalar*data[2][1];
    matOut[3][1] = scalar*data[3][1];
    matOut[4][1] = scalar*data[4][1];
    matOut[5][1] = scalar*data[5][1];
    matOut[6][1] = scalar*data[6][1];
    matOut[7][1] = scalar*data[7][1];
    matOut[8][1] = scalar*data[8][1];
    matOut[9][1] = scalar*data[9][1];

    matOut[0][2] = scalar*data[0][2];
    matOut[1][2] = scalar*data[1][2];
    matOut[2][2] = scalar*data[2][2];
    matOut[3][2] = scalar*data[3][2];
    matOut[4][2] = scalar*data[4][2];
    matOut[5][2] = scalar*data[5][2];
    matOut[6][2] = scalar*data[6][2];
    matOut[7][2] = scalar*data[7][2];
    matOut[8][2] = scalar*data[8][2];
    matOut[9][2] = scalar*data[9][2];

    return matOut;
}


Matrix10F10 Matrix10F10::operator+(const Matrix10F10 &matIn) const
{
    Matrix10F10 matOut;

    matOut[0][0] = data[0][0] + matIn[0][0];
    matOut[1][0] = data[1][0] + matIn[1][0];
    matOut[2][0] = data[2][0] + matIn[2][0];
    matOut[3][0] = data[3][0] + matIn[3][0];
    matOut[4][0] = data[4][0] + matIn[4][0];
    matOut[5][0] = data[5][0] + matIn[5][0];
    matOut[6][0] = data[6][0] + matIn[6][0];
    matOut[7][0] = data[7][0] + matIn[7][0];
    matOut[8][0] = data[8][0] + matIn[8][0];
    matOut[9][0] = data[9][0] + matIn[9][0];

    matOut[0][1] = data[0][1] + matIn[0][1];
    matOut[1][1] = data[1][1] + matIn[1][1];
    matOut[2][1] = data[2][1] + matIn[2][1];
    matOut[3][1] = data[3][1] + matIn[3][1];
    matOut[4][1] = data[4][1] + matIn[4][1];
    matOut[5][1] = data[5][1] + matIn[5][1];
    matOut[6][1] = data[6][1] + matIn[6][1];
    matOut[7][1] = data[7][1] + matIn[7][1];
    matOut[8][1] = data[8][1] + matIn[8][1];
    matOut[9][1] = data[9][1] + matIn[9][1];

    matOut[0][2] = data[0][2] + matIn[0][2];
    matOut[1][2] = data[1][2] + matIn[1][2];
    matOut[2][2] = data[2][2] + matIn[2][2];
    matOut[3][2] = data[3][2] + matIn[3][2];
    matOut[4][2] = data[4][2] + matIn[4][2];
    matOut[5][2] = data[5][2] + matIn[5][2];
    matOut[6][2] = data[6][2] + matIn[6][2];
    matOut[7][2] = data[7][2] + matIn[7][2];
    matOut[8][2] = data[8][2] + matIn[8][2];
    matOut[9][2] = data[9][2] + matIn[9][2];

    return matOut;
}


//Matrix10F10 Matrix10F10::RotationMatrix(Matrix10F10::Rotation_t boardRotation)
//{
//	Matrix10F10 matrix;
//
//	float roll  = M_DEG_TO_RAD_F * (float)Matrix10F10::RotLookup[boardRotation].roll;
//	float pitch = M_DEG_TO_RAD_F * (float)Matrix10F10::RotLookup[boardRotation].pitch;
//	float yaw   = M_DEG_TO_RAD_F * (float)Matrix10F10::RotLookup[boardRotation].yaw;
//
//	return Matrix10F10::FromEuler(roll, pitch, yaw);
//}
//
//
//
//Matrix10F10 Matrix10F10::FromEuler(float roll, float pitch, float yaw)
//{
//	Matrix10F10 matrix;
//
//	float cp = cosf(pitch);
//	float sp = sinf(pitch);
//	float cr = cosf(roll);
//	float sr = sinf(roll);
//	float cy = cosf(yaw);
//	float sy = sinf(yaw);
//
//	matrix[0][0] = cp * cy;
//	matrix[0][1] = (sr * sp * cy) - (cr * sy);
//	matrix[0][2] = (cr * sp * cy) + (sr * sy);
//	matrix[1][0] = cp * sy;
//	matrix[1][1] = (sr * sp * sy) + (cr * cy);
//	matrix[1][2] = (cr * sp * sy) - (sr * cy);
//	matrix[2][0] = -sp;
//	matrix[2][1] = sr * cp;
//	matrix[2][2] = cr * cp;
//
//	return matrix;
//}
//
//
//
//Vector10F Matrix10F10::ToEuler(void) const
//{
//	Vector10F euler;
//	euler[1] = asinf(-data[2][0]);
//
//	if (fabsf(euler[1] - M_PI_2) < 1.0e-3f) {
//		euler[0] = 0.0f;
//		euler[2] = atan2f(data[1][2] - data[0][1], data[0][2] + data[1][1]) + euler[0];
//
//	} else if (fabsf(euler[1] + M_PI_2) < 1.0e-3f) {
//		euler[0] = 0.0f;
//		euler[2] = atan2f(data[1][2] - data[0][1], data[0][2] + data[1][1]) - euler[0];
//
//	} else {
//		euler[0] = atan2f(data[2][1], data[2][2]);
//		euler[2] = atan2f(data[1][0], data[0][0]);
//	}
//
//	return euler;
//}
