#include <Novice.h>
#define _USE_MATH_DEFIENS
#include <cmath>
#include <assert.h>
#include <Novice.h>
#include "Matrix4x4.h"
#include "Vector2.h"
#include "Vector3.h"
#include <math.h>
#include <algorithm>

const char kWindowTitle[] = "MT4";
static const int kRowHeight = 20;
static const int kColumnWidth = 80;


struct Quaternion final {
	float x;
	float y;
	float z;
	float w;
};

Matrix4x4 Add(const Matrix4x4& m1, const Matrix4x4& m2) {

	Matrix4x4 result;

	for (int y = 0; y < 4; y++) {
		for (int x = 0; x < 4; x++) {
			result.m[y][x] = m1.m[y][x] + m2.m[y][x];
		}
	}
	return result;
}
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2)
{
	Matrix4x4 result;

	result.m[0][0] = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] + m1.m[0][2] * m2.m[2][0] + m1.m[0][3] * m2.m[3][0];
	result.m[0][1] = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] + m1.m[0][2] * m2.m[2][1] + m1.m[0][3] * m2.m[3][1];
	result.m[0][2] = m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] + m1.m[0][2] * m2.m[2][2] + m1.m[0][3] * m2.m[3][2];
	result.m[0][3] = m1.m[0][0] * m2.m[0][3] + m1.m[0][1] * m2.m[1][3] + m1.m[0][2] * m2.m[2][3] + m1.m[0][3] * m2.m[3][3];

	result.m[1][0] = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] + m1.m[1][2] * m2.m[2][0] + m1.m[1][3] * m2.m[3][0];
	result.m[1][1] = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] + m1.m[1][2] * m2.m[2][1] + m1.m[1][3] * m2.m[3][1];
	result.m[1][2] = m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] + m1.m[1][2] * m2.m[2][2] + m1.m[1][3] * m2.m[3][2];
	result.m[1][3] = m1.m[1][0] * m2.m[0][3] + m1.m[1][1] * m2.m[1][3] + m1.m[1][2] * m2.m[2][3] + m1.m[1][3] * m2.m[3][3];

	result.m[2][0] = m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] + m1.m[2][2] * m2.m[2][0] + m1.m[2][3] * m2.m[3][0];
	result.m[2][1] = m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] + m1.m[2][2] * m2.m[2][1] + m1.m[2][3] * m2.m[3][1];
	result.m[2][2] = m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] + m1.m[2][2] * m2.m[2][2] + m1.m[2][3] * m2.m[3][2];
	result.m[2][3] = m1.m[2][0] * m2.m[0][3] + m1.m[2][1] * m2.m[1][3] + m1.m[2][2] * m2.m[2][3] + m1.m[2][3] * m2.m[3][3];

	result.m[3][0] = m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] + m1.m[3][2] * m2.m[2][0] + m1.m[3][3] * m2.m[3][0];
	result.m[3][1] = m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] + m1.m[3][2] * m2.m[2][1] + m1.m[3][3] * m2.m[3][1];
	result.m[3][2] = m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] + m1.m[3][2] * m2.m[2][2] + m1.m[3][3] * m2.m[3][2];
	result.m[3][3] = m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] + m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3];

	return result;
}

Matrix4x4 Multiply(float scalar, const Matrix4x4& m)
{

	Matrix4x4 result;

	for (int y = 0; y < 4; y++) {
		for (int x = 0; x < 4; x++) {
			result.m[y][x] = scalar * m.m[y][x];
		}
	}
	return result;
}
//拡大縮小行列
Matrix4x4 MakeScaleMatrix(const Vector3& scale)
{
	Matrix4x4 result;

	result.m[0][0] = scale.x;
	result.m[0][1] = 0.0f;
	result.m[0][2] = 0.0f;
	result.m[0][3] = 0.0f;

	result.m[1][0] = 0.0f;
	result.m[1][1] = scale.y;
	result.m[1][2] = 0.0f;
	result.m[1][3] = 0.0f;

	result.m[2][0] = 0.0f;
	result.m[2][1] = 0.0f;
	result.m[2][2] = scale.z;
	result.m[2][3] = 0.0f;

	result.m[3][0] = 0.0f;
	result.m[3][1] = 0.0f;
	result.m[3][2] = 0.0f;
	result.m[3][3] = 1.0f;

	return result;
}

//  X軸回転
Matrix4x4 MakeRotateXMatrix(float radian) {
	Matrix4x4 result;

	result.m[0][0] = 1;
	result.m[0][1] = 0;
	result.m[0][2] = 0;
	result.m[0][3] = 0;

	result.m[1][0] = 0;
	result.m[1][1] = std::cos(radian);
	result.m[1][2] = std::sin(radian);
	result.m[1][3] = 0;

	result.m[2][0] = 0;
	result.m[2][1] = -std::sin(radian);
	result.m[2][2] = std::cos(radian);
	result.m[2][3] = 0;

	result.m[3][0] = 0;
	result.m[3][1] = 0;
	result.m[3][2] = 0;
	result.m[3][3] = 1;

	return result;
}
// Y軸回転
Matrix4x4 MakeRotateYMatrix(float radian) {
	Matrix4x4 result;

	result.m[0][0] = std::cos(radian);
	result.m[0][1] = 0;
	result.m[0][2] = -std::sin(radian);
	result.m[0][3] = 0;

	result.m[1][0] = 0;
	result.m[1][1] = 1;
	result.m[1][2] = 0;
	result.m[1][3] = 0;

	result.m[2][0] = std::sin(radian);
	result.m[2][1] = 0;
	result.m[2][2] = std::cos(radian);
	result.m[2][3] = 0;

	result.m[3][0] = 0;
	result.m[3][1] = 0;
	result.m[3][2] = 0;
	result.m[3][3] = 1;

	return result;
}
// Z軸回転
Matrix4x4 MakeRotateZMatrix(float radian) {
	Matrix4x4 result;

	result.m[0][0] = std::cos(radian);
	result.m[0][1] = std::sin(radian);
	result.m[0][2] = 0;
	result.m[0][3] = 0;

	result.m[1][0] = -std::sin(radian);
	result.m[1][1] = std::cos(radian);
	result.m[1][2] = 0;
	result.m[1][3] = 0;

	result.m[2][0] = 0;
	result.m[2][1] = 0;
	result.m[2][2] = 1;
	result.m[2][3] = 0;

	result.m[3][0] = 0;
	result.m[3][1] = 0;
	result.m[3][2] = 0;
	result.m[3][3] = 1;

	return result;
}
// 回転行列
Matrix4x4 MakeRotateMatrix(const Vector3& radian) {
	Matrix4x4 result;
	Matrix4x4 RotateX = MakeRotateXMatrix(radian.x);
	Matrix4x4 RotateY = MakeRotateYMatrix(radian.y);
	Matrix4x4 RotateZ = MakeRotateZMatrix(radian.z);

	result = Multiply(RotateX, Multiply(RotateY, RotateZ));

	return result;
}

Matrix4x4 MakeTranslateMatrix(const Vector3& translate)
{
	Matrix4x4 result;

	result.m[0][0] = 1.0f;
	result.m[0][1] = 0.0f;
	result.m[0][2] = 0.0f;
	result.m[0][3] = 0.0f;

	result.m[1][0] = 0.0f;
	result.m[1][1] = 1.0f;
	result.m[1][2] = 0.0f;
	result.m[1][3] = 0.0f;

	result.m[2][0] = 0.0f;
	result.m[2][1] = 0.0f;
	result.m[2][2] = 1.0f;
	result.m[2][3] = 0.0f;

	result.m[3][0] = translate.x;
	result.m[3][1] = translate.y;
	result.m[3][2] = translate.z;
	result.m[3][3] = 1.0f;

	return result;
}

//アフィン行列
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rot, const Vector3& translate)
{
	Matrix4x4 result;

	// スケーリング行列＊回転行列＊平行移動行列
	result = Multiply(MakeScaleMatrix(scale), Multiply(MakeRotateMatrix(rot), MakeTranslateMatrix(translate)));

	return result;
}

//単位行列の作成
Matrix4x4 MakeIdentity4x4() {

	Matrix4x4 result;

	for (int y = 0; y < 4; y++) {
		for (int x = 0; x < 4; x++) {
			if (y == x) {
				result.m[y][x] = 1.0f;
			}
			else {
				result.m[y][x] = 0.0f;
			}
		}
	}

	return result;
}

Matrix4x4 MakeRotateAxisAngle(const Vector3& axis, float angle)
{
	//S
	Matrix4x4 matrixS = MakeIdentity4x4();
	matrixS.m[0][0] = std::cosf(angle);
	matrixS.m[1][1] = std::cosf(angle);
	matrixS.m[2][2] = std::cosf(angle);

	//P
	Matrix4x4 matrixP = MakeIdentity4x4();
	matrixP.m[0][0] = axis.x * axis.x;
	matrixP.m[0][1] = axis.x * axis.y;
	matrixP.m[0][2] = axis.x * axis.z;
	matrixP.m[1][0] = axis.y * axis.x;
	matrixP.m[1][1] = axis.y * axis.y;
	matrixP.m[1][2] = axis.y * axis.z;
	matrixP.m[2][0] = axis.z * axis.x;
	matrixP.m[2][1] = axis.z * axis.y;
	matrixP.m[2][2] = axis.z * axis.z;
	matrixP.m[3][3] = 0.0f;
	matrixP = Multiply((1.0f - std::cosf(angle)), matrixP);

	//C
	Matrix4x4 matrixC = MakeIdentity4x4();
	matrixC.m[0][0] = 0.0f;
	matrixC.m[0][1] = -axis.z;
	matrixC.m[0][2] = axis.y;
	matrixC.m[1][0] = axis.z;
	matrixC.m[1][1] = 0.0f;
	matrixC.m[1][2] = -axis.x;
	matrixC.m[2][0] = -axis.y;
	matrixC.m[2][1] = axis.x;
	matrixC.m[2][2] = 0.0f;
	matrixC.m[3][3] = 0.0f;
	matrixC = Multiply(-std::sinf(angle), matrixC);

	// result
	Matrix4x4 resultMatrix = Add(Add(matrixS, matrixP), matrixC);

	return resultMatrix;
}

Vector3 Multiply(float scalar, const Vector3& v) {

	Vector3 result;

	result = { scalar * v.x, scalar * v.y, scalar * v.z };

	return result;

}

//内積
float Dot(const Vector3& v1, const Vector3& v2) {

	float result;

	result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

	return result;

}
Vector3 Cross(const Vector3& v1, const Vector3& v2) {

	Vector3 result = { v1.y * v2.z - v1.z * v2.y,v1.z * v2.x - v1.x * v2.z,v1.x * v2.y - v1.y * v2.x, };

	return result;

}

// 長さ(ノルム)
float Length(const Vector3& v) {
	float result;
	result = powf(v.x, 2.0) + powf(v.y, 2.0) + powf(v.z, 2.0);

	return sqrtf(result);
};


Matrix4x4 Add(const Matrix4x4& m1, const Matrix4x4& m2) {

	Matrix4x4 result;

	for (int y = 0; y < 4; y++) {
		for (int x = 0; x < 4; x++) {
			result.m[y][x] = m1.m[y][x] + m2.m[y][x];
		}
	}
	return result;
}
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2)
{
	Matrix4x4 result;

	result.m[0][0] = m1.m[0][0] * m2.m[0][0] + m1.m[0][1] * m2.m[1][0] + m1.m[0][2] * m2.m[2][0] + m1.m[0][3] * m2.m[3][0];
	result.m[0][1] = m1.m[0][0] * m2.m[0][1] + m1.m[0][1] * m2.m[1][1] + m1.m[0][2] * m2.m[2][1] + m1.m[0][3] * m2.m[3][1];
	result.m[0][2] = m1.m[0][0] * m2.m[0][2] + m1.m[0][1] * m2.m[1][2] + m1.m[0][2] * m2.m[2][2] + m1.m[0][3] * m2.m[3][2];
	result.m[0][3] = m1.m[0][0] * m2.m[0][3] + m1.m[0][1] * m2.m[1][3] + m1.m[0][2] * m2.m[2][3] + m1.m[0][3] * m2.m[3][3];

	result.m[1][0] = m1.m[1][0] * m2.m[0][0] + m1.m[1][1] * m2.m[1][0] + m1.m[1][2] * m2.m[2][0] + m1.m[1][3] * m2.m[3][0];
	result.m[1][1] = m1.m[1][0] * m2.m[0][1] + m1.m[1][1] * m2.m[1][1] + m1.m[1][2] * m2.m[2][1] + m1.m[1][3] * m2.m[3][1];
	result.m[1][2] = m1.m[1][0] * m2.m[0][2] + m1.m[1][1] * m2.m[1][2] + m1.m[1][2] * m2.m[2][2] + m1.m[1][3] * m2.m[3][2];
	result.m[1][3] = m1.m[1][0] * m2.m[0][3] + m1.m[1][1] * m2.m[1][3] + m1.m[1][2] * m2.m[2][3] + m1.m[1][3] * m2.m[3][3];

	result.m[2][0] = m1.m[2][0] * m2.m[0][0] + m1.m[2][1] * m2.m[1][0] + m1.m[2][2] * m2.m[2][0] + m1.m[2][3] * m2.m[3][0];
	result.m[2][1] = m1.m[2][0] * m2.m[0][1] + m1.m[2][1] * m2.m[1][1] + m1.m[2][2] * m2.m[2][1] + m1.m[2][3] * m2.m[3][1];
	result.m[2][2] = m1.m[2][0] * m2.m[0][2] + m1.m[2][1] * m2.m[1][2] + m1.m[2][2] * m2.m[2][2] + m1.m[2][3] * m2.m[3][2];
	result.m[2][3] = m1.m[2][0] * m2.m[0][3] + m1.m[2][1] * m2.m[1][3] + m1.m[2][2] * m2.m[2][3] + m1.m[2][3] * m2.m[3][3];

	result.m[3][0] = m1.m[3][0] * m2.m[0][0] + m1.m[3][1] * m2.m[1][0] + m1.m[3][2] * m2.m[2][0] + m1.m[3][3] * m2.m[3][0];
	result.m[3][1] = m1.m[3][0] * m2.m[0][1] + m1.m[3][1] * m2.m[1][1] + m1.m[3][2] * m2.m[2][1] + m1.m[3][3] * m2.m[3][1];
	result.m[3][2] = m1.m[3][0] * m2.m[0][2] + m1.m[3][1] * m2.m[1][2] + m1.m[3][2] * m2.m[2][2] + m1.m[3][3] * m2.m[3][2];
	result.m[3][3] = m1.m[3][0] * m2.m[0][3] + m1.m[3][1] * m2.m[1][3] + m1.m[3][2] * m2.m[2][3] + m1.m[3][3] * m2.m[3][3];

	return result;
}

Matrix4x4 Multiply(float scalar, const Matrix4x4& m)
{

	Matrix4x4 result;

	for (int y = 0; y < 4; y++) {
		for (int x = 0; x < 4; x++) {
			result.m[y][x] = scalar * m.m[y][x];
		}
	}
	return result;
}
//拡大縮小行列
Matrix4x4 MakeScaleMatrix(const Vector3& scale)
{
	Matrix4x4 result;

	result.m[0][0] = scale.x;
	result.m[0][1] = 0.0f;
	result.m[0][2] = 0.0f;
	result.m[0][3] = 0.0f;

	result.m[1][0] = 0.0f;
	result.m[1][1] = scale.y;
	result.m[1][2] = 0.0f;
	result.m[1][3] = 0.0f;

	result.m[2][0] = 0.0f;
	result.m[2][1] = 0.0f;
	result.m[2][2] = scale.z;
	result.m[2][3] = 0.0f;

	result.m[3][0] = 0.0f;
	result.m[3][1] = 0.0f;
	result.m[3][2] = 0.0f;
	result.m[3][3] = 1.0f;

	return result;
}

//  X軸回転
Matrix4x4 MakeRotateXMatrix(float radian) {
	Matrix4x4 result;

	result.m[0][0] = 1;
	result.m[0][1] = 0;
	result.m[0][2] = 0;
	result.m[0][3] = 0;

	result.m[1][0] = 0;
	result.m[1][1] = std::cos(radian);
	result.m[1][2] = std::sin(radian);
	result.m[1][3] = 0;

	result.m[2][0] = 0;
	result.m[2][1] = -std::sin(radian);
	result.m[2][2] = std::cos(radian);
	result.m[2][3] = 0;

	result.m[3][0] = 0;
	result.m[3][1] = 0;
	result.m[3][2] = 0;
	result.m[3][3] = 1;

	return result;
}
// Y軸回転
Matrix4x4 MakeRotateYMatrix(float radian) {
	Matrix4x4 result;

	result.m[0][0] = std::cos(radian);
	result.m[0][1] = 0;
	result.m[0][2] = -std::sin(radian);
	result.m[0][3] = 0;

	result.m[1][0] = 0;
	result.m[1][1] = 1;
	result.m[1][2] = 0;
	result.m[1][3] = 0;

	result.m[2][0] = std::sin(radian);
	result.m[2][1] = 0;
	result.m[2][2] = std::cos(radian);
	result.m[2][3] = 0;

	result.m[3][0] = 0;
	result.m[3][1] = 0;
	result.m[3][2] = 0;
	result.m[3][3] = 1;

	return result;
}
// Z軸回転
Matrix4x4 MakeRotateZMatrix(float radian) {
	Matrix4x4 result;

	result.m[0][0] = std::cos(radian);
	result.m[0][1] = std::sin(radian);
	result.m[0][2] = 0;
	result.m[0][3] = 0;

	result.m[1][0] = -std::sin(radian);
	result.m[1][1] = std::cos(radian);
	result.m[1][2] = 0;
	result.m[1][3] = 0;

	result.m[2][0] = 0;
	result.m[2][1] = 0;
	result.m[2][2] = 1;
	result.m[2][3] = 0;

	result.m[3][0] = 0;
	result.m[3][1] = 0;
	result.m[3][2] = 0;
	result.m[3][3] = 1;

	return result;
}
// 回転行列
Matrix4x4 MakeRotateMatrix(const Vector3& radian) {
	Matrix4x4 result;
	Matrix4x4 RotateX = MakeRotateXMatrix(radian.x);
	Matrix4x4 RotateY = MakeRotateYMatrix(radian.y);
	Matrix4x4 RotateZ = MakeRotateZMatrix(radian.z);

	result = Multiply(RotateX, Multiply(RotateY, RotateZ));

	return result;
}

Matrix4x4 MakeTranslateMatrix(const Vector3& translate)
{
	Matrix4x4 result;

	result.m[0][0] = 1.0f;
	result.m[0][1] = 0.0f;
	result.m[0][2] = 0.0f;
	result.m[0][3] = 0.0f;

	result.m[1][0] = 0.0f;
	result.m[1][1] = 1.0f;
	result.m[1][2] = 0.0f;
	result.m[1][3] = 0.0f;

	result.m[2][0] = 0.0f;
	result.m[2][1] = 0.0f;
	result.m[2][2] = 1.0f;
	result.m[2][3] = 0.0f;

	result.m[3][0] = translate.x;
	result.m[3][1] = translate.y;
	result.m[3][2] = translate.z;
	result.m[3][3] = 1.0f;

	return result;
}

//アフィン行列
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rot, const Vector3& translate)
{
	Matrix4x4 result;

	// スケーリング行列＊回転行列＊平行移動行列
	result = Multiply(MakeScaleMatrix(scale), Multiply(MakeRotateMatrix(rot), MakeTranslateMatrix(translate)));

	return result;
}

//単位行列の作成
Matrix4x4 MakeIdentity4x4() {

	Matrix4x4 result;

	for (int y = 0; y < 4; y++) {
		for (int x = 0; x < 4; x++) {
			if (y == x) {
				result.m[y][x] = 1.0f;
			}
			else {
				result.m[y][x] = 0.0f;
			}
		}
	}

	return result;
}

Matrix4x4 MakeRotateAxisAngle(const Vector3& axis, float angle)
{
	//S
	Matrix4x4 matrixS = MakeIdentity4x4();
	matrixS.m[0][0] = std::cosf(angle);
	matrixS.m[1][1] = std::cosf(angle);
	matrixS.m[2][2] = std::cosf(angle);

	//P
	Matrix4x4 matrixP = MakeIdentity4x4();
	matrixP.m[0][0] = axis.x * axis.x;
	matrixP.m[0][1] = axis.x * axis.y;
	matrixP.m[0][2] = axis.x * axis.z;
	matrixP.m[1][0] = axis.y * axis.x;
	matrixP.m[1][1] = axis.y * axis.y;
	matrixP.m[1][2] = axis.y * axis.z;
	matrixP.m[2][0] = axis.z * axis.x;
	matrixP.m[2][1] = axis.z * axis.y;
	matrixP.m[2][2] = axis.z * axis.z;
	matrixP.m[3][3] = 0.0f;
	matrixP = Multiply((1.0f - std::cosf(angle)), matrixP);

	//C
	Matrix4x4 matrixC = MakeIdentity4x4();
	matrixC.m[0][0] = 0.0f;
	matrixC.m[0][1] = -axis.z;
	matrixC.m[0][2] = axis.y;
	matrixC.m[1][0] = axis.z;
	matrixC.m[1][1] = 0.0f;
	matrixC.m[1][2] = -axis.x;
	matrixC.m[2][0] = -axis.y;
	matrixC.m[2][1] = axis.x;
	matrixC.m[2][2] = 0.0f;
	matrixC.m[3][3] = 0.0f;
	matrixC = Multiply(-std::sinf(angle), matrixC);

	// result
	Matrix4x4 resultMatrix = Add(Add(matrixS, matrixP), matrixC);

	return resultMatrix;
}

Vector3 Multiply(float scalar, const Vector3& v) {

	Vector3 result;

	result = { scalar * v.x, scalar * v.y, scalar * v.z };

	return result;

}

//内積
float Dot(const Vector3& v1, const Vector3& v2) {

	float result;

	result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

	return result;

}
Vector3 Cross(const Vector3& v1, const Vector3& v2) {

	Vector3 result = { v1.y * v2.z - v1.z * v2.y,v1.z * v2.x - v1.x * v2.z,v1.x * v2.y - v1.y * v2.x, };

	return result;

}

// 長さ(ノルム)
float Length(const Vector3& v) {
	float result;
	result = powf(v.x, 2.0) + powf(v.y, 2.0) + powf(v.z, 2.0);

	return sqrtf(result);
};

// 正規化
Vector3 Normalize(const Vector3& v) {
	Vector3 result;
	float x;
	x = Length(v);
	result.x = v.x / x;
	result.y = v.y / x;
	result.z = v.z / x;
	return result;
}

Matrix4x4 DirectionToDirection(const Vector3& from, const Vector3& to) {
	Vector3 normal = Normalize(Cross(from, to));
	
	Vector3 minusTo = Multiply(-1.0f, to);
	Matrix4x4 result = MakeIdentity4x4();
	if ((from.x == minusTo.x &&
		from.y == minusTo.y &&
		from.z == minusTo.z)) {
		if (from.x != 0.0f || from.y != 0.0f) {
			normal = { from.y, -from.x, 0.0f };
		}
		else if (from.x != 0.0f || from.z != 0.0f) {
			normal = { from.z, 0.0f, -from.x };
		}
	}
	float cos = Dot(from, to);
	float sin = Length(Cross(from, to));

	result.m[0][0] = normal.x * normal.x * (1.0f - cos) + cos;
	result.m[0][1] = normal.x * normal.y * (1.0f - cos) + normal.z * sin;
	result.m[0][2] = normal.x * normal.z * (1.0f - cos) - normal.y * sin;

	result.m[1][0] = normal.x * normal.y * (1.0f - cos) - normal.z * sin;
	result.m[1][1] = normal.y * normal.y * (1.0f - cos) + cos;
	result.m[1][2] = normal.y * normal.z * (1.0f - cos) + normal.x * sin;

	result.m[2][0] = normal.x * normal.z * (1.0f - cos) + normal.y * sin;
	result.m[2][1] = normal.y * normal.z * (1.0f - cos) - normal.x * sin;
	result.m[2][2] = normal.z * normal.z * (1.0f - cos) + cos;

	return result;
}

Vector3 Add(const Vector3& v1, const Vector3& v2) {

	Vector3 result;

	result = { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };

	return result;

}

Quaternion Multiply(const Quaternion& lhs, const Quaternion& rhs)
{

	Quaternion result = {};

	Vector3 q = { lhs.x, lhs.y, lhs.z, };
	Vector3 r = { rhs.x, rhs.y, rhs.z, };
	Vector3 resultVector = Add(Cross(q, r),Add(Multiply(rhs.w, q),Multiply(lhs.w, r)));

	result.x = resultVector.x;
	result.y = resultVector.y;
	result.z = resultVector.z;
	result.w = lhs.w * rhs.w - Dot(q, r);

	return result;
}

Quaternion IdentityQuaternion()
{
	Quaternion result = { 0.0f,0.0f,0.0f,1.0f };
	return result;
}

Quaternion Conjugate(const Quaternion& quaternion)
{
	Quaternion result = { -quaternion.x, -quaternion.y, -quaternion.z, quaternion.w };
	return result;
}

float Norm(const Quaternion& quaternion)
{
	float result = sqrtf(quaternion.x * quaternion.x +
		quaternion.y * quaternion.y +
		quaternion.z * quaternion.z +
		quaternion.w * quaternion.w);

	return result;
}

Quaternion Normalize(const Quaternion& quaternion)
{
	Quaternion result = {};
	float norm = Norm(quaternion);

	if (quaternion.x != 0.0) {
		result.x = quaternion.x / norm;
	}
	else {
		result.x = 0.0f;
	}

	if (quaternion.y != 0.0) {
		result.y = quaternion.y / norm;
	}
	else {
		result.y = 0.0f;
	}

	if (quaternion.z != 0.0) {
		result.z = quaternion.z / norm;
	}
	else {
		result.z = 0.0f;
	}

	if (quaternion.w != 0.0) {
		result.w = quaternion.w / norm;
	}
	else {
		result.w = 0.0f;
	}

	return result;

}

Quaternion Inverse(const Quaternion& quaternion)
{

	Quaternion result = {};
	Quaternion conjugate = Conjugate(quaternion);

	float norm = Norm(quaternion);

	if (norm != 0.0f) {
		result.x = conjugate.x / (norm * norm);
		result.y = conjugate.y / (norm * norm);
		result.z = conjugate.z / (norm * norm);
		result.w = conjugate.w / (norm * norm);
	}

	return result;
}
// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	Quaternion q1 = { 2.0f,3.0f,4.0f,1.0f };
	Quaternion q2 = { 1.0f,3.0f,5.0f,2.0f };
	Quaternion identity = IdentityQuaternion();
	Quaternion conj = Conjugate(q1);
	Quaternion inv = Inverse(q1);
	Quaternion normal = Normalize(q1);
	//Quaternion mul1 = Multiply(q1, q2);
	//Quaternion mul2 = Multiply(q2, q1);
	//float norm = Norm(q1);

	// ウィンドウの×ボタンが押されるまでループ
	while (Novice::ProcessMessage() == 0) {
		// フレームの開始
		Novice::BeginFrame();

		// キー入力を受け取る
		memcpy(preKeys, keys, 256);
		Novice::GetHitKeyStateAll(keys);

		///
		/// ↓更新処理ここから
		///

		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///
		Novice::ScreenPrintf(20, 20, "%f,%f,%f,%f :Identity", identity.x,identity.y,identity.z,identity.w);
		Novice::ScreenPrintf(20, 40, "%f,%f,%f,%f : Conjugate", conj.x,conj.y,conj.z,conj.w);
		Novice::ScreenPrintf(20, 60, "%f,%f,%f,%f : Inverse", inv.x,inv.y,inv.z,inv.w);
		Novice::ScreenPrintf(20, 80, "%f,%f,%f,%f : Normalize",normal.z,normal.y,normal.z,normal.w);
		Novice::ScreenPrintf(20, 100, "%f,%f,%f,%f : Multiply(q1, q2)",mul1.x,mul1.y,mul1.z,mul1.w);
		Novice::ScreenPrintf(20, 120, "%f,%f,%f,%f : Multiply(q2, q1)",mul2.x,mul2.y,mul2.z,mul2);
		Novice::ScreenPrintf(20, 140, "%f : Norm",norm);


		///
		/// ↑描画処理ここまで
		///

		// フレームの終了
		Novice::EndFrame();

		// ESCキーが押されたらループを抜ける
		if (preKeys[DIK_ESCAPE] == 0 && keys[DIK_ESCAPE] != 0) {
			break;
		}
	}

	// ライブラリの終了
	Novice::Finalize();
	return 0;
}
