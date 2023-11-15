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

struct Quaternion final {
	//float x;
	//float y;
	//float z;
	Vector3 v;
	float w;
};

Vector2 Add(const Vector2& v1, const Vector2& v2) {
	Vector2 v;
	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	return v;
}

Vector2 Subtruct(const Vector2& v1, const Vector2& v2) {
	Vector2 v;
	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	return v;
}

Vector2 Multiply(float scalar, const Vector2& v) {
	Vector2 returnV;
	returnV.x = v.x * scalar;
	returnV.y = v.y * scalar;
	return returnV;
}

Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 v;
	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	v.z = v1.z + v2.z;
	return v;
}

Vector3 Subtruct(const Vector3& v1, const Vector3& v2) {
	Vector3 v;
	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	v.z = v1.z - v2.z;
	return v;
}

Vector3 Multiply(float scalar, const Vector3& v) {
	Vector3 returnV;
	returnV.x = v.x * scalar;
	returnV.y = v.y * scalar;
	returnV.z = v.z * scalar;
	return returnV;
}

float Dot(const Vector3& v1, const Vector3& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }

float Length(const Vector3& v) { return sqrtf(Dot(v, v)); }

Vector3 Normalize(const Vector3& v) {
	float length = Length(v);
	if (length == 0) {
		return v;
	}
	return Multiply((1.0f / length), v);
}

Vector3 Cross(const Vector3& v1, const Vector3& v2)
{
	return { v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x };
}

Vector2 Lerp(const Vector2& v1, const Vector2& v2, float t) { return v1 + Multiply(t, v2 - v1); }

Vector3 Lerp(const Vector3& v1, const Vector3& v2, float t)
{
	return v1 + Multiply(t, v2 - v1);
}


Vector3 Slerp(const Vector3& v1, const Vector3& v2, float t)
{
	Vector3 a = Normalize(v1), b = Normalize(v2);
	float s = (1.0f - t) * Length(a) + t * Length(b);
	Vector3 e1, e2;
	e1 = float(1.0f / Length(a)) * a;
	e2 = float(1.0f / Length(b)) * b;

	float dot = std::clamp(Dot(a, b), 0.0f, 1.0f);
	float theta = std::acos(dot/*/( Length(a)*Length(b))*/);
	if (theta == 0.0f)
	{
		return Lerp(a, b, t);
	}
	return s * ((std::sinf((1.0f - t) * theta) / std::sinf(theta)) * a + (std::sinf(t * theta) / std::sinf(theta)) * b);
}


Quaternion Multiply(const Quaternion& q, const Quaternion& r) {
	Quaternion result;
	result.w = q.w * r.w - Dot(q.v, r.v);
	result.v = Cross(q.v, r.v) + r.w * q.v + q.w * r.v;
	return result;
}

Quaternion IdentityQuaternion() {
	return Quaternion{ 1.0f,0,0,0 };
}

Quaternion Conjugate(const Quaternion& q) {
	return Quaternion{ (-1.0f * q.v),q.w };
}

float Norm(const Quaternion& q) {
	return std::sqrtf(q.w * q.w + q.v.x * q.v.x + q.v.y * q.v.y + q.v.z * q.v.z);
}

Quaternion Normalize(const Quaternion& q) {
	Quaternion result = q;
	float norm = Norm(q);
	result.w /= norm;
	result.v.x /= norm;
	result.v.y /= norm;
	result.v.z /= norm;
	return result;
}

Quaternion Inverse(const Quaternion& q) {
	Quaternion result = Conjugate(q);
	float length = Norm(q) * Norm(q);
	result.w /= length;
	result.v.x /= length;
	result.v.y /= length;
	result.v.z /= length;
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
	Quaternion mul1 = Multiply(q1, q2);
	Quaternion mul2 = Multiply(q2, q1);
	float norm = Norm(q1);

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
