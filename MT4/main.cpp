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

const char kWindowTitle[] = "学籍番号";

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
