#include <Novice.h>
#include <cmath>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <imgui.h>

const char kWindowTitle[] = "LD2B_01_オノ_フミカ_MT3_";

// 3次元ベクトルの構造体
struct Vector3 { float x, y, z; };

// 4ｘ4行列の構造体
struct Matrix4x4 {
	float m[4][4];
};

// 球
struct Sphere {
	Vector3 center; // 中心点
	float radius; // 半径
};

// 円錐振り子
struct ConicalPendulm {
	// アンカーポイント
	Vector3 anchor;
	// ヒモの長さ
	float length;
	// 円錐の頂角の半分
	float halfApexAngle;
	// 現在の角度
	float angle;
	// 角速度ω
	float anglarVelocity;
};

// 加算
Vector3 Add(const Vector3& v1, const Vector3& v2) {
	Vector3 result = {};

	result.x = v1.x + v2.x;
	result.y = v1.y + v2.y;
	result.z = v1.z + v2.z;

	return result;
};
// 減産
Vector3 Subtract(const Vector3& v1, const Vector3& v2) {
	Vector3 result = {};

	result.x = v1.x - v2.x;
	result.y = v1.y - v2.y;
	result.z = v1.z - v2.z;

	return result;
};
// 割り算
Vector3 Division(float scalar, const Vector3& v) {
	Vector3 result = {};

	result.x = v.x / scalar;
	result.y = v.y / scalar;
	result.z = v.z / scalar;

	return result;
}
// スカラー倍
Vector3 Multiply(float scalar, const Vector3& v) {
	Vector3 result = {};

	result.x = v.x * scalar;
	result.y = v.y * scalar;
	result.z = v.z * scalar;

	return result;
};
// 内席
float Dot(const Vector3& v1, const Vector3& v2) {
	float result = {};
	result = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);

	return result;
};
// 長さ（ノルム）
float Length(const Vector3& v) {
	float result = {};

	result = sqrtf(Dot(v, v));

	return result;
};
// 正規化
Vector3 Normalize(const Vector3& v) {
	Vector3 result = {};

	result.x = (v.x / Length(v));
	result.y = (v.y / Length(v));
	result.z = (v.z / Length(v));

	return result;
};
// 座標変換
Vector3 Transform(const Vector3& vector, const Matrix4x4& matrix) {
	// ｗ = 1 がデカルト座標系であるので（ｘ、ｙ、z,1）のベクトルとしてmatrixとの積をとる
	Vector3 result = {};

	result.x =
		vector.x * matrix.m[0][0] +
		vector.y * matrix.m[1][0] +
		vector.z * matrix.m[2][0] +
		1.0f * matrix.m[3][0];

	result.y =
		vector.x * matrix.m[0][1] +
		vector.y * matrix.m[1][1] +
		vector.z * matrix.m[2][1] +
		1.0f * matrix.m[3][1];

	result.z =
		vector.x * matrix.m[0][2] +
		vector.y * matrix.m[1][2] +
		vector.z * matrix.m[2][2] +
		1.0f * matrix.m[3][2];

	float w =
		vector.x * matrix.m[0][3] +
		vector.y * matrix.m[1][3] +
		vector.z * matrix.m[2][3] +
		1.0f * matrix.m[3][3];

	// ベクトルに対して基本的な操作を行う行列でｗが0になることは無い
	assert(w != 0.0f);
	// ｗ＝1がデカルト座標系であるので、ｗ除算することで同次座標をデカルト座標に戻す
	result.x /= w;
	result.y /= w;
	result.z /= w;

	return result;
};
// X軸回転行列
Matrix4x4 MakeRotateXMatrix(float radian) {
	Matrix4x4 result = {};

	result.m[0][0] = 1;
	result.m[1][1] = std::cos(radian);
	result.m[1][2] = std::sin(radian);
	result.m[2][1] = -std::sin(radian);
	result.m[2][2] = std::cos(radian);
	result.m[3][3] = 1;

	return result;
};
// Y軸回転行列
Matrix4x4 MakeRotateYMatrix(float radian) {
	Matrix4x4 result = {};

	result.m[0][0] = std::cos(radian);
	result.m[0][2] = -std::sin(radian);
	result.m[1][1] = 1;
	result.m[2][0] = std::sin(radian);
	result.m[2][2] = std::cos(radian);
	result.m[3][3] = 1;

	return result;
};
// Z軸回転行列
Matrix4x4 MakeRotateZMatrix(float radian) {
	Matrix4x4 result = {};

	result.m[0][0] = std::cos(radian);
	result.m[0][1] = std::sin(radian);
	result.m[1][0] = -std::sin(radian);;
	result.m[1][1] = std::cos(radian);
	result.m[2][2] = 1;
	result.m[3][3] = 1;

	return result;
};
// 行列の積
Matrix4x4 Multiply(const Matrix4x4& m1, const Matrix4x4& m2) {
	Matrix4x4 result = {};

	for (int x = 0; x < 4; x++) {
		for (int y = 0; y < 4; y++) {
			for (int a = 0; a < 4; a++) {
				result.m[x][y] += m1.m[x][a] * m2.m[a][y];
			}
		}
	}

	return result;
};
// アフィン変換
Matrix4x4 MakeAffineMatrix(const Vector3& scale, const Vector3& rotate, const Vector3& translate) {
	Matrix4x4 result = {};

	Matrix4x4 rotateXMatrix = MakeRotateXMatrix(rotate.x);
	Matrix4x4 rotateYMatrix = MakeRotateYMatrix(rotate.y);
	Matrix4x4 rotateZMatrix = MakeRotateZMatrix(rotate.z);
	Matrix4x4 rotateXYZMatrix = Multiply(rotateXMatrix, Multiply(rotateYMatrix, rotateZMatrix));


	result.m[0][0] = scale.x * rotateXYZMatrix.m[0][0];
	result.m[0][1] = scale.x * rotateXYZMatrix.m[0][1];
	result.m[0][2] = scale.x * rotateXYZMatrix.m[0][2];

	result.m[1][0] = scale.y * rotateXYZMatrix.m[1][0];
	result.m[1][1] = scale.y * rotateXYZMatrix.m[1][1];
	result.m[1][2] = scale.y * rotateXYZMatrix.m[1][2];

	result.m[2][0] = scale.z * rotateXYZMatrix.m[2][0];
	result.m[2][1] = scale.z * rotateXYZMatrix.m[2][1];
	result.m[2][2] = scale.z * rotateXYZMatrix.m[2][2];

	result.m[3][0] = translate.x;
	result.m[3][1] = translate.y;
	result.m[3][2] = translate.z;
	result.m[3][3] = 1;

	return result;
};
// 逆行列
Matrix4x4 Inverse(const Matrix4x4& m) {
	Matrix4x4 result = {};
	float determinant = 0;

	// 行列式 |A| を求める
	determinant =
		(m.m[0][0] * m.m[1][1] * m.m[2][2] * m.m[3][3]) + (m.m[0][0] * m.m[1][2] * m.m[2][3] * m.m[3][1]) + (m.m[0][0] * m.m[1][3] * m.m[2][1] * m.m[3][2]) -
		(m.m[0][0] * m.m[1][3] * m.m[2][2] * m.m[3][1]) - (m.m[0][0] * m.m[1][2] * m.m[2][1] * m.m[3][3]) - (m.m[0][0] * m.m[1][1] * m.m[2][3] * m.m[3][2]) -
		(m.m[0][1] * m.m[1][0] * m.m[2][2] * m.m[3][3]) - (m.m[0][2] * m.m[1][0] * m.m[2][3] * m.m[3][1]) - (m.m[0][3] * m.m[1][0] * m.m[2][1] * m.m[3][2]) +
		(m.m[0][3] * m.m[1][0] * m.m[2][2] * m.m[3][1]) + (m.m[0][2] * m.m[1][0] * m.m[2][1] * m.m[3][3]) + (m.m[0][1] * m.m[1][0] * m.m[2][3] * m.m[3][2]) +
		(m.m[0][1] * m.m[1][2] * m.m[2][0] * m.m[3][3]) + (m.m[0][2] * m.m[1][3] * m.m[2][0] * m.m[3][1]) + (m.m[0][3] * m.m[1][1] * m.m[2][0] * m.m[3][2]) -
		(m.m[0][3] * m.m[1][2] * m.m[2][0] * m.m[3][1]) - (m.m[0][2] * m.m[1][1] * m.m[2][0] * m.m[3][3]) - (m.m[0][1] * m.m[1][3] * m.m[2][0] * m.m[3][2]) -
		(m.m[0][1] * m.m[1][2] * m.m[2][3] * m.m[3][0]) - (m.m[0][2] * m.m[1][3] * m.m[2][1] * m.m[3][0]) - (m.m[0][3] * m.m[1][1] * m.m[2][2] * m.m[3][0]) +
		(m.m[0][3] * m.m[1][2] * m.m[2][1] * m.m[3][0]) + (m.m[0][2] * m.m[1][1] * m.m[2][3] * m.m[3][0]) + (m.m[0][1] * m.m[1][3] * m.m[2][2] * m.m[3][0]);

	// 逆行列を求める
	result.m[0][0] = (1 / determinant) *
		((m.m[1][1] * m.m[2][2] * m.m[3][3]) + (m.m[1][2] * m.m[2][3] * m.m[3][1]) + (m.m[1][3] * m.m[2][1] * m.m[3][2]) -
			(m.m[1][3] * m.m[2][2] * m.m[3][1]) - (m.m[1][2] * m.m[2][1] * m.m[3][3]) - (m.m[1][1] * m.m[2][3] * m.m[3][2]));
	result.m[0][1] = (1 / determinant) *
		(-(m.m[0][1] * m.m[2][2] * m.m[3][3]) - (m.m[0][2] * m.m[2][3] * m.m[3][1]) - (m.m[0][3] * m.m[2][1] * m.m[3][2]) +
			(m.m[0][3] * m.m[2][2] * m.m[3][1]) + (m.m[0][2] * m.m[2][1] * m.m[3][3]) + (m.m[0][1] * m.m[2][3] * m.m[3][2]));
	result.m[0][2] = (1 / determinant) *
		((m.m[0][1] * m.m[1][2] * m.m[3][3]) + (m.m[0][2] * m.m[1][3] * m.m[3][1]) + (m.m[0][3] * m.m[1][1] * m.m[3][2]) -
			(m.m[0][3] * m.m[1][2] * m.m[3][1]) - (m.m[0][2] * m.m[1][1] * m.m[3][3]) - (m.m[0][1] * m.m[1][3] * m.m[3][2]));
	result.m[0][3] = (1 / determinant) *
		(-(m.m[0][1] * m.m[1][2] * m.m[2][3]) - (m.m[0][2] * m.m[1][3] * m.m[2][1]) - (m.m[0][3] * m.m[1][1] * m.m[2][2]) +
			(m.m[0][3] * m.m[1][2] * m.m[2][1]) + (m.m[0][2] * m.m[1][1] * m.m[2][3]) + (m.m[0][1] * m.m[1][3] * m.m[2][2]));

	result.m[1][0] = (1 / determinant) *
		(-(m.m[1][0] * m.m[2][2] * m.m[3][3]) - (m.m[1][2] * m.m[2][3] * m.m[3][0]) - (m.m[1][3] * m.m[2][0] * m.m[3][2]) +
			(m.m[1][3] * m.m[2][2] * m.m[3][0]) + (m.m[1][2] * m.m[2][0] * m.m[3][3]) + (m.m[1][0] * m.m[2][3] * m.m[3][2]));
	result.m[1][1] = (1 / determinant) *
		((m.m[0][0] * m.m[2][2] * m.m[3][3]) + (m.m[0][2] * m.m[2][3] * m.m[3][0]) + (m.m[0][3] * m.m[2][0] * m.m[3][2]) -
			(m.m[0][3] * m.m[2][2] * m.m[3][0]) - (m.m[0][2] * m.m[2][0] * m.m[3][3]) - (m.m[0][0] * m.m[2][3] * m.m[3][2]));
	result.m[1][2] = (1 / determinant) *
		(-(m.m[0][0] * m.m[1][2] * m.m[3][3]) - (m.m[0][2] * m.m[1][3] * m.m[3][0]) - (m.m[0][3] * m.m[1][0] * m.m[3][2]) +
			(m.m[0][3] * m.m[1][2] * m.m[3][0]) + (m.m[0][2] * m.m[1][0] * m.m[3][3]) + (m.m[0][0] * m.m[1][3] * m.m[3][2]));
	result.m[1][3] = (1 / determinant) *
		((m.m[0][0] * m.m[1][2] * m.m[2][3]) + (m.m[0][2] * m.m[1][3] * m.m[2][0]) + (m.m[0][3] * m.m[1][0] * m.m[2][2]) -
			(m.m[0][3] * m.m[1][2] * m.m[2][0]) - (m.m[0][2] * m.m[1][0] * m.m[2][3]) - (m.m[0][0] * m.m[1][3] * m.m[2][2]));

	result.m[2][0] = (1 / determinant) *
		((m.m[1][0] * m.m[2][1] * m.m[3][3]) + (m.m[1][1] * m.m[2][3] * m.m[3][0]) + (m.m[1][3] * m.m[2][0] * m.m[3][1]) -
			(m.m[1][3] * m.m[2][1] * m.m[3][0]) - (m.m[1][1] * m.m[2][0] * m.m[3][3]) - (m.m[1][0] * m.m[2][3] * m.m[3][1]));
	result.m[2][1] = (1 / determinant) *
		(-(m.m[0][0] * m.m[2][1] * m.m[3][3]) - (m.m[0][1] * m.m[2][3] * m.m[3][0]) - (m.m[0][3] * m.m[2][0] * m.m[3][1]) +
			(m.m[0][3] * m.m[2][1] * m.m[3][0]) + (m.m[0][1] * m.m[2][0] * m.m[3][3]) + (m.m[0][0] * m.m[2][3] * m.m[3][1]));
	result.m[2][2] = (1 / determinant) *
		((m.m[0][0] * m.m[1][1] * m.m[3][3]) + (m.m[0][1] * m.m[1][3] * m.m[3][0]) + (m.m[0][3] * m.m[1][0] * m.m[3][1]) -
			(m.m[0][3] * m.m[1][1] * m.m[3][0]) - (m.m[0][1] * m.m[1][0] * m.m[3][3]) - (m.m[0][0] * m.m[1][3] * m.m[3][1]));
	result.m[2][3] = (1 / determinant) *
		(-(m.m[0][0] * m.m[1][1] * m.m[2][3]) - (m.m[0][1] * m.m[1][3] * m.m[2][0]) - (m.m[0][3] * m.m[1][0] * m.m[2][1]) +
			(m.m[0][3] * m.m[1][1] * m.m[2][0]) + (m.m[0][1] * m.m[1][0] * m.m[2][3]) + (m.m[0][0] * m.m[1][3] * m.m[2][1]));

	result.m[3][0] = (1 / determinant) *
		(-(m.m[1][0] * m.m[2][1] * m.m[3][2]) - (m.m[1][1] * m.m[2][2] * m.m[3][0]) - (m.m[1][2] * m.m[2][0] * m.m[3][1]) +
			(m.m[1][2] * m.m[2][1] * m.m[3][0]) + (m.m[1][1] * m.m[2][0] * m.m[3][2]) + (m.m[1][0] * m.m[2][2] * m.m[3][1]));
	result.m[3][1] = (1 / determinant) *
		((m.m[0][0] * m.m[2][1] * m.m[3][2]) + (m.m[0][1] * m.m[2][2] * m.m[3][0]) + (m.m[0][2] * m.m[2][0] * m.m[3][1]) -
			(m.m[0][2] * m.m[2][1] * m.m[3][0]) - (m.m[0][1] * m.m[2][0] * m.m[3][2]) - (m.m[0][0] * m.m[2][2] * m.m[3][1]));
	result.m[3][2] = (1 / determinant) *
		(-(m.m[0][0] * m.m[1][1] * m.m[3][2]) - (m.m[0][1] * m.m[1][2] * m.m[3][0]) - (m.m[0][2] * m.m[1][0] * m.m[3][1]) +
			(m.m[0][2] * m.m[1][1] * m.m[3][0]) + (m.m[0][1] * m.m[1][0] * m.m[3][2]) + (m.m[0][0] * m.m[1][2] * m.m[3][1]));
	result.m[3][3] = (1 / determinant) *
		((m.m[0][0] * m.m[1][1] * m.m[2][2]) + (m.m[0][1] * m.m[1][2] * m.m[2][0]) + (m.m[0][2] * m.m[1][0] * m.m[2][1]) -
			(m.m[0][2] * m.m[1][1] * m.m[2][0]) - (m.m[0][1] * m.m[1][0] * m.m[2][2]) - (m.m[0][0] * m.m[1][2] * m.m[2][1]));

	return result;
};
// 透視投影行列
Matrix4x4 MakePerspectiveFovMatrix(float fovY, float aspectRatio, float nearClip, float farClip) {
	Matrix4x4 result = {};
	float cot = 1 / (std::tan(fovY / 2));

	result.m[0][0] = (1 / aspectRatio) * cot;
	result.m[1][1] = cot;
	result.m[2][2] = farClip / (farClip - nearClip);
	result.m[2][3] = 1;
	result.m[3][2] = -nearClip * (farClip / (farClip - nearClip));

	return result;
};
// ビューポート変換行列
Matrix4x4 MakeViewportMatrix(float left, float top, float width, float height, float minDepth, float maxDepth) {
	Matrix4x4 result = {};

	result.m[0][0] = width / 2;
	result.m[1][1] = -(height / 2);
	result.m[2][2] = maxDepth - minDepth;
	result.m[3][0] = left + (width / 2);
	result.m[3][1] = top + (height / 2);
	result.m[3][2] = minDepth;
	result.m[3][3] = 1;

	return result;
};

// グリッド
void DrawGrid(const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewportMatrix) {
	const float kGridHalfWidth = 2.0f;// 半分
	const uint32_t kSubdivision = 10;// 分割数
	const float kGridEvery = (kGridHalfWidth * 2.0f) / float(kSubdivision);// 1つ分の長さ

	Vector3 zLineStart;
	Vector3 zLineEnd;
	Vector3 xLineStart;
	Vector3 xLineEnd;

	// 奥から手前の線を順々に引いていく
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		// 視点と頂点を求める
		zLineStart = Vector3(xIndex * kGridEvery - kGridHalfWidth, 0.0f, kGridHalfWidth);
		zLineEnd = Vector3(xIndex * kGridEvery - kGridHalfWidth, 0.0f, -kGridHalfWidth);

		// 変換
		Matrix4x4 worldMatrixStart = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, zLineStart);
		Matrix4x4 worldMatrixEnd = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, zLineEnd);

		Matrix4x4 wvpMatrixStart = Multiply(worldMatrixStart, viewProjectionMatrix);
		Matrix4x4 wvpMatrixEnd = Multiply(worldMatrixEnd, viewProjectionMatrix);

		Vector3 nbcVertexStart = Transform(Vector3{}, wvpMatrixStart);
		Vector3 screenStartPoint = Transform(nbcVertexStart, viewportMatrix);
		Vector3 nbcVertexEnd = Transform(Vector3{}, wvpMatrixEnd);
		Vector3 screenEndPoint = Transform(nbcVertexEnd, viewportMatrix);

		// 描画
		if (xIndex == 5) {
			Novice::DrawLine((int)screenStartPoint.x, (int)screenStartPoint.y, (int)screenEndPoint.x, (int)screenEndPoint.y, BLACK);
		}
		else {
			Novice::DrawLine((int)screenStartPoint.x, (int)screenStartPoint.y, (int)screenEndPoint.x, (int)screenEndPoint.y, 0xAAAAAAFF);
		}
	}

	// 左から右
	for (uint32_t xIndex = 0; xIndex <= kSubdivision; ++xIndex) {
		// 視点と頂点を求める
		xLineStart = Vector3(kGridHalfWidth, 0.0f, xIndex * kGridEvery - kGridHalfWidth);
		xLineEnd = Vector3(-kGridHalfWidth, 0.0f, xIndex * kGridEvery - kGridHalfWidth);

		// 変換
		Matrix4x4 worldMatrixStart = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, xLineStart);
		Matrix4x4 worldMatrixEnd = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, xLineEnd);

		Matrix4x4 wvpMatrixStart = Multiply(worldMatrixStart, viewProjectionMatrix);
		Matrix4x4 wvpMatrixEnd = Multiply(worldMatrixEnd, viewProjectionMatrix);

		Vector3 nbcVertexStart = Transform(Vector3{}, wvpMatrixStart);
		Vector3 screenStartPoint = Transform(nbcVertexStart, viewportMatrix);
		Vector3 nbcVertexEnd = Transform(Vector3{}, wvpMatrixEnd);
		Vector3 screenEndPoint = Transform(nbcVertexEnd, viewportMatrix);

		// 描画
		if (xIndex == 5) {
			Novice::DrawLine((int)screenStartPoint.x, (int)screenStartPoint.y, (int)screenEndPoint.x, (int)screenEndPoint.y, BLACK);
		}
		else {
			Novice::DrawLine((int)screenStartPoint.x, (int)screenStartPoint.y, (int)screenEndPoint.x, (int)screenEndPoint.y, 0xAAAAAAFF);
		}
	}

}
// 球
void DrawSphere(const Sphere& sphere, const Matrix4x4& viewProjectionMatrix, const Matrix4x4& viewPortMatrix, uint32_t color) {
	float pi = float(M_PI);
	const uint32_t kSudivision = 10; //分割数
	const float kLatEvery = pi / kSudivision; //緯度分割一つ分の角度
	const float kLonEvery = (pi * 2) / kSudivision; //緯度分割二つ分の角度
	// world座標系でのa,b,cを求める
	Vector3 a, b, c;

	// 緯度の方向に分割
	for (uint32_t latindex = 0; latindex < kSudivision; ++latindex) {
		sphere;
		float lat = -pi / 2.0f + kLatEvery * latindex; //現在の緯度

		// 経度の方向に分割
		for (uint32_t lonIndex = 0; lonIndex < kSudivision; ++lonIndex) {
			float lon = lonIndex * kLonEvery; //現在の経度
			// ワールド座標系での頂点を求める
			a = {
				sphere.center.x + sphere.radius * (cosf(lat) * cosf(lon)),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * (cosf(lat) * sinf(lon))
			};
			b = {
				sphere.center.x + sphere.radius * (cosf(lat + kLatEvery) * cosf(lon)),
				sphere.center.y + sphere.radius * sinf(lat + kLatEvery),
				sphere.center.z + sphere.radius * (cosf(lat + kLatEvery) * sin(lon))
			};
			c = {
				sphere.center.x + sphere.radius * (cosf(lat) * cosf(lon + kLonEvery)),
				sphere.center.y + sphere.radius * sinf(lat),
				sphere.center.z + sphere.radius * (cos(lat) * sinf(lon + kLonEvery))
			};

			// 変換
			Matrix4x4 worldMatrixA = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, a);
			Matrix4x4 worldMatrixB = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, b);
			Matrix4x4 worldMatrixC = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, { 0.0f,0.0f,0.0f }, c);

			Matrix4x4 wvpMatrixA = Multiply(worldMatrixA, viewProjectionMatrix);
			Matrix4x4 wvpMatrixB = Multiply(worldMatrixB, viewProjectionMatrix);
			Matrix4x4 wvpMatrixC = Multiply(worldMatrixC, viewProjectionMatrix);

			Vector3 nbcVertexA = Transform(Vector3{}, wvpMatrixA);
			Vector3 nbcVertexB = Transform(Vector3{}, wvpMatrixB);
			Vector3 nbcVertexC = Transform(Vector3{}, wvpMatrixC);

			Vector3 screenA = Transform(nbcVertexA, viewPortMatrix);
			Vector3 screenB = Transform(nbcVertexB, viewPortMatrix);
			Vector3 screenC = Transform(nbcVertexC, viewPortMatrix);

			// 描画
			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenB.x, (int)screenB.y, color);
			Novice::DrawLine((int)screenA.x, (int)screenA.y, (int)screenC.x, (int)screenC.y, color);


		}
	}
}

// Windowsアプリでのエントリーポイント(main関数)
int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int) {

	// ライブラリの初期化
	Novice::Initialize(kWindowTitle, 1280, 720);

	// キー入力結果を受け取る箱
	char keys[256] = {0};
	char preKeys[256] = {0};

	Vector3 cameraRotate{ 0.26f,0.0f,0.0f };
	Vector3 cameraTranslate{ 0.0f,1.9f,-6.49f };
	int kWindowWidth = 1280;
	int kWindowHeight = 720;

	// 円錐振り子
	ConicalPendulm conicalPendulm;
	conicalPendulm.anchor = { 0.0f,1.0f,0.0f };
	conicalPendulm.length = 0.8f;
	conicalPendulm.halfApexAngle = 0.7f;
	conicalPendulm.angle = 0.0f;
	conicalPendulm.anglarVelocity = 0.0f;
	// 振り子の先端
	Vector3 position = {};

	// 球
	Sphere sphere = {};

	// デルタタイム
	float deltaTime = 1.0f / 60.0f;
	// 揺れているかのフラグ
	bool isRotate = false;

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

		if (ImGui::Button("Start")) {
			isRotate = true;
		}

		if (isRotate) {
			// 角速度を計算する
			conicalPendulm.anglarVelocity = std::sqrt(9.8f / (conicalPendulm.length * std::cos(conicalPendulm.halfApexAngle)));
			conicalPendulm.angle += conicalPendulm.anglarVelocity * deltaTime;
			float radius = std::sin(conicalPendulm.halfApexAngle) * conicalPendulm.length;
			float heigth = std::cos(conicalPendulm.halfApexAngle) * conicalPendulm.length;
			position.x = conicalPendulm.anchor.x + std::cos(conicalPendulm.angle) * radius;
			position.y = conicalPendulm.anchor.y - heigth;
			position.z = conicalPendulm.anchor.z - std::sin(conicalPendulm.angle) * radius;
		}

		// カメラ
		Matrix4x4 cameraWorldMatrix = MakeAffineMatrix({ 1.0f,1.0f,1.0f }, cameraRotate, cameraTranslate);
		Matrix4x4 viewMatrix = Inverse(cameraWorldMatrix);
		Matrix4x4 projectionMatrix = MakePerspectiveFovMatrix(0.45f, float(kWindowWidth) / float(kWindowHeight), 0.1f, 100.0f);
		Matrix4x4 viewProjectionMatrix = Multiply(viewMatrix, projectionMatrix);
		Matrix4x4 viewportMatrix = MakeViewportMatrix(0, 0, float(kWindowWidth), float(kWindowHeight), 0.0f, 1.0f);

		sphere = { position,0.05f };

		// 座標変換
		Vector3 screenPosition = Transform(Transform(sphere.center, viewProjectionMatrix), viewportMatrix);
		Vector3 screenAnchor = Transform(Transform(conicalPendulm.anchor, viewProjectionMatrix), viewportMatrix);


		///
		/// ↑更新処理ここまで
		///

		///
		/// ↓描画処理ここから
		///

		// グリッド描画
		DrawGrid(viewProjectionMatrix, viewportMatrix);

		// 振り子のボール
		DrawSphere(sphere, viewProjectionMatrix, viewportMatrix, WHITE);
		// 振り子のヒモ
		Novice::DrawLine((int)screenAnchor.x, (int)screenAnchor.y, (int)screenPosition.x, (int)screenPosition.y, WHITE);


		ImGui::DragFloat("length", &conicalPendulm.length, 0.01f);
		ImGui::DragFloat("halfApexAngle", &conicalPendulm.halfApexAngle, 0.01f);

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
