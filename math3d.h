#include <math.h>

#ifndef MATH3D_H
#define MATH3D_H

#define PI 3.14159265f

struct Matrix;

float DegToRad(float deg);
float RadToDeg(float rad);

struct Vector
{
	Vector() {}
	const float& operator[](int i) const { return (&x)[i]; }
	float& operator[](int i) { return (&x)[i]; }
	
	float x;
	float y;
	float z;
	float w;
};

float Dot2(const Vector& a, const Vector& b);
float Dot3(const Vector& a, const Vector& b);
float Dot4(const Vector& a, const Vector& b);

struct Vector2 : public Vector
{
	Vector2() {}
	Vector2(float xIn, float yIn) { x = xIn; y = yIn; }
	Vector2& operator+=(const Vector2& v);
	void Set(float xIn, float yIn);
};

Vector2 operator-(const Vector2& v);
Vector2 operator-(const Vector2& a, const Vector2& b);
Vector2 operator+(const Vector2& a, const Vector2& b);
Vector2 operator*(const Vector2& v, float factor);
Vector2 operator*(float factor, const Vector2& v);
float operator*(const Vector2& a, const Vector2& b);
float Len(const Vector2& v);
Vector2 UnitVec(const Vector2& v);

struct Vector3 : public Vector
{
	Vector3() {}
	Vector3(const Vector2& v, float zIn) { x = v.x; y = v.y; z = zIn; }
	Vector3(float xIn, float yIn, float zIn) { x = xIn; y = yIn; z = zIn; }
	void Set(float xIn, float yIn, float zIn);
	void Normalize();
	void Mul3x4(const Matrix& m);
};

Vector3 operator-(const Vector3& v);
Vector3 operator-(const Vector3& a, const Vector3& b);
Vector3 operator+(const Vector3& a, const Vector3& b);
Vector3 operator*(const Vector3& v, float factor);
Vector3 operator*(float factor, const Vector3& v);
float operator*(const Vector3& a, const Vector3& b);
float Len(const Vector3& v);
Vector3 UnitVec(const Vector3& v);
Vector3 operator%(const Vector3& a, const Vector3& b);	// cross product

struct Vector4 : public Vector
{
	Vector4() {}
	Vector4(const Vector3& v, float wIn) { x = v.x; y = v.y; z = v.z; w = wIn; }
	Vector4(float xIn, float yIn, float zIn, float wIn) { x = xIn; y = yIn; z = zIn; w = wIn; }
	void Set(float xIn, float yIn, float zIn, float wIn);
};

Vector4 operator-(const Vector4& v);
Vector4 operator-(const Vector4& a, const Vector4& b);
Vector4 operator+(const Vector4& a, const Vector4& b);
Vector4 operator*(const Vector4& v, float factor);
Vector4 operator*(float factor, const Vector4& v);
float operator*(const Vector4& a, const Vector4& b);
float Len(const Vector4& v);
Vector4 UnitVec(const Vector4& v);

struct Quat : public Vector4
{
	Quat() {}
	Quat(float xIn, float yIn, float zIn, float wIn) { x = xIn; y = yIn; z = zIn; w = wIn; }
	Quat(const Vector3& axis, float angle);
	
	Vector3 Rotate(const Vector3& v) const;
};

Quat operator~(const Quat& v);	// conjugate
Quat operator-(const Quat& v);
Quat operator-(const Quat& a, const Quat& b);
Quat operator+(const Quat& a, const Quat& b);
Quat operator*(const Quat& v, float factor);
Quat operator*(float factor, const Quat& v);
Quat operator*(const Quat& a, const Quat& b);
float Len(const Quat& v);
Quat UnitVec(const Quat& v);
Quat QuatExp(const Quat& q);

// column major
struct Matrix
{
	Matrix() {}
	Matrix(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis, const Vector3& trans);
	Matrix(const Vector4& row0In, const Vector4& row1In, const Vector4& row2In, const Vector4& row3In);
	Matrix(const Quat& q);
	Matrix(const Quat& q, const Vector3& trans);
	Matrix(const Vector3& scale, const Quat& rot, const Vector3& trans);
	Matrix(const Vector3 axis, float angle);
	
	void MakeIdent();
	void MakeProjection(float fovy, float aspect, float near, float far);
	void MakeLookAt(const Vector3& eye, const Vector3& target, const Vector3& up);
	
	Vector3 GetXAxis() const;
	Vector3 GetYAxis() const;
	Vector3 GetZAxis() const;
	Vector3 GetTrans() const;
	
	void SetXAxis(const Vector3 xAxis);
	void SetYAxis(const Vector3 yAxis);
	void SetZAxis(const Vector3 zAxis);
	void SetTrans(const Vector3 trans);
	
	void SetScale(float scale);
	void SetScale(const Vector3 scale);
	
	float GetElem(int r, int c) const;
	
	Vector4 GetCol(int c) const;
	
	Vector4 row0;
	Vector4 row1;
	Vector4 row2;
	Vector4 row3;
};

Matrix operator*(const Matrix& a, const Matrix& b);
Vector3 Mul3x3(const Matrix& m, const Vector3& v);
Vector3 Mul3x4(const Matrix& m, const Vector3& v);
Vector4 Mul4x4(const Matrix& m, const Vector4& v);
Matrix Transpose(const Matrix& m);
Matrix OrthonormalInverse(const Matrix& m);

#endif
