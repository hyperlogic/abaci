/*
The MIT License

Copyright (c) 2008 Anthony J. Thibault

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#ifndef ABACI_H
#define ABACI_H

#include <math.h>

#ifdef ABACI_NAMESPACE
namespace ABACI_NAMESPACE
{
#endif

#define PI 3.14159265f

struct Matrix;

float DegToRad(float deg);
float RadToDeg(float rad);
float Clamp(float value, float min, float max);
// Limits angle between -PI & PI using modulus arithmetic
float LimitPi(float theta);
// Limits angle between zero & PI using modulus arithmetic
float Mod2Pi(float theta);
bool FuzzyEqual(float rhs, float lhs, float epsilon = 0.0001f);

struct VectorBase
{
	VectorBase() {}

	const float& operator[](int i) const;
	float& operator[](int i);

	const float& X() const;
	const float& Y() const;
	const float& Z() const;
	const float& W() const;

	float& X();
	float& Y();
	float& Z();
	float& W();

	float data[4];
};

float Dot2(const VectorBase& a, const VectorBase& b);
float Dot3(const VectorBase& a, const VectorBase& b);
float Dot4(const VectorBase& a, const VectorBase& b);

struct Vector2 : public VectorBase
{
	Vector2() {}
	Vector2(float xIn, float yIn);
	void Set(float xIn, float yIn);
	void SetZero();
};

Vector2 operator-(const Vector2& v);
Vector2 operator-(const Vector2& a, const Vector2& b);
Vector2 operator+(const Vector2& a, const Vector2& b);
Vector2 operator*(const Vector2& v, float factor);
Vector2 operator*(float factor, const Vector2& v);
float operator*(const Vector2& a, const Vector2& b);  // dot product
Vector2 CompMul(const Vector2& a, const Vector2& b);
Vector2 operator/(const Vector2& v, float denominator);
Vector2 operator/(float numerator, const Vector2& v);
Vector2 CompDiv(const Vector2& a, const Vector2& b);
float Len(const Vector2& v);
Vector2 UnitVec(const Vector2& v);

struct Vector3 : public VectorBase
{
	Vector3() {}
	Vector3(const Vector2& v, float zIn);
	Vector3(float xIn, float yIn, float zIn);
	void Set(float xIn, float yIn, float zIn);
	void SetZero();
};

Vector3 operator-(const Vector3& v);
Vector3 operator-(const Vector3& a, const Vector3& b);
Vector3 operator+(const Vector3& a, const Vector3& b);
Vector3 operator*(const Vector3& v, float factor);
Vector3 operator*(float factor, const Vector3& v);
float operator*(const Vector3& a, const Vector3& b); // dot product
Vector3 CompMul(const Vector3& a, const Vector3& b);
Vector3 operator/(const Vector3& v, float denominator);
Vector3 operator/(float numerator, const Vector3& v);
Vector3 CompDiv(const Vector3& a, const Vector3& b);
float Len(const Vector3& v);
Vector3 UnitVec(const Vector3& v);
Vector3 operator%(const Vector3& a, const Vector3& b);	// cross product

struct Vector4 : public VectorBase
{
	Vector4() {}
	Vector4(const Vector3& v, float wIn);
	Vector4(float xIn, float yIn, float zIn, float wIn);
	void Set(float xIn, float yIn, float zIn, float wIn);
	void SetZero();
};

Vector4 operator-(const Vector4& v);
Vector4 operator-(const Vector4& a, const Vector4& b);
Vector4 operator+(const Vector4& a, const Vector4& b);
Vector4 operator*(const Vector4& v, float factor);
Vector4 operator*(float factor, const Vector4& v);
float operator*(const Vector4& a, const Vector4& b); // dot product
Vector4 CompMul(const Vector4& a, const Vector4& b);
Vector4 operator/(const Vector4& v, float denominator);
Vector4 operator/(float numerator, const Vector4& v);
Vector4 CompDiv(const Vector4& a, const Vector4& b);
float Len(const Vector4& v);
Vector4 UnitVec(const Vector4& v);

struct Quat : public Vector4
{
	Quat() {}
	Quat(float xIn, float yIn, float zIn, float wIn);
	Quat(const Vector3& axis, float angle);
};

Vector3 Rotate(const Quat& q, const Vector3& v);
Quat operator~(const Quat& v);	// conjugate
Quat operator-(const Quat& v);
Quat operator-(const Quat& a, const Quat& b);
Quat operator+(const Quat& a, const Quat& b);
Quat operator*(const Quat& v, float factor);
Quat operator*(float factor, const Quat& v);
Quat operator*(const Quat& a, const Quat& b);  // quat multiply
Quat CompMul(const Quat& a, const Quat& b);
Quat operator/(const Quat& v, float denominator);
Quat operator/(float numerator, const Quat& v);
Quat CompDiv(const Quat& a, const Quat& b);
float Len(const Quat& v);
Quat UnitVec(const Quat& v);
Quat QuatExp(const Quat& q);
Quat QuatLog(const Quat& q);

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
	void MakeFrustum(float fovy, float aspect, float nearVal, float farVal);
	void MakeOrtho(float left, float right, float bottom, float top, float nearVal, float farVal);
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
	
	const float& Elem(int r, int c) const;
	float& Elem(int r, int c);
	
	Vector4 GetCol(int c) const;

	Quat GetQuat() const;
	
	Vector4 row0;
	Vector4 row1;
	Vector4 row2;
	Vector4 row3;
};

Matrix operator*(const Matrix& a, const Matrix& b);
Matrix operator+(const Matrix& a, const Matrix& b);
Matrix operator-(const Matrix& a, const Matrix& b);
Vector3 Transform3x3(const Matrix& m, const Vector3& v);
Vector3 Transform3x4(const Matrix& m, const Vector3& v);
Vector4 Transform4x4(const Matrix& m, const Vector4& v);
Matrix Transpose(const Matrix& m);
bool Inverse(const Matrix& m, Matrix& result);
Matrix OrthonormalInverse(const Matrix& m);
void PrintMatrix(const Matrix& m);

struct Complex
{
	Complex() {}
	Complex(float realIn, float imagIn);
	Complex conj() const;
	float len() const;

	float real;
	float imag;
};

Complex operator+(const Complex& a, const Complex& b);
Complex operator-(const Complex& a, const Complex& b);
Complex operator-(const Complex& a);
Complex operator+(const Complex& a);
Complex operator*(const Complex& a, const Complex& b);
Complex operator/(const Complex& a, const Complex& b);
Complex operator*(float f, const Complex& c);
Complex operator*(Complex& c, float f);
Complex expi(float imag);
float dot(const Complex& a, const Complex& b);

// inlines
#include "abaciinlines.cpp"

#ifdef ABACI_NAMESPACE
} // namespace
#endif

#endif
