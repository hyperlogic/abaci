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
struct Vector3;
struct Vector4;

// Convert from degrees to radians
float DegToRad(float deg);

// Convert from radians to degrees.
float RadToDeg(float rad);

// Clamps value between min & max.
float Clamp(float value, float min, float max);

// Limits angle between -PI & PI using modulus arithmetic
float LimitPi(float theta);

// Limits angle between zero & PI using modulus arithmetic
float Mod2Pi(float theta);

// Fuzzy comparison between two float values.
bool FuzzyEqual(float rhs, float lhs, float epsilon = 0.0001f);


//////////////////////////////////////////////////////
struct Vector2
{
	// Uninitialized by default.
	Vector2() {}

	// Construct from two floats
	Vector2(float xIn, float yIn);

	// Set from two floats
	void Set(float xIn, float yIn);

	// Sets all elements to zero.
	void SetZero();

	// Returns a vector with same direction but unit length.
	Vector2 Unit() const;

	// Returns vector length.
	float Len() const;

	// Returns length squared.
	float LenSq() const;

	float x;
	float y;
};

// Dot product of two vectors.
float Dot(const Vector2& a, const Vector2& b);

// Linear interpolation between two vectors
Vector2 Lerp(const Vector2& a, const Vector2& b, float t);

// Unary minus
Vector2 operator-(const Vector2& a);

// Vector subtraction.
Vector2 operator-(const Vector2& a, const Vector2& b);

// Vector addition.
Vector2 operator+(const Vector2& a, const Vector2& b);

// Multplies all elements of a vector by a scalar.
Vector2 operator*(const Vector2& v, float scalar);

// Multplies all elements of a vector by a scalar.
Vector2 operator*(float factor, const Vector2& v);

// Vector multiplication
Vector2 operator*(const Vector2& a, const Vector2& b);

// Divides all elements of a vector by a scalar.
Vector2 operator/(const Vector2& v, float denominator);

// Multiplies a scalar to the reciprical of all elements in a vector.
Vector2 operator/(float numerator, const Vector2& v);

// Vector division.
Vector2 operator/(const Vector2& a, const Vector2& b);


//////////////////////////////////////////////////////

struct Vector3
{
	// Uninitialized by default.
	Vector3() {}

	// Construct from three floats
	Vector3(float xIn, float yIn, float zIn);

	// Set from three floats
	void Set(float xIn, float yIn, float zIn);

	// Set all elements to zero
	void SetZero();

	// Returns a vector with same direction but unit length.
	Vector3 Unit() const;

	// Returns vector length.
	float Len() const;

	// Returns length squared.
	float LenSq() const;

	float x;
	float y;
	float z;
};

// Dot product of two vectors.
float Dot(const Vector3& a, const Vector3& b);

// Cross product of two vectors.
Vector3 Cross(const Vector3& a, const Vector3& b);

// Linear interpolation between two vectors
Vector3 Lerp(const Vector3& a, const Vector3& a, float t);

// Unary minus.
Vector3 operator-(const Vector3& a);

// Vector subtraction.
Vector3 operator-(const Vector3& a, const Vector3& b);

// Vector addition.
Vector3 operator+(const Vector3& a, const Vector3& b);

// Multplies all elements of a vector by a scalar.
Vector3 operator*(const Vector3& v, float factor);

// Multplies all elements of a vector by a scalar.
Vector3 operator*(float factor, const Vector3& v);

// Vector multiplication
Vector3 operator*(const Vector3& a, const Vector3& b);

// Divides all elements of a vector by a scalar.
Vector3 operator/(const Vector3& v, float denominator);

// Multiplies a scalar to the reciprical of all elements in a vector.
Vector3 operator/(float numerator, const Vector3& v);

// Vector division.
Vector3 operator/(const Vector3& a, const Vector3& b);


//////////////////////////////////////////////////////

struct Vector4
{
	// Uninitialized by default.
	Vector4() {}

	// Construct from four floats.
	Vector4(float xIn, float yIn, float zIn, float wIn);

	// Set from four floats.
	void Set(float xIn, float yIn, float zIn, float wIn);

	// Set all elements to zero.
	void SetZero();

	// Returns a vector with same direction but unit length.
	Vector4 Unit() const;

	// Returns vector length.
	float Len() const;

	// Returns length squared.
	float LenSq() const;

	// const array accessor
	float operator[](int i) const;

	// array accessor
	float& operator[](int i);

	float x;
	float y;
	float z;
	float w;
};

// Dot product of two vectors.
float Dot(const Vector4& a, const Vector4& b);

// Linear interpolation between two vectors
Vector4 Lerp(const Vector4& a, const Vector4& b, float t);

// Unary minus.
Vector4 operator-(const Vector4& v);

// Vector subtraction.
Vector4 operator-(const Vector4& a, const Vector4& b);

// Vector addition.
Vector4 operator+(const Vector4& a, const Vector4& b);

// Multplies all elements of a vector by a scalar.
Vector4 operator*(const Vector4& v, float factor);

// Multplies all elements of a vector by a scalar.
Vector4 operator*(float factor, const Vector4& v);

// Vector multiplication.
Vector4 operator*(const Vector4& a, const Vector4& b);

// Divides all elements of a vector by a scalar.
Vector4 operator/(const Vector4& v, float denominator);

// Multiplies a scalar to the reciprical of all elements in a vector.
Vector4 operator/(float numerator, const Vector4& v);

// Vector division.
Vector4 operator/(const Vector4& a, const Vector4& b);


//////////////////////////////////////////////////////

struct Quat
{
	// Create from axis and angle
	static Quat AxisAngle(const Vector3& axis, float angle);

	// Uninitialized by default.
	Quat() {}

	// Construct from four floats.
	Quat(float iIn, float jIn, float kIn, float rIn);

	// Set all elements to zero.
	void SetZero();

	// Returns this quaternion normalized.
	Quat Unit() const;

	// Returns quat length.
	float Len() const;

	// Returns quat length squared.
	float LenSq() const;

	// Rotate a vector.
	Vector3 Rotate(const Vector3& v) const;

	// Quaternian exponential
	Quat Exp() const;

	// Quaternian logarithm
	Quat Log() const;

	float i;
	float j;
	float k;
	float r;  // real part
};

// Dot product
float Dot(const Quat& a, const Quat& b);

// Quaternion conjugate
Quat operator~(const Quat& v);

// Unary minus.
Quat operator-(const Quat& v);

// Quaternion subtraction.
Quat operator-(const Quat& a, const Quat& b);

// Quaternion addition.
Quat operator+(const Quat& a, const Quat& b);

// Quaternion multplication.
Quat operator*(const Quat& a, const Quat& b);


//////////////////////////////////////////////////////

// column major
struct Matrix
{
	// Create a Matrix from three principle axes and a translation.
	static Matrix Axes(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis, const Vector3& trans = Vector3(0,0,0));

	// Create a Matrix from four row vectors.
	static Matrix Rows(const Vector4& row0In, const Vector4& row1In, const Vector4& row2In, const Vector4& row3In);

	// Create a Matrix from a Quat and a translation.
	static Matrix QuatTrans(const Quat& q, const Vector3& trans);

	// Create a Matrix from a Scale vector, a Quat and a translation.
	static Matrix ScaleQuatTrans(const Vector3& scale, const Quat& rot, const Vector3& trans);

	// Create a Matrix from a rotation represented by an axis and an angle.
	static Matrix AxisAngle(const Vector3 axis, float angle);

	// Create an identity Matrix.
	static Matrix Identity();

	// Create a persective projection Matrix.
	static Matrix Frustum(float fovy, float aspect, float nearVal, float farVal);

	// Create an orthograpic projection Matrix.
	static Matrix Ortho(float left, float right, float bottom, float top, float nearVal, float farVal);

	// Create a look at matrix. (x is forward)
	static Matrix LookAt(const Vector3& eye, const Vector3& target, const Vector3& up);
	
	// Uninitialized by default.
	Matrix() {}
	
	// Axes accessors
	Vector3 GetXAxis() const;
	Vector3 GetYAxis() const;
	Vector3 GetZAxis() const;
	Vector3 GetTrans() const;
	
	void SetXAxis(const Vector3 xAxis);
	void SetYAxis(const Vector3 yAxis);
	void SetZAxis(const Vector3 zAxis);
	void SetTrans(const Vector3 trans);
	
	// Multiplies by uniform scale.
	void SetScale(float scale);

	// Multiplies by non-uniform scale.
	void SetScale(const Vector3 scale);
	
	// Element accessors
	const float& Elem(int r, int c) const;
	float& Elem(int r, int c);
	
	// Column accessor
	Vector4 GetCol(int c) const;

	// Returns the rotation component of this Matrix.
	Quat GetQuat() const;

	// Multiply the 3x3 component of this Matrix with a column vector.
	Vector3 Mul3x3(const Vector3& v) const;

	// Multiply the 3x4 component of this Matrix with a column vector. (w component of vector is 1.0)
	Vector3 Mul3x4(const Vector3& v) const;

	// Multiply this Matrix with a column vector.
	Vector4 Mul4x4(const Vector4& v) const;

	// Returns the transpose of this Matrix
	Matrix Transpose() const;

	// If the 3x3 portion of this Matrix is Orthogonal (i.e. columns are orthogonal unit vectors)
	// this will return the Inverse of that matrix.
	Matrix OrthoInverse() const;

	// Full 4x4 Matrix Inverse, returns Identity if matrix has no inverse.
	Matrix FullInverse() const;

	Vector4 row0;
	Vector4 row1;
	Vector4 row2;
	Vector4 row3;
};

// Matrix addition
Matrix operator+(const Matrix& a, const Matrix& b);

// Matrix subtraction
Matrix operator-(const Matrix& a, const Matrix& b);

// Matrix multiplication
Matrix operator*(const Matrix& a, const Matrix& b);

// Full 4x4 Matrix inverse, returns false if Matrix has no inverse.
bool FullInverse(const Matrix& m, Matrix& result);

// Print to stdout.
void PrintMatrix(const Matrix& m);


//////////////////////////////////////////////////////

struct Complex
{
	// Uninitialized by default.
	Complex() {}

	// Construct from two floats
	Complex(float rIn, float iIn);

	// Length
	float Len() const;

	// Square of length
	float LenSq() const;

	float r;
	float i;
};

// Dot product
float Dot(const Complex& a, const Complex& b);

// Complex conjugate
Complex operator~(const Complex& a);

// Unary minus.
Complex operator-(const Complex& a);

// Unary plus.
Complex operator+(const Complex& a);

// Complex addition.
Complex operator+(const Complex& a, const Complex& b);

// Complex subtraction.
Complex operator-(const Complex& a, const Complex& b);

// Complex multiplication.
Complex operator*(const Complex& a, const Complex& b);

// Complex division.
Complex operator/(const Complex& a, const Complex& b);

// Multiplication by a real number.
Complex operator*(float f, const Complex& c);

// Multiplication by a real number.
Complex operator*(Complex& c, float f);

// e ^ 0 + xi
Complex ExpI(float x);

// inlines
#include "abaciinlines.cpp"

#ifdef ABACI_NAMESPACE
} // namespace
#endif

#endif
