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

#define PI M_PI

// Convert from degrees to radians
template <typename Scalar>
Scalar DegToRad(Scalar deg);

// Convert from radians to degrees.
template <typename Scalar>
Scalar RadToDeg(Scalar rad);

// Clamps value between min & max.
template <typename Scalar>
Scalar Clamp(Scalar value, Scalar min, Scalar max);

// Limits angle between -PI & PI using modulus arithmetic
template <typename Scalar>
Scalar LimitPi(Scalar theta);

// Limits angle between zero & 2 PI using modulus arithmetic
template <typename Scalar>
Scalar Mod2Pi(Scalar theta);

// Fuzzy comparison between two Scalar values.
template <typename Scalar>
bool FuzzyEqual(Scalar rhs, Scalar lhs, Scalar epsilon = 0.0001);

template <typename Scalar>
Scalar Lerp(Scalar a, Scalar b, Scalar t);

// forward declare
template <typename Scalar> class Complex;

//////////////////////////////////////////////////////

template <typename Scalar>
struct Vector2
{
	// Uninitialized by default.
	Vector2() {}

	// Construct from two Scalars
	Vector2(Scalar xIn, Scalar yIn);

	// Construct from complex
	Vector2(const Complex<Scalar>& complexIn);

	// Set from two Scalars
	void Set(Scalar xIn, Scalar yIn);

	// Set from complex
	void Set(const Complex<Scalar>& complexIn);

	// Sets all elements to zero.
	void SetZero();

	// const array accessor
	Scalar operator[](int i) const;

	// array accessor
	Scalar& operator[](int i);

	// Returns a vector with same direction but unit length.
	Vector2 Unit() const;

	// Returns vector length.
	Scalar Len() const;

	// Returns length squared.
	Scalar LenSq() const;

	Scalar x;
	Scalar y;
};

typedef Vector2<float> Vector2f;
typedef Vector2<double> Vector2d;

// Dot product of two vectors.
template <typename Scalar>
Scalar Dot(const Vector2<Scalar>& a, const Vector2<Scalar>& b);

// Linear interpolation between two vectors
template <typename Scalar>
Vector2<Scalar> Lerp(const Vector2<Scalar>& a, const Vector2<Scalar>& b, Scalar t);

// Unary minus
template <typename Scalar>
Vector2<Scalar> operator-(const Vector2<Scalar>& a);

// Vector subtraction.
template <typename Scalar>
Vector2<Scalar> operator-(const Vector2<Scalar>& a, const Vector2<Scalar>& b);

// Vector addition.
template <typename Scalar>
Vector2<Scalar> operator+(const Vector2<Scalar>& a, const Vector2<Scalar>& b);

// Vector decrement.
template <typename Scalar>
Vector2<Scalar>& operator-=(Vector2<Scalar>& a, const Vector2<Scalar>& b);

// Vector increment.
template <typename Scalar>
Vector2<Scalar>& operator+=(Vector2<Scalar>& a, const Vector2<Scalar>& b);

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
Vector2<Scalar> operator*(const Vector2<Scalar>& v, Scalar scalar);

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
Vector2<Scalar> operator*(Scalar factor, const Vector2<Scalar>& v);

// Vector multiplication
template <typename Scalar>
Vector2<Scalar> operator*(const Vector2<Scalar>& a, const Vector2<Scalar>& b);

// Divides all elements of a vector by a scalar.
template <typename Scalar>
Vector2<Scalar> operator/(const Vector2<Scalar>& v, Scalar denominator);

// Multiplies a scalar to the reciprical of all elements in a vector.
template <typename Scalar>
Vector2<Scalar> operator/(Scalar numerator, const Vector2<Scalar>& v);

// Vector division.
template <typename Scalar>
Vector2<Scalar> operator/(const Vector2<Scalar>& a, const Vector2<Scalar>& b);


//////////////////////////////////////////////////////

template <typename Scalar>
struct Vector3
{
	// Uninitialized by default.
	Vector3() {}

	// Construct from three Scalars
	Vector3(Scalar xIn, Scalar yIn, Scalar zIn);

	// Set from three Scalars
	void Set(Scalar xIn, Scalar yIn, Scalar zIn);

	// Set all elements to zero
	void SetZero();

	// const array accessor
	Scalar operator[](int i) const;

	// array accessor
	Scalar& operator[](int i);

	// Returns a vector with same direction but unit length.
	Vector3 Unit() const;

	// Returns vector length.
	Scalar Len() const;

	// Returns length squared.
	Scalar LenSq() const;

	Scalar x;
	Scalar y;
	Scalar z;
};

typedef Vector3<float> Vector3f;
typedef Vector3<double> Vector3d;

// Dot product of two vectors.
template <typename Scalar>
Scalar Dot(const Vector3<Scalar>& a, const Vector3<Scalar>& b);

// Cross product of two vectors.
template <typename Scalar>
Vector3<Scalar> Cross(const Vector3<Scalar>& a, const Vector3<Scalar>& b);

// Linear interpolation between two vectors
template <typename Scalar>
Vector3<Scalar> Lerp(const Vector3<Scalar>& a, const Vector3<Scalar>& b, Scalar t);

// Unary minus.
template <typename Scalar>
Vector3<Scalar> operator-(const Vector3<Scalar>& a);

// Vector subtraction.
template <typename Scalar>
Vector3<Scalar> operator-(const Vector3<Scalar>& a, const Vector3<Scalar>& b);

// Vector addition.
template <typename Scalar>
Vector3<Scalar> operator+(const Vector3<Scalar>& a, const Vector3<Scalar>& b);

// Vector addition.
template <typename Scalar>
Vector3<Scalar> operator+(const Vector3<Scalar>& a, const Vector3<Scalar>& b);

// Vector decrement.
template <typename Scalar>
Vector3<Scalar>& operator-=(Vector3<Scalar>& a, const Vector3<Scalar>& b);

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
Vector3<Scalar> operator*(const Vector3<Scalar>& v, Scalar factor);

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
Vector3<Scalar> operator*(Scalar factor, const Vector3<Scalar>& v);

// Vector multiplication
template <typename Scalar>
Vector3<Scalar> operator*(const Vector3<Scalar>& a, const Vector3<Scalar>& b);

// Divides all elements of a vector by a scalar.
template <typename Scalar>
Vector3<Scalar> operator/(const Vector3<Scalar>& v, Scalar denominator);

// Multiplies a scalar to the reciprical of all elements in a vector.
template <typename Scalar>
Vector3<Scalar> operator/(Scalar numerator, const Vector3<Scalar>& v);

// Vector division.
template <typename Scalar>
Vector3<Scalar> operator/(const Vector3<Scalar>& a, const Vector3<Scalar>& b);


//////////////////////////////////////////////////////

template <typename Scalar>
struct Vector4
{
	// Uninitialized by default.
	Vector4() {}

	// Construct from four Scalars.
	Vector4(Scalar xIn, Scalar yIn, Scalar zIn, Scalar wIn);

	// Set from four Scalars.
	void Set(Scalar xIn, Scalar yIn, Scalar zIn, Scalar wIn);

	// Set all elements to zero.
	void SetZero();

	// Returns a vector with same direction but unit length.
	Vector4 Unit() const;

	// Returns vector length.
	Scalar Len() const;

	// Returns length squared.
	Scalar LenSq() const;

	// const array accessor
	Scalar operator[](int i) const;

	// array accessor
	Scalar& operator[](int i);

	Scalar x;
	Scalar y;
	Scalar z;
	Scalar w;
};

typedef Vector4<float> Vector4f;
typedef Vector4<double> Vector4d;

// Dot product of two vectors.
template <typename Scalar>
Scalar Dot(const Vector4<Scalar>& a, const Vector4<Scalar>& b);

// Linear interpolation between two vectors
template <typename Scalar>
Vector4<Scalar> Lerp(const Vector4<Scalar>& a, const Vector4<Scalar>& b, Scalar t);

// Unary minus.
template <typename Scalar>
Vector4<Scalar> operator-(const Vector4<Scalar>& v);

// Vector subtraction.
template <typename Scalar>
Vector4<Scalar> operator-(const Vector4<Scalar>& a, const Vector4<Scalar>& b);

// Vector addition.
template <typename Scalar>
Vector4<Scalar> operator+(const Vector4<Scalar>& a, const Vector4<Scalar>& b);

// Vector addition.
template <typename Scalar>
Vector4<Scalar> operator+(const Vector4<Scalar>& a, const Vector4<Scalar>& b);

// Vector decrement.
template <typename Scalar>
Vector4<Scalar>& operator-=(Vector4<Scalar>& a, const Vector4<Scalar>& b);

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
Vector4<Scalar> operator*(const Vector4<Scalar>& v, Scalar factor);

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
Vector4<Scalar> operator*(Scalar factor, const Vector4<Scalar>& v);

// Vector multiplication.
template <typename Scalar>
Vector4<Scalar> operator*(const Vector4<Scalar>& a, const Vector4<Scalar>& b);

// Divides all elements of a vector by a scalar.
template <typename Scalar>
Vector4<Scalar> operator/(const Vector4<Scalar>& v, Scalar denominator);

// Multiplies a scalar to the reciprical of all elements in a vector.
template <typename Scalar>
Vector4<Scalar> operator/(Scalar numerator, const Vector4<Scalar>& v);

// Vector division.
template <typename Scalar>
Vector4<Scalar> operator/(const Vector4<Scalar>& a, const Vector4<Scalar>& b);


//////////////////////////////////////////////////////

template <typename Scalar>
struct Quat
{
	// Create from axis and angle
	static Quat AxisAngle(const Vector3<Scalar>& axis, Scalar angle);

	// Uninitialized by default.
	Quat() {}

	// Construct from four Scalars.
	Quat(Scalar iIn, Scalar jIn, Scalar kIn, Scalar rIn);

	// Set from four Scalars.
	void Set(Scalar iIn, Scalar jIn, Scalar kIn, Scalar rIn);

	// Set all elements to zero.
	void SetZero();

	// Returns this quaternion normalized.
	Quat Unit() const;

	// Returns quat length.
	Scalar Len() const;

	// Returns quat length squared.
	Scalar LenSq() const;

	// Rotate a vector.
	Vector3<Scalar> Rotate(const Vector3<Scalar>& v) const;

	// Quaternian exponential
	Quat Exp() const;

	// Quaternian logarithm
	Quat Log() const;

	Scalar i;
	Scalar j;
	Scalar k;
	Scalar r;  // real part
};

typedef Quat<float> Quatf;
typedef Quat<double> Quatd;

// Dot product
template <typename Scalar>
Scalar Dot(const Quat<Scalar>& a, const Quat<Scalar>& b);

// Quaternion conjugate
template <typename Scalar>
Quat<Scalar> operator~(const Quat<Scalar>& v);

// Unary minus.
template <typename Scalar>
Quat<Scalar> operator-(const Quat<Scalar>& v);

// Quaternion subtraction.
template <typename Scalar>
Quat<Scalar> operator-(const Quat<Scalar>& a, const Quat<Scalar>& b);

// Quaternion addition.
template <typename Scalar>
Quat<Scalar> operator+(const Quat<Scalar>& a, const Quat<Scalar>& b);

// Quaternion multplication.
template <typename Scalar>
Quat<Scalar> operator*(const Quat<Scalar>& a, const Quat<Scalar>& b);


//////////////////////////////////////////////////////

// column major
template <typename Scalar>
struct Matrix
{
	// Create a Matrix from three principle axes and a translation.
	static Matrix Axes(const Vector3<Scalar>& xAxis, const Vector3<Scalar>& yAxis, const Vector3<Scalar>& zAxis, const Vector3<Scalar>& trans = Vector3<Scalar>(0,0,0));

	// Create a Matrix from four row vectors.
	static Matrix Rows(const Vector4<Scalar>& row0In, const Vector4<Scalar>& row1In, const Vector4<Scalar>& row2In, const Vector4<Scalar>& row3In);

	// Create a Matrix from a Quat and a translation.
	static Matrix QuatTrans(const Quat<Scalar>& q, const Vector3<Scalar>& trans);

	// Create a Matrix from a Scale vector, a Quat and a translation.
	static Matrix ScaleQuatTrans(const Vector3<Scalar>& scale, const Quat<Scalar>& rot, const Vector3<Scalar>& trans);

	// Create a Matrix from a rotation represented by an axis and an angle.
	static Matrix AxisAngle(const Vector3<Scalar>& axis, Scalar angle);

	// Create an identity Matrix.
	static Matrix Identity();

	// Create a persective projection Matrix.
	static Matrix Frustum(Scalar fovy, Scalar aspect, Scalar nearVal, Scalar farVal);

	// Create an orthograpic projection Matrix.
	static Matrix Ortho(Scalar left, Scalar right, Scalar bottom, Scalar top, Scalar nearVal, Scalar farVal);

	// Create a look at matrix. (x is forward)
	static Matrix LookAt(const Vector3<Scalar>& eye, const Vector3<Scalar>& target, const Vector3<Scalar>& up);
	
	// Uninitialized by default.
	Matrix() {}
	
	// Axes accessors
	Vector3<Scalar> GetXAxis() const;
	Vector3<Scalar> GetYAxis() const;
	Vector3<Scalar> GetZAxis() const;
	Vector3<Scalar> GetTrans() const;
	
	void SetXAxis(const Vector3<Scalar>& xAxis);
	void SetYAxis(const Vector3<Scalar>& yAxis);
	void SetZAxis(const Vector3<Scalar>& zAxis);
	void SetTrans(const Vector3<Scalar>& trans);
	
	// Multiplies by uniform scale.
	void SetScale(Scalar scale);

	// Multiplies by non-uniform scale.
	void SetScale(const Vector3<Scalar>& scale);
	
	// Element accessors
	const Scalar& Elem(int r, int c) const;
	Scalar& Elem(int r, int c);
	
	// Column accessor
	Vector4<Scalar> GetCol(int c) const;

	// Returns the rotation component of this Matrix.
	Quat<Scalar> GetQuat() const;

	// Multiply the 3x3 component of this Matrix with a column vector.
	Vector3<Scalar> Mul3x3(const Vector3<Scalar>& v) const;

	// Multiply the 3x4 component of this Matrix with a column vector. (w component of vector is 1.0)
	Vector3<Scalar> Mul3x4(const Vector3<Scalar>& v) const;

	// Multiply this Matrix with a column vector.
	Vector4<Scalar> Mul4x4(const Vector4<Scalar>& v) const;

	// Returns the transpose of this Matrix
	Matrix Transpose() const;

	// If the 3x3 portion of this Matrix is Orthogonal (i.e. columns are orthogonal unit vectors)
	// this will return the Inverse of that matrix.
	Matrix OrthoInverse() const;

	// Full 4x4 Matrix Inverse, returns Identity if matrix has no inverse.
	Matrix FullInverse() const;

	Vector4<Scalar> row0;
	Vector4<Scalar> row1;
	Vector4<Scalar> row2;
	Vector4<Scalar> row3;
};

typedef Matrix<float> Matrixf;
typedef Matrix<double> Matrixd;

// Matrix addition
template <typename Scalar>
Matrix<Scalar> operator+(const Matrix<Scalar>& a, const Matrix<Scalar>& b);

// Matrix subtraction
template <typename Scalar>
Matrix<Scalar> operator-(const Matrix<Scalar>& a, const Matrix<Scalar>& b);

// Matrix multiplication
template <typename Scalar>
Matrix<Scalar> operator*(const Matrix<Scalar>& a, const Matrix<Scalar>& b);

// Full 4x4 Matrix inverse, returns false if Matrix has no inverse.
template <typename Scalar>
bool FullInverse(const Matrix<Scalar>& m, Matrix<Scalar>& result);

// Print to stdout.
template <typename Scalar>
void PrintMatrix(const Matrix<Scalar>& m);


//////////////////////////////////////////////////////

template <typename Scalar>
struct Complex
{
	// Uninitialized by default.
	Complex() {}

	// Construct from two floats
	Complex(Scalar rIn, Scalar iIn);

	// Construct from a Vector2
	Complex(const Vector2<Scalar>& vector2In);

	// Length
	Scalar Len() const;

	// Square of length
	Scalar LenSq() const;

	Scalar r;
	Scalar i;
};

typedef Complex<float> Complexf;
typedef Complex<double> Complexd;

// Dot product
template <typename Scalar>
Scalar Dot(const Complex<Scalar>& a, const Complex<Scalar>& b);

// Complex conjugate
template <typename Scalar>
Complex<Scalar> operator~(const Complex<Scalar>& a);

// Unary minus.
template <typename Scalar>
Complex<Scalar> operator-(const Complex<Scalar>& a);

// Unary plus.
template <typename Scalar>
Complex<Scalar> operator+(const Complex<Scalar>& a);

// Complex addition.
template <typename Scalar>
Complex<Scalar> operator+(const Complex<Scalar>& a, const Complex<Scalar>& b);

// Complex subtraction.
template <typename Scalar>
Complex<Scalar> operator-(const Complex<Scalar>& a, const Complex<Scalar>& b);

// Complex multiplication.
template <typename Scalar>
Complex<Scalar> operator*(const Complex<Scalar>& a, const Complex<Scalar>& b);

// Complex division.
template <typename Scalar>
Complex<Scalar> operator/(const Complex<Scalar>& a, const Complex<Scalar>& b);

// Multiplication by a real number.
template <typename Scalar>
Complex<Scalar> operator*(Scalar scalar, const Complex<Scalar>& c);

// Multiplication by a real number.
template <typename Scalar>
Complex<Scalar> operator*(const Complex<Scalar>& c, Scalar scalar);

// e ^ 0 + xi
template <typename Scalar>
Complex<Scalar> ExpI(Scalar x);

// Square root
template <typename Scalar>
Complex<Scalar> sqrt(const Complex<Scalar>& z);

// Exponent
template <typename Scalar>
Complex<Scalar> exp(const Complex<Scalar>& z);

// Natural Logarithm
template <typename Scalar>
Complex<Scalar> log(const Complex<Scalar>& z);


// inlines
#include "abaciinlines.cpp"

#ifdef ABACI_NAMESPACE
} // namespace
#endif

#endif
