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

#ifdef ABACI_H

#include <stdio.h>

// Convert from degrees to radians
template <typename Scalar>
inline Scalar DegToRad(Scalar deg)
{ 
	return deg * (PI / 180.0);
}

// Convert from radians to degrees.
template <typename Scalar>
inline Scalar RadToDeg(Scalar rad)
{
	return rad * (180.0 / PI);
}

// Clamps value between min & max.
template <typename Scalar>
inline Scalar Clamp(Scalar value, Scalar min, Scalar max)
{
	Scalar result = value;
	if (value > max)
		result = max;
	else if (value < min)
		result = min;
	return result; 
}

// Limits angle between -PI & PI using modulus arithmetic
template <typename Scalar>
inline Scalar LimitPi(Scalar theta)
{
	return fmod(theta + PI, 2 * PI) - PI;
}

// Limits angle between zero & PI using modulus arithmetic
template <typename Scalar>
inline Scalar Mod2Pi(Scalar theta)
{
	return fmod(theta, 2 * PI);
}

// Fuzzy comparison between two Scalar values.
template <typename Scalar>
inline bool FuzzyEqual(Scalar rhs, Scalar lhs, Scalar epsilon)
{
	return fabs(rhs - lhs) <= epsilon;
}

////////////////////////////////////////////////////////////////////////////////
// Vector2

// Construct from two Scalars
template <typename Scalar>
inline Vector2<Scalar>::Vector2(Scalar xIn, Scalar yIn) : x(xIn), y(yIn) {}

// Set from two Scalars
template <typename Scalar>
inline void Vector2<Scalar>::Set(Scalar xIn, Scalar yIn)
{
	x = xIn;
	y = yIn;
}

// Sets all elements to zero.
template <typename Scalar>
inline void Vector2<Scalar>::SetZero()
{
	x = y = 0;
}

// Returns a vector with same direction but unit length.
template <typename Scalar>
inline Vector2<Scalar> Vector2<Scalar>::Unit() const
{
	return *this / Len();
}

// Returns vector length.
template <typename Scalar>
inline Scalar Vector2<Scalar>::Len() const
{
	return sqrt(Dot(*this, *this));
}

// Returns length squared.
template <typename Scalar>
inline Scalar Vector2<Scalar>::LenSq() const
{
	return Dot(*this, *this);
}

// Dot product of two vectors.
template <typename Scalar>
inline Scalar Dot(const Vector2<Scalar>& a, const Vector2<Scalar>& b)
{
	return a.x * b.x + a.y * b.y;
}

// Linear interpolation between two vectors
template <typename Scalar>
inline Vector2<Scalar> Lerp(const Vector2<Scalar>& a, const Vector2<Scalar>& b, Scalar t)
{
    return a + (b - a) * t;
}

// Unary minus
template <typename Scalar>
inline Vector2<Scalar> operator-(const Vector2<Scalar>& a)
{
	return Vector2<Scalar>(-a.x, -a.y);
}

// Vector subtraction.
template <typename Scalar>
inline Vector2<Scalar> operator-(const Vector2<Scalar>& a, const Vector2<Scalar>& b)
{
	return Vector2<Scalar>(a.x - b.x, a.y - b.y);
}

// Vector addition.
template <typename Scalar>
inline Vector2<Scalar> operator+(const Vector2<Scalar>& a, const Vector2<Scalar>& b)
{
	return Vector2<Scalar>(a.x + b.x, a.y + b.y);
}

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
inline Vector2<Scalar> operator*(const Vector2<Scalar>& v, Scalar scalar)
{
	return Vector2<Scalar>(scalar * v.x, scalar * v.y);
}

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
inline Vector2<Scalar> operator*(Scalar scalar, const Vector2<Scalar>& v)
{
	return Vector2<Scalar>(scalar * v.x, scalar * v.y);
}

// Vector multiplication.
template <typename Scalar>
inline Vector2<Scalar> operator*(const Vector2<Scalar>& a, const Vector2<Scalar>& b)
{
	return Vector2<Scalar>(a.x * b.x, a.y * b.y);
}

// Divides all elements of a vector by a scalar.
template <typename Scalar>
inline Vector2<Scalar> operator/(const Vector2<Scalar>& v, Scalar denominator)
{
	return Vector2<Scalar>(v.x / denominator, v.y / denominator);
}

// Multiplies a scalar to the reciprical of all elements in a vector.
template <typename Scalar>
inline Vector2<Scalar> operator/(Scalar numerator, const Vector2<Scalar>& v)
{
	return Vector2<Scalar>(numerator / v.x, numerator / v.y);
}

// Vector division.
template <typename Scalar>
inline Vector2<Scalar> operator/(const Vector2<Scalar>& a, const Vector2<Scalar>& b)
{
	return Vector2<Scalar>(a.x / b.x, a.y / b.y);
}

////////////////////////////////////////////////////////////////////////////////
// Vector3

// Construct from three Scalars
template <typename Scalar>
inline Vector3<Scalar>::Vector3(Scalar xIn, Scalar yIn, Scalar zIn) : x(xIn), y(yIn), z(zIn) {}

// Set from three Scalars
template <typename Scalar>
inline void Vector3<Scalar>::Set(Scalar xIn, Scalar yIn, Scalar zIn)
{
	x = xIn;
	y = yIn;
	z = zIn;
}

// Set all elements to zero
template <typename Scalar>
inline void Vector3<Scalar>::SetZero()
{
	x = y = z = 0;
}

// Returns a vector with same direction but unit length.
template <typename Scalar>
inline Vector3<Scalar> Vector3<Scalar>::Unit() const
{
	return (*this) / Len();
}

// Returns vector length.
template <typename Scalar>
inline Scalar Vector3<Scalar>::Len() const
{
	return sqrt(Dot(*this, *this));
}

// Returns vector length squared.
template <typename Scalar>
inline Scalar Vector3<Scalar>::LenSq() const
{
	return Dot(*this, *this);
}

// Dot product of two vectors.
template <typename Scalar>
inline Scalar Dot(const Vector3<Scalar>& a, const Vector3<Scalar>& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Cross product of two vectors.
template <typename Scalar>
inline Vector3<Scalar> Cross(const Vector3<Scalar>& a, const Vector3<Scalar>& b)
{
	return Vector3<Scalar>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

// Linear interpolation between two vectors
template <typename Scalar>
inline Vector3<Scalar> Lerp(const Vector3<Scalar>& a, const Vector3<Scalar>& b, Scalar t)
{
    return a + (b - a) * t;
}

// Unary minus
template <typename Scalar>
inline Vector3<Scalar> operator-(const Vector3<Scalar>& a)
{
	return Vector3<Scalar>(-a.x, -a.y, -a.z);
}

// Vector subtraction.
template <typename Scalar>
inline Vector3<Scalar> operator-(const Vector3<Scalar>& a, const Vector3<Scalar>& b)
{
	return Vector3<Scalar>(a.x - b.x, a.y - b.y, a.z - b.z);
}

// Vector addition.
template <typename Scalar>
inline Vector3<Scalar> operator+(const Vector3<Scalar>& a, const Vector3<Scalar>& b)
{
	return Vector3<Scalar>(a.x + b.x, a.y + b.y, a.z + b.z);
}

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
inline Vector3<Scalar> operator*(const Vector3<Scalar>& v, Scalar factor)
{
	return Vector3<Scalar>(v.x * factor, v.y * factor, v.z * factor);
}

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
inline Vector3<Scalar> operator*(Scalar factor, const Vector3<Scalar>& v)
{
	return Vector3<Scalar>(v.x * factor, v.y * factor, v.z * factor);
}

// Vector multiplication
template <typename Scalar>
inline Vector3<Scalar> operator*(const Vector3<Scalar>& a, const Vector3<Scalar>& b)
{
	return Vector3<Scalar>(a.x * b.x, a.y * b.y, a.z * b.z);
}

// Divides all elements of a vector by a scalar.
template <typename Scalar>
inline Vector3<Scalar> operator/(const Vector3<Scalar>& v, Scalar denominator)
{
	return Vector3<Scalar>(v.x / denominator, v.y / denominator, v.z / denominator);
}

// Multiplies a scalar to the reciprical of all elements in a vector.
template <typename Scalar>
inline Vector3<Scalar> operator/(Scalar numerator, const Vector3<Scalar>& v)
{
	return Vector3<Scalar>(numerator / v.x, numerator / v.y, numerator / v.z);
}

// Vector division.
template <typename Scalar>
inline Vector3<Scalar> operator/(const Vector3<Scalar>& a, const Vector3<Scalar>& b)
{
	return Vector3<Scalar>(a.x / b.x, a.y / b.y, a.z / b.z);
}

////////////////////////////////////////////////////////////////////////////////
// Vector4

// Construct from four Scalars.
template <typename Scalar>
inline Vector4<Scalar>::Vector4(Scalar xIn, Scalar yIn, Scalar zIn, Scalar wIn) : x(xIn), y(yIn), z(zIn), w(wIn) {}

// Set from four Scalars.
template <typename Scalar>
inline void Vector4<Scalar>::Set(Scalar xIn, Scalar yIn, Scalar zIn, Scalar wIn)
{
	x = xIn;
	y = yIn;
	z = zIn;
	w = wIn;
}

// Set all elements to zero.
template <typename Scalar>
inline void Vector4<Scalar>::SetZero()
{
	x = y = z = w = 0;
}

// Returns a vector with same direction but unit length.
template <typename Scalar>
inline Vector4<Scalar> Vector4<Scalar>::Unit() const
{
	return (*this) / Len();
}

// Returns vector length.
template <typename Scalar>
inline Scalar Vector4<Scalar>::Len() const
{
	return sqrt(Dot(*this, *this));
}

// Returns length squared.
template <typename Scalar>
inline Scalar Vector4<Scalar>::LenSq() const
{
	return Dot(*this, *this);
}

// const array accessor
template <typename Scalar>
inline Scalar Vector4<Scalar>::operator[](int i) const
{
	return *(&x + i);
}

// array accessor
template <typename Scalar>
inline Scalar& Vector4<Scalar>::operator[](int i)
{
	return *(&x + i);
}

// Dot product of two vectors.
template <typename Scalar>
inline Scalar Dot(const Vector4<Scalar>& a, const Vector4<Scalar>& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

// Linear interpolation between two vectors
template <typename Scalar>
inline Vector4<Scalar> Lerp(const Vector4<Scalar>& a, const Vector4<Scalar>& b, Scalar t)
{
	return a + (b - a) * t;
}

// Unary minus.
template <typename Scalar>
inline Vector4<Scalar> operator-(const Vector4<Scalar>& v)
{
	return Vector4<Scalar>(-v.x, -v.y, -v.z, -v.w);
}

// Vector subtraction.
template <typename Scalar>
inline Vector4<Scalar> operator-(const Vector4<Scalar>& a, const Vector4<Scalar>& b)
{
	return Vector4<Scalar>(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

// Vector addition.
template <typename Scalar>
inline Vector4<Scalar> operator+(const Vector4<Scalar>& a, const Vector4<Scalar>& b)
{
	return Vector4<Scalar>(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
inline Vector4<Scalar> operator*(const Vector4<Scalar>& v, Scalar factor)
{
	return Vector4<Scalar>(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

// Multplies all elements of a vector by a scalar.
template <typename Scalar>
inline Vector4<Scalar> operator*(Scalar factor, const Vector4<Scalar>& v)
{
	return Vector4<Scalar>(factor * v.x, factor * v.y, factor * v.z, factor * v.w);	
}

// Vector multiplication.
template <typename Scalar>
inline Vector4<Scalar> operator*(const Vector4<Scalar>& a, const Vector4<Scalar>& b)
{
	return Vector4<Scalar>(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);	
}

// Divides all elements of a vector by a scalar.
template <typename Scalar>
inline Vector4<Scalar> operator/(const Vector4<Scalar>& v, Scalar denominator)
{
	return Vector4<Scalar>(v.x / denominator, v.y / denominator, v.z / denominator, v.w / denominator);
}

// Multiplies a scalar to the reciprical of all elements in a vector.
template <typename Scalar>
inline Vector4<Scalar> operator/(Scalar numerator, const Vector4<Scalar>& v)
{
	return Vector4<Scalar>(numerator / v.x, numerator / v.y, numerator / v.z, numerator / v.w);
}

// Vector division.
template <typename Scalar>
inline Vector4<Scalar> operator/(const Vector4<Scalar>& a, const Vector4<Scalar>& b)
{
	return Vector4<Scalar>(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}

////////////////////////////////////////////////////////////////////////////////
// Quat

// Create from axis and angle
template <typename Scalar>
inline Quat<Scalar> Quat<Scalar>::AxisAngle(const Vector3<Scalar>& axis, Scalar angle)
{
	Vector3<Scalar> n = axis.Unit() * static_cast<Scalar>(sin(angle/2.0f));
	Scalar w = cos(angle/2.0f);
	return Quat<Scalar>(n.x, n.y, n.z, w);
}

// Construct from four Scalars.
template <typename Scalar>
inline Quat<Scalar>::Quat(Scalar iIn, Scalar jIn, Scalar kIn, Scalar rIn) : i(iIn), j(jIn), k(kIn), r(rIn) {}

// Set all elements to zero.
template <typename Scalar>
inline void Quat<Scalar>::SetZero()
{
	i = j = k = r = 0;
}

// Returns this quaternion normalized.
template <typename Scalar>
inline Quat<Scalar> Quat<Scalar>::Unit() const
{
	Scalar len = Len();
	return Quat<Scalar>(i / len, j / len, k / len, r / len);
}

// Returns quat length.
template <typename Scalar>
inline Scalar Quat<Scalar>::Len() const
{
	return sqrt(Dot(*this, *this));
}

// Returns quat length squared.
template <typename Scalar>
inline Scalar Quat<Scalar>::LenSq() const
{
	return Dot(*this, *this);
}

// Rotate a vector.
template <typename Scalar>
inline Vector3<Scalar> Quat<Scalar>::Rotate(const Vector3<Scalar>& v) const
{
	Quat<Scalar> r = (*this) * Quat<Scalar>(v.x, v.y, v.z, 0.0f) * ~(*this);
	return Vector3<Scalar>(r.i, r.j, r.k);
}

// Quaternian exponential
template <typename Scalar>
Quat<Scalar> Quat<Scalar>::Exp() const
{
	Scalar angle = Vector3<Scalar>(i, j, k).Len();
	Vector3<Scalar> n;
	if (angle > 0.0001f)
		n = Vector3<Scalar>(i, j, k).Unit() * static_cast<Scalar>(sin(angle / 2));
	else
		n.Set(0,0,0);

	return Quat<Scalar>(n.x, n.y, n.z, cos(angle / 2));
}

// Quaternian logarithm
template <typename Scalar>
Quat<Scalar> Quat<Scalar>::Log() const
{
	Scalar cos_a = r;
	if (cos_a > 1) cos_a = 1;
	if (cos_a < -1) cos_a = -1;

    Scalar sin_a = sqrt(1 - cos_a * cos_a);

    if (fabs(sin_a) < 0.0005)
		sin_a = 1;
	else
		sin_a = 1/sin_a;

    Scalar angle = 2 * acos(cos_a);
	Quat<Scalar> log;
    log.i = i * sin_a * angle;
    log.j = j * sin_a * angle;
    log.k = k * sin_a * angle;
	log.r = 0;
	return log;
}

// Dot product
template <typename Scalar>
inline Scalar Dot(const Quat<Scalar>& a, const Quat<Scalar>& b)
{
	return a.i * b.i + a.j * b.j + a.k * b.k + a.r * b.r;
}

// Quaternion conjugate
template <typename Scalar>
inline Quat<Scalar> operator~(const Quat<Scalar>& v)
{
	return Quat<Scalar>(-v.i, -v.j, -v.k, v.r);
}

// Unary minus.
template <typename Scalar>
inline Quat<Scalar> operator-(const Quat<Scalar>& v)
{
	return Quat<Scalar>(-v.i, -v.j, -v.k, -v.r);
}

// Quaternion subtraction.
template <typename Scalar>
inline Quat<Scalar> operator-(const Quat<Scalar>& a, const Quat<Scalar>& b)
{
	return Quat<Scalar>(a.i - b.i, a.j - b.j, a.k - b.k, a.r - b.r);
}

// Quaternion addition.
template <typename Scalar>
inline Quat<Scalar> operator+(const Quat<Scalar>& a, const Quat<Scalar>& b)
{
	return Quat<Scalar>(a.i + b.i, a.j + b.j, a.k + b.k, a.r + b.r);
}

// Quaternion multplication.
template <typename Scalar>
inline Quat<Scalar> operator*(const Quat<Scalar>& q1, const Quat<Scalar>& q2)
{
	return Quat<Scalar>( q1.i * q2.r + q1.j * q2.k - q1.k * q2.j + q1.r * q2.i,
						-q1.i * q2.k + q1.j * q2.r + q1.k * q2.i + q1.r * q2.j,
						 q1.i * q2.j - q1.j * q2.i + q1.k * q2.r + q1.r * q2.k,
						-q1.i * q2.i - q1.j * q2.j - q1.k * q2.k + q1.r * q2.r);
}

////////////////////////////////////////////////////////////////////////////////
// Matrix helpers

template <typename Scalar>
inline Scalar Dot3(const Vector3<Scalar>& a, const Vector3<Scalar>& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

template <typename Scalar>
inline Scalar Dot3(const Vector3<Scalar>& a, const Vector4<Scalar>& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

template <typename Scalar>
inline Scalar Dot3(const Vector4<Scalar>& a, const Vector3<Scalar>& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

template <typename Scalar>
inline Scalar Dot3(const Vector4<Scalar>& a, const Vector4<Scalar>& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

////////////////////////////////////////////////////////////////////////////////
// Matrix

// Create Matrix from three principle axes and a translation.
template <typename Scalar>
inline Matrix<Scalar> Matrix<Scalar>::Axes(const Vector3<Scalar>& xAxis, const Vector3<Scalar>& yAxis, const Vector3<Scalar>& zAxis, const Vector3<Scalar>& trans)
{
	Matrix<Scalar> m;
	m.SetXAxis(xAxis);
	m.SetYAxis(yAxis);
	m.SetZAxis(zAxis);
	m.SetTrans(trans);
	return m;
}

// Create a Matrix from four row vectors.
template <typename Scalar>
inline Matrix<Scalar> Matrix<Scalar>::Rows(const Vector4<Scalar>& row0In, const Vector4<Scalar>& row1In, const Vector4<Scalar>& row2In, const Vector4<Scalar>& row3In)
{
	Matrix<Scalar> m;
	m.row0 = row0In;
	m.row1 = row1In;
	m.row2 = row2In;
	m.row3 = row3In;
	return m;
}

// Create a Matrix from a Quat and a translation.
template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::QuatTrans(const Quat<Scalar>& q, const Vector3<Scalar>& trans)
{
	Matrix<Scalar> m;
	m.SetXAxis(q.Rotate(Vector3<Scalar>(1, 0, 0)));
	m.SetYAxis(q.Rotate(Vector3<Scalar>(0, 1, 0)));
	m.SetZAxis(q.Rotate(Vector3<Scalar>(0, 0, 1)));
	m.SetTrans(Vector3<Scalar>(0, 0, 0));
	return m;
}

// Create a Matrix from a Scale vector, a Quat and a translation.
template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::ScaleQuatTrans(const Vector3<Scalar>& scale, const Quat<Scalar>& q, const Vector3<Scalar>& trans)
{
	Matrix<Scalar> m;
	m.SetXAxis(q.Rotate(Vector3<Scalar>(scale.x, 0, 0)));
	m.SetYAxis(q.Rotate(Vector3<Scalar>(0, scale.y, 0)));
	m.SetZAxis(q.Rotate(Vector3<Scalar>(0, 0, scale.z)));
	m.SetTrans(trans);
	return m;
}

// Create a Matrix from a rotation represented by an axis and an angle.
template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::AxisAngle(const Vector3<Scalar>& axis, Scalar angle)
{
	Matrix<Scalar> m;
	Quat<Scalar> q = Quat<Scalar>::AxisAngle(axis, angle);
	m.SetXAxis(q.Rotate(Vector3<Scalar>(1, 0, 0)));
	m.SetYAxis(q.Rotate(Vector3<Scalar>(0, 1, 0)));
	m.SetZAxis(q.Rotate(Vector3<Scalar>(0, 0, 1)));
	m.SetTrans(Vector3<Scalar>(0, 0, 0));
	return m;
}

// Create an identity Matrix.
template <typename Scalar>
inline Matrix<Scalar> Matrix<Scalar>::Identity()
{
	Matrix<Scalar> m;
	m.row0.Set(1, 0, 0, 0);
	m.row1.Set(0, 1, 0, 0);
	m.row2.Set(0, 0, 1, 0);
	m.row3.Set(0, 0, 0, 1);
	return m;
}

// Create a persective transformation Matrix.
template <typename Scalar>
inline Matrix<Scalar> Matrix<Scalar>::Frustum(Scalar fovy, Scalar aspect, Scalar nearVal, Scalar farVal)
{
	Matrix<Scalar> m;
	Scalar f = 1.0f / tan(fovy / 2);
	m.row0.Set(f / aspect, 0, 0, 0);
	m.row1.Set(0, f, 0, 0);
	m.row2.Set(0, 0, (farVal + nearVal) / (nearVal - farVal), (2 * farVal * nearVal) / (nearVal - farVal));
	m.row3.Set(0, 0, -1, 0);
	return m;
}

// Create an orthograpic projection Matrix.
template <typename Scalar>
inline Matrix<Scalar> Matrix<Scalar>::Ortho(Scalar left, Scalar right, Scalar bottom, Scalar top, Scalar nearVal, Scalar farVal)
{
	Matrix<Scalar> m;
	Scalar tx = -(right + left / right - left);
	Scalar ty = -(top + bottom / top - bottom);
	Scalar tz = -(farVal + nearVal / farVal - nearVal);
	m.row0.Set(2.0f / right - left, 0, 0, tx);
	m.row1.Set(0, 2.0f / top - bottom, 0, ty);
	m.row2.Set(0, 0, -2.0f / farVal - nearVal, tz);
	m.row3.Set(0, 0, 0, 1);
	return m;
}

// Create a look at matrix. (x is forward)
template <typename Scalar>
inline Matrix<Scalar> Matrix<Scalar>::LookAt(const Vector3<Scalar>& eye, const Vector3<Scalar>& target, const Vector3<Scalar>& up)
{
	Matrix<Scalar> m;
	Vector3<Scalar> x = (eye - target).Unit();
	Vector3<Scalar> u = up.Unit();
	Vector3<Scalar> z = Cross(x, u);
	Vector3<Scalar> y = Cross(z, x);
	m.SetXAxis(x);
	m.SetYAxis(y);
	m.SetZAxis(z);
	m.SetTrans(eye);
	return m;
}

// Axes accessors
template <typename Scalar>
inline Vector3<Scalar> Matrix<Scalar>::GetXAxis() const
{
	return Vector3<Scalar>(row0.x, row1.x, row2.x);
}

template <typename Scalar>
inline Vector3<Scalar> Matrix<Scalar>::GetYAxis() const
{
	return Vector3<Scalar>(row0.y, row1.y, row2.y);
}

template <typename Scalar>
inline Vector3<Scalar> Matrix<Scalar>::GetZAxis() const
{
	return Vector3<Scalar>(row0.z, row1.z, row2.z);
}

template <typename Scalar>
inline Vector3<Scalar> Matrix<Scalar>::GetTrans() const
{
	return Vector3<Scalar>(row0.w, row1.w, row2.w);
}

template <typename Scalar>
inline void Matrix<Scalar>::SetXAxis(const Vector3<Scalar>& xAxis)
{
	row0.x = xAxis.x;
	row1.x = xAxis.y;
	row2.x = xAxis.z;
	row3.x = 0.0f;
}

template <typename Scalar>
inline void Matrix<Scalar>::SetYAxis(const Vector3<Scalar>& yAxis)
{
	row0.y = yAxis.x;
	row1.y = yAxis.y;
	row2.y = yAxis.z;
	row3.y = 0.0f;
}

template <typename Scalar>
inline void Matrix<Scalar>::SetZAxis(const Vector3<Scalar>& zAxis)
{
	row0.z = zAxis.x;
	row1.z = zAxis.y;
	row2.z = zAxis.z;
	row3.z = 0.0f;
}

template <typename Scalar>
inline void Matrix<Scalar>::SetTrans(const Vector3<Scalar>& trans)
{
	row0.w = trans.x;
	row1.w = trans.y;
	row2.w = trans.z;
	row3.w = 1.0f;
}

// Multiplies by uniform scale.
template <typename Scalar>
inline void Matrix<Scalar>::SetScale(Scalar scale)
{
	SetXAxis(GetXAxis() * scale);
	SetYAxis(GetYAxis() * scale);
	SetZAxis(GetZAxis() * scale);	
}

// Multiplies by non-uniform scale.
template <typename Scalar>
inline void Matrix<Scalar>::SetScale(const Vector3<Scalar>& scale)
{
	SetXAxis(GetXAxis() * scale.x);
	SetYAxis(GetYAxis() * scale.y);
	SetZAxis(GetZAxis() * scale.z);
}

// Element accessors
template <typename Scalar>
inline const Scalar& Matrix<Scalar>::Elem(int r, int c) const
{
	return ((Scalar*)&row0)[r * 4 + c];
}

template <typename Scalar>
inline Scalar& Matrix<Scalar>::Elem(int r, int c)
{
	return ((Scalar*)&row0)[r * 4 + c];
}

// Column accessor
template <typename Scalar>
inline Vector4<Scalar> Matrix<Scalar>::GetCol(int c) const
{
	return Vector4<Scalar>(row0[c], row1[c], row2[c], row3[c]);
}

// Returns the rotation component of this Matrix.
template <typename Scalar>
Quat<Scalar> Matrix<Scalar>::GetQuat() const
{
	// TODO: add to unit test!

	int x = 0, y = 1, z = 2;
	// create a matrix with no scale
	Matrix<Scalar> m = Axes(GetXAxis().Unit(), GetYAxis().Unit(), GetZAxis().Unit(), Vector3<Scalar>(0,0,0));
	Scalar trace = m.Elem(0,0) + m.Elem(1,1) + m.Elem(2,2);
	Vector4<Scalar> q;
	if (trace > -1.0f)
	{
		int i = x, j = y, k = z;
		if (m.Elem(y,y) > m.Elem(x,x)) { i = y; j = z; k = x; }
		if (m.Elem(z,z) > m.Elem(i,i)) { i = z; j = x; k = y; }
		Scalar r = sqrt(m.Elem(i,i) - m.Elem(j,j) - m.Elem(k,k) + 1);
		q[i] = r / 2;
		q[j] = (m.Elem(i,j) + m.Elem(j,i)) / (2 * r);
		q[k] = (m.Elem(k,i) + m.Elem(i,k)) / (2 * r);
		q[3] = (m.Elem(k,j) + m.Elem(j,k)) / (2 * r);
	}
    else
	{
		q.w = sqrt(1 + trace) / 2;
		q.x = (row2[1] - row1[2]) / (4 * q.w);
		q.y = (row0[2] - row2[0]) / (4 * q.w);
		q.z = (row1[0] - row0[1]) / (4 * q.w);
	}
	return Quat<Scalar>(q.x, q.y, q.z, q.w);
}

// Multiply the 3x3 component of this Matrix with a column vector.
template <typename Scalar>
inline Vector3<Scalar> Matrix<Scalar>::Mul3x3(const Vector3<Scalar>& v) const
{
	return Vector3<Scalar>(Dot3(row0, v), Dot3(row1, v), Dot3(row2, v));
}

// Multiply the 3x4 component of this Matrix with a column vector. (w component of vector is 1.0)
template <typename Scalar>
inline Vector3<Scalar> Matrix<Scalar>::Mul3x4(const Vector3<Scalar>& v) const
{
	return Vector3<Scalar>(Dot3(row0, v) + row0.w, Dot3(row1, v) + row1.w, Dot3(row2, v) + row2.w);
}

// Multiply this Matrix with a column vector.
template <typename Scalar>
inline Vector4<Scalar> Matrix<Scalar>::Mul4x4(const Vector4<Scalar>& v) const
{
	return Vector4<Scalar>(Dot(row0, v), Dot(row1, v), Dot(row2, v), Dot(row3, v));
}

// Returns the transpose of this Matrix
template <typename Scalar>
inline Matrix<Scalar> Matrix<Scalar>::Transpose() const
{
	return Rows(GetCol(0), GetCol(1), GetCol(2), GetCol(3));
}

// If the 3x3 portion of this Matrix is Orthogonal (i.e. columns are orthogonal unit vectors)
// this will return the Inverse of that matrix.
template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::OrthoInverse() const
{
	Matrix<Scalar> r(*this);
	r.SetTrans(Vector3<Scalar>(0, 0, 0));
	r = r.Transpose();
	r.SetTrans(-r.Mul3x4(GetTrans()));
	return r;
}

// Full 4x4 Matrix Inverse, returns Identity if matrix has no inverse.
template <typename Scalar>
Matrix<Scalar> Matrix<Scalar>::FullInverse() const
{
	Matrix<Scalar> m;
	if (::FullInverse((*this), m))
		return m;
	else
		return Identity();
}

// Matrix addition
template <typename Scalar>
Matrix<Scalar> operator+(const Matrix<Scalar>& a, const Matrix<Scalar>& b)
{
	return Matrix<Scalar>::Rows(a.row0 + b.row0, a.row1 + b.row1, a.row2 + b.row2, a.row3 + b.row2);
}

// Matrix subtraction
template <typename Scalar>
Matrix<Scalar> operator-(const Matrix<Scalar>& a, const Matrix<Scalar>& b)
{
	return Matrix<Scalar>::Rows(a.row0 - b.row0, a.row1 - b.row1, a.row2 - b.row2, a.row3 - b.row2);
}

// Matrix multiplication
template <typename Scalar>
Matrix<Scalar> operator*(const Matrix<Scalar>& a, const Matrix<Scalar>& b)
{
	Matrix<Scalar> bt = b.Transpose();
	return Matrix<Scalar>::Rows( Vector4<Scalar>(Dot(a.row0, bt.row0), Dot(a.row0, bt.row1), 
												 Dot(a.row0, bt.row2), Dot(a.row0, bt.row3)), 
								 Vector4<Scalar>(Dot(a.row1, bt.row0), Dot(a.row1, bt.row1), 
												 Dot(a.row1, bt.row2), Dot(a.row1, bt.row3)),
								 Vector4<Scalar>(Dot(a.row2, bt.row0), Dot(a.row2, bt.row1), 
												 Dot(a.row2, bt.row2), Dot(a.row2, bt.row3)),
								 Vector4<Scalar>(Dot(a.row3, bt.row0), Dot(a.row3, bt.row1), 
												 Dot(a.row3, bt.row2), Dot(a.row3, bt.row3)));
}

template <typename T>
void AbaciSwap(T& a, T& b)
{
	T temp = a;
	a = b;
	b = temp;
}

// Full 4x4 Matrix inverse, returns false if Matrix has no inverse.
template <typename Scalar>
bool FullInverse(const Matrix<Scalar>& m, Matrix<Scalar>& result)
{
    // Gaussian-Jordan Elimination, TODO: There are more numerically stable ways to do this...

	Scalar temp[4][8];
	Scalar* row[4];

	// initialize the r pointers
	for (int r = 0; r < 4; ++r)
		row[r] = temp[r];

	// initialize the augmented temp matrix. 
	// the first four columns are from m
	temp[0][0] = m.row0.x; temp[0][1] = m.row0.y; temp[0][2] = m.row0.z; temp[0][3] = m.row0.w;
	temp[1][0] = m.row1.x; temp[1][1] = m.row1.y; temp[1][2] = m.row1.z; temp[1][3] = m.row1.w;
	temp[2][0] = m.row2.x; temp[2][1] = m.row2.y; temp[2][2] = m.row2.z; temp[2][3] = m.row2.w;
	temp[3][0] = m.row3.x; temp[3][1] = m.row3.y; temp[3][2] = m.row3.z; temp[3][3] = m.row3.w;
	// the second four are identity
	temp[0][4] = 1; temp[0][5] = 0; temp[0][6] = 0; temp[0][7] = 0;
	temp[1][4] = 0; temp[1][5] = 1; temp[1][6] = 0; temp[1][7] = 0;
	temp[2][4] = 0; temp[2][5] = 0; temp[2][6] = 1; temp[2][7] = 0;
	temp[3][4] = 0; temp[3][5] = 0; temp[3][6] = 0; temp[3][7] = 1;

	// bubble up row with largest leading number (partial pivot)
	if (fabs(row[0][0]) < fabs(row[1][0]))
		AbaciSwap(row[0], row[1]);
	if (fabs(row[0][0]) < fabs(row[2][0]))
		AbaciSwap(row[0], row[2]);
	if (fabs(row[0][0]) < fabs(row[3][0]))
		AbaciSwap(row[0], row[3]);

	if (fabs(row[0][0]) < 0.00001)  // column is all zeros, there is no inverse.
	{		
		return false;
	}

	// mult row[0] by 1/row[0][0].  To introduce a leading 1.
	Scalar s = 1 / row[0][0];
	for (int c = 0; c < 8; ++c)	row[0][c] *= s;

	// add multiples of top row to lower rows so that all entries below leading 1 become zeros.
	for (int r = 1; r < 4; ++r)
	{
		Scalar s = row[r][0];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[0][c];
	}

	// move row with largest leading number 
	if (fabs(row[1][1]) < fabs(row[2][1]))
		AbaciSwap(row[1], row[2]);
	if (fabs(row[1][1]) < fabs(row[3][1]))
		AbaciSwap(row[1], row[3]);

	if (fabs(row[1][1]) < 0.00001)  // column is all zeros, there is no inverse.
	{
		return false;
	}

	// mult row[1] by 1/row[1][1].  To introduce a leading 1.
	s = 1 / row[1][1];
	for (int c = 0; c < 8; ++c)	row[1][c] *= s;
	
	// add multiples of top row to lower rows so that all entries below leading 1 become zeros.
	for (int r = 2; r < 4; ++r)
	{
		Scalar s = row[r][1];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[1][c];
	}

	// move row with largest leading number 
	if (fabs(row[2][2]) < fabs(row[3][2]))
		AbaciSwap(row[2], row[3]);

	if (fabs(row[2][2]) < 0.00001)  // column is all zeros, there is no inverse.
	{
		return false;
	}

	// mult row[2] by 1/row[2][2].  To introduce a leading 1.
	s = 1.0f / row[2][2];
	for (int c = 0; c < 8; ++c)	row[2][c] *= s;
	
	// add multiples of top row to lower rows so that all entries below leading 1 become zeros.
	for (int r = 3; r < 4; ++r)
	{
		Scalar s = row[r][2];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[2][c];
	}

	// mult row[3] by 1/row[3][3]. To introduce a leading 1.
	s = 1.0f / row[3][3];
	for (int c = 0; c < 8; ++c)	row[3][c] *= s;

	// at this point row matrix should be in row-echelon form.

	// add multiples of row[3] to above rows, to zero out that column
	for (int r = 0; r < 3; ++r)
	{
		Scalar s = row[r][3];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[3][c];
	}

	// add multiples of row[2] to above rows, to zero out that column
	for (int r = 0; r < 2; ++r)
	{
		Scalar s = row[r][2];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[2][c];
	}

	// add multiples of row[1] to above row, to zero out that column
	for (int r = 0; r < 1; ++r)
	{
		Scalar s = row[r][1];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[1][c];
	}

	// init result
	result.row0.x = row[0][4]; result.row0.y = row[0][5]; result.row0.z = row[0][6]; result.row0.w = row[0][7];
	result.row1.x = row[1][4]; result.row1.y = row[1][5]; result.row1.z = row[1][6]; result.row1.w = row[1][7];
	result.row2.x = row[2][4]; result.row2.y = row[2][5]; result.row2.z = row[2][6]; result.row2.w = row[2][7];
	result.row3.x = row[3][4]; result.row3.y = row[3][5]; result.row3.z = row[3][6]; result.row3.w = row[3][7];

	return true;	
}


// Print to stdout.
template <typename Scalar>
void PrintMatrix(const Matrix<Scalar>& m)
{
	printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.row0.x, m.row0.y, m.row0.z, m.row0.w);
	printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.row1.x, m.row1.y, m.row1.z, m.row1.w);
	printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.row2.x, m.row2.y, m.row2.z, m.row2.w);
	printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.row3.x, m.row3.y, m.row3.z, m.row3.w);
}

////////////////////////////////////////////////////////////////////////////////
// Complex

// Construct from two floats
template <typename Scalar>
inline Complex<Scalar>::Complex(Scalar rIn, Scalar iIn) : r(rIn), i(iIn) {}

// Length
template <typename Scalar>
inline Scalar Complex<Scalar>::Len() const
{
	return sqrt(Dot(*this, *this));
}

// Square of length
template <typename Scalar>
inline Scalar Complex<Scalar>::LenSq() const
{
	return Dot(*this, *this);
}

// Dot product
template <typename Scalar>
inline Scalar Dot(const Complex<Scalar>& a, const Complex<Scalar>& b)
{
	return a.r * b.r + b.i * b.i;
}

// Complex conjugate
template <typename Scalar>
inline Complex<Scalar> operator~(const Complex<Scalar>& a)
{
	return Complex<Scalar>(a.r, -a.i);
}

// Unary minus.
template <typename Scalar>
inline Complex<Scalar> operator-(const Complex<Scalar>& a)
{
	return Complex<Scalar>(-a.r, -a.i);
}

// Unary plus.
template <typename Scalar>
inline Complex<Scalar> operator+(const Complex<Scalar>& a)
{
	return Complex<Scalar>(a.r, a.i);
}

// Complex addition.
template <typename Scalar>
inline Complex<Scalar> operator+(const Complex<Scalar>& a, const Complex<Scalar>& b)
{
	return Complex<Scalar>(a.r + b.r, a.i + b.i);
}

// Complex subtraction.
template <typename Scalar>
inline Complex<Scalar> operator-(const Complex<Scalar>& a, const Complex<Scalar>& b)
{
	return Complex<Scalar>(a.r - b.r, a.i - b.i);
}

// Complex multiplication.
template <typename Scalar>
inline Complex<Scalar> operator*(const Complex<Scalar>& a, const Complex<Scalar>& b)
{
	Scalar aa = a.r;
	Scalar bb = a.i;
	Scalar cc = b.r;
	Scalar dd = b.i;

	return Complex<Scalar>(aa * cc - (bb * dd), aa * dd + bb * cc);
}

// Multiplication by a real number.
template <typename Scalar>
inline Complex<Scalar> operator*(Scalar scalar, const Complex<Scalar>& c)
{
	return Complex<Scalar>(scalar, 0) * c;
}

// Multiplication by a real number.
template <typename Scalar>
inline Complex<Scalar> operator*(Complex<Scalar>& c, Scalar scalar)
{
	return c * Complex<Scalar>(scalar, 0);
}

// Complex division.
template <typename Scalar>
inline Complex<Scalar> operator/(const Complex<Scalar>& a, const Complex<Scalar>& b)
{
	Scalar aa = a.r;
	Scalar bb = a.i;
	Scalar cc = b.r;
	Scalar dd = b.i;
	Scalar denom = cc * cc + dd * dd;

	return Complex<Scalar>((aa * cc + bb * dd) / denom, (bb * cc - aa * dd) / denom);
}

// e ^ 0 + xi
template <typename Scalar>
inline Complex<Scalar> ExpI(Scalar x)
{
	return Complex<Scalar>(cos(x), sin(x));
}

#endif


	
