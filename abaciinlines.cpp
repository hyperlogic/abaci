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

// Convert from degrees to radians
inline float DegToRad(float deg)
{ 
	return deg * (PI / 180.0f);
}

// Convert from radians to degrees.
inline float RadToDeg(float rad)
{
	return rad * (180.0f / PI);
}

// Clamps value between min & max.
inline float Clamp(float value, float min, float max)
{
	float result = value;
	if (value > max)
		result = max;
	else if (value < min)
		result = min;
	return result; 
}

// Limits angle between -PI & PI using modulus arithmetic
inline float LimitPi(float theta)
{
	return fmod(theta + PI, 2.0f * PI) - PI;
}

// Limits angle between zero & PI using modulus arithmetic
inline float Mod2Pi(float theta)
{
	return fmod(theta, 2.0f * PI);
}

// Fuzzy comparison between two float values.
inline bool FuzzyEqual(float rhs, float lhs, float epsilon)
{
	return fabs(rhs - lhs) <= epsilon;
}

// Construct from two floats
inline Vector2::Vector2(float xIn, float yIn) 
{ 
	x = xIn; 
	y = yIn; 
}

// Set from two floats
inline void Vector2::Set(float xIn, float yIn)
{
	x = xIn;
	y = yIn;
}

// Sets all elements to zero.
inline void Vector2::SetZero()
{
	x = 0;
	y = 0;
}

// Returns a vector with same direction but unit length.
inline Vector2 Vector2::Unit() const
{
	return *this / Len();
}

// Returns vector length.
inline float Vector2::Len() const
{
	return sqrt(Dot(*this, *this));
}

// Returns length squared.
inline float Vector2::LenSq() const
{
	return Dot(*this, *this);
}

// Dot product of two vectors.
inline float Dot(const Vector2& a, const Vector2& b)
{
	return a.x * b.x + a.y * b.y;
}

// Linear interpolation between two vectors
inline Vector2 Lerp(const Vector2& a, const Vector2& b, float t)
{
    return a + (b - a) * t;
}

// Unary minus
inline Vector2 operator-(const Vector2& a)
{
	return Vector2(-a.x, -a.y);
}

// Vector subtraction.
inline Vector2 operator-(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x - b.x, a.y - b.y);
}

// Vector addition.
inline Vector2 operator+(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x + b.x, a.y + b.y);
}

// Multplies all elements of a vector by a scalar.
inline Vector2 operator*(const Vector2& v, float scalar)
{
	return Vector2(scalar * v.x, scalar * v.y);
}

// Multplies all elements of a vector by a scalar.
inline Vector2 operator*(float scalar, const Vector2& v)
{
	return Vector2(scalar * v.x, scalar * v.y);
}

// Vector multiplication.
inline Vector2 operator*(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x * b.x, a.y * b.y);
}

// Divides all elements of a vector by a scalar.
inline Vector2 operator/(const Vector2& v, float denominator)
{
	return Vector2(v.x / denominator, v.y / denominator);
}

// Multiplies a scalar to the reciprical of all elements in a vector.
inline Vector2 operator/(float numerator, const Vector2& v)
{
	return Vector2(numerator / v.x, numerator / v.y);
}

// Vector division.
inline Vector2 operator/(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x / b.x, a.y / b.y);
}

// Construct from a Vector2 and a float
inline Vector3::Vector3(const Vector2& v, float zIn) 
{ 
	x = v.x; 
	y = v.y; 
	z = zIn; 
}

// Construct from three floats
inline Vector3::Vector3(float xIn, float yIn, float zIn) 
{ 
	x = xIn; 
	y = yIn; 
	z = zIn; 
}

// Set from three floats
inline void Vector3::Set(float xIn, float yIn, float zIn)
{
	x = xIn;
	y = yIn;
	z = zIn;
}

// Set all elements to zero
inline void Vector3::SetZero()
{
	x = 0;
	y = 0;
	z = 0;
}

// Returns a vector with same direction but unit length.
inline Vector3 Vector3::Unit() const
{
	return (*this) / Len();
}

// Returns vector length.
inline float Vector3::Len() const
{
	return sqrt(Dot(*this, *this));
}

// Returns vector length squared.
inline float Vector3::LenSq() const
{
	return Dot(*this, *this);
}

// Dot product of two vectors.
inline float Dot(const Vector3& a, const Vector3& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Cross product of two vectors.
inline Vector3 Cross(const Vector3& a, const Vector3& b)
{
	return Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

// Linear interpolation between two vectors
inline Vector3 Lerp(const Vector3& a, const Vector3& b, float t)
{
    return a + (b - a) * t;
}

// Unary minus
inline Vector3 operator-(const Vector3& a)
{
	return Vector3(-a.x, -a.y, -a.z);
}

// Vector subtraction.
inline Vector3 operator-(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

// Vector addition.
inline Vector3 operator+(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

// Multplies all elements of a vector by a scalar.
inline Vector3 operator*(const Vector3& v, float factor)
{
	return Vector3(v.x * factor, v.y * factor, v.z * factor);
}

// Multplies all elements of a vector by a scalar.
inline Vector3 operator*(float factor, const Vector3& v)
{
	return Vector3(v.x * factor, v.y * factor, v.z * factor);
}

// Vector multiplication
inline Vector3 operator*(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x * b.x, a.y * b.y, a.z * b.z);
}

// Divides all elements of a vector by a scalar.
inline Vector3 operator/(const Vector3& v, float denominator)
{
	return Vector3(v.x / denominator, v.y / denominator, v.z / denominator);
}

// Multiplies a scalar to the reciprical of all elements in a vector.
inline Vector3 operator/(float numerator, const Vector3& v)
{
	return Vector3(numerator / v.x, numerator / v.y, numerator / v.z);
}

// Vector division.
inline Vector3 operator/(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x / b.x, a.y / b.y, a.z / b.z);
}

// Construct from a Vector3 and a float.
inline Vector4::Vector4(const Vector3& v, float wIn) 
{ 
	x = v.x; 
	y = v.y; 
	z = v.z; 
    w = wIn; 
}

// Construct from four floats.
inline Vector4::Vector4(float xIn, float yIn, float zIn, float wIn) 
{ 
	x = xIn; 
	y = yIn; 
	z = zIn; 
	w = wIn; 
}

// Set from four floats.
inline void Vector4::Set(float xIn, float yIn, float zIn, float wIn)
{
	x = xIn;
	y = yIn;
	z = zIn;
	w = wIn;
}

// Set all elements to zero.
inline void Vector4::SetZero()
{
	x = 0;
	y = 0;
	z = 0;
	w = 0;
}

// Returns a vector with same direction but unit length.
inline Vector4 Vector4::Unit() const
{
	return (*this) / Len();
}

// Returns vector length.
inline float Vector4::Len() const
{
	return sqrt(Dot(*this, *this));
}

// Returns length squared.
inline float Vector4::LenSq() const
{
	return Dot(*this, *this);
}

// const array accessor
inline float Vector4::operator[](int i) const
{
	return *(&x + i);
}

// array accessor
inline float& Vector4::operator[](int i)
{
	return *(&x + i);
}

// Dot product of two vectors.
inline float Dot(const Vector4& a, const Vector4& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

// Linear interpolation between two vectors
inline Vector4 Lerp(const Vector4& a, const Vector4& b, float t)
{
	return a + (b - a) * t;
}

// Unary minus.
inline Vector4 operator-(const Vector4& v)
{
	return Vector4(-v.x, -v.y, -v.z, -v.w);
}

// Vector subtraction.
inline Vector4 operator-(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

// Vector addition.
inline Vector4 operator+(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

// Multplies all elements of a vector by a scalar.
inline Vector4 operator*(const Vector4& v, float factor)
{
	return Vector4(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

// Multplies all elements of a vector by a scalar.
inline Vector4 operator*(float factor, const Vector4& v)
{
	return Vector4(factor * v.x, factor * v.y, factor * v.z, factor * v.w);	
}

// Vector multiplication.
inline Vector4 operator*(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);	
}

// Divides all elements of a vector by a scalar.
inline Vector4 operator/(const Vector4& v, float denominator)
{
	return Vector4(v.x / denominator, v.y / denominator, v.z / denominator, v.w / denominator);
}

// Multiplies a scalar to the reciprical of all elements in a vector.
inline Vector4 operator/(float numerator, const Vector4& v)
{
	return Vector4(numerator / v.x, numerator / v.y, numerator / v.z, numerator / v.w);
}

// Vector division.
inline Vector4 operator/(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}

// Create from axis and angle
inline Quat Quat::AxisAngle(const Vector3& axis, float angle)
{
	Vector3 n = axis.Unit() * sin(angle/2.0f);
	float w = cos(angle/2.0f);
	return Quat(n.x, n.y, n.z, w);
}

// Construct from four floats.
inline Quat::Quat(float iIn, float jIn, float kIn, float rIn)
{
	i = iIn;
	j = jIn;
	k = kIn;
	r = rIn;
}

// Set all elements to zero.
inline void Quat::SetZero()
{
	i = j = k = r = 0.0f;
}

// Returns this quaternion normalized.
inline Quat Quat::Unit() const
{
	float len = Len();
	return Quat(i / len, j / len, k / len, r / len);
}

// Returns quat length.
inline float Quat::Len() const
{
	return sqrt(Dot(*this, *this));
}

// Returns quat length squared.
inline float Quat::LenSq() const
{
	return Dot(*this, *this);
}

// Rotate a vector.
inline Vector3 Quat::Rotate(const Vector3& v) const
{
	Quat r = (*this) * Quat(v.x, v.y, v.z, 0.0f) * ~(*this);
	return Vector3(r.i, r.j, r.k);
}

// Dot product
inline float Dot(const Quat& a, const Quat& b)
{
	return a.i * b.i + a.j * b.j + a.k * b.k + a.r * b.r;
}

// Quaternion conjugate
inline Quat operator~(const Quat& v)
{
	return Quat(-v.i, -v.j, -v.k, v.r);
}

// Unary minus.
inline Quat operator-(const Quat& v)
{
	return Quat(-v.i, -v.j, -v.k, -v.r);
}

// Quaternion subtraction.
inline Quat operator-(const Quat& a, const Quat& b)
{
	return Quat(a.i - b.i, a.j - b.j, a.k - b.k, a.r - b.r);
}

// Quaternion addition.
inline Quat operator+(const Quat& a, const Quat& b)
{
	return Quat(a.i + b.i, a.j + b.j, a.k + b.k, a.r + b.r);
}

// Quaternion multplication.
inline Quat operator*(const Quat& q1, const Quat& q2)
{
	return Quat( q1.i * q2.r + q1.j * q2.k - q1.k * q2.j + q1.r * q2.i,
				-q1.i * q2.k + q1.j * q2.r + q1.k * q2.i + q1.r * q2.j,
				 q1.i * q2.j - q1.j * q2.i + q1.k * q2.r + q1.r * q2.k,
				-q1.i * q2.i - q1.j * q2.j - q1.k * q2.k + q1.r * q2.r);
}

inline float Dot3(const Vector3& a, const Vector3& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

inline float Dot3(const Vector3& a, const Vector4& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

inline float Dot3(const Vector4& a, const Vector3& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

inline float Dot3(const Vector4& a, const Vector4& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

inline Matrix::Matrix(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis, const Vector3& trans)
{
	SetXAxis(xAxis);
	SetYAxis(yAxis);
	SetZAxis(zAxis);
	SetTrans(trans);
}

inline Matrix::Matrix(const Vector4& row0In, const Vector4& row1In, const Vector4& row2In, const Vector4& row3In)
{
	row0 = row0In;
	row1 = row1In;
	row2 = row2In;
	row3 = row3In;
}

inline void Matrix::MakeIdent()
{
	row0.Set(1.0f, 0.0f, 0.0f, 0.0f);
	row1.Set(0.0f, 1.0f, 0.0f, 0.0f);
	row2.Set(0.0f, 0.0f, 1.0f, 0.0f);
	row3.Set(0.0f, 0.0f, 0.0f, 1.0f);
}

inline void Matrix::MakeFrustum(float fovy, float aspect, float nearVal, float farVal)
{
	float f = 1.0f / tan(fovy / 2.0f);
	row0.Set(f / aspect, 0, 0, 0);
	row1.Set(0, f, 0, 0);
	row2.Set(0, 0, (farVal + nearVal) / (nearVal - farVal), (2.0f * farVal * nearVal) / (nearVal - farVal));
	row3.Set(0, 0, -1, 0);
}

inline void Matrix::MakeOrtho(float left, float right, float bottom, float top, float nearVal, float farVal)
{
	float tx = -(right + left / right - left);
	float ty = -(top + bottom / top - bottom);
	float tz = -(farVal + nearVal / farVal - nearVal);
	row0.Set(2.0f / right - left, 0, 0, tx);
	row1.Set(0, 2.0f / top - bottom, 0, ty);
	row2.Set(0, 0, -2.0f / farVal - nearVal, tz);
	row3.Set(0, 0, 0, 1);
}

inline void Matrix::MakeLookAt(const Vector3& eye, const Vector3& target, const Vector3& up)
{
	Vector3 x = (eye - target).Unit();
	Vector3 u = up.Unit();
	Vector3 z = Cross(x, u);
	Vector3 y = Cross(z, x);
	SetXAxis(x);
	SetYAxis(y);
	SetZAxis(z);
	SetTrans(eye);
}

inline Vector3 Matrix::GetXAxis() const
{
	return Vector3(row0.x, row1.x, row2.x);
}

inline Vector3 Matrix::GetYAxis() const
{
	return Vector3(row0.y, row1.y, row2.y);
}

inline Vector3 Matrix::GetZAxis() const
{
	return Vector3(row0.z, row1.z, row2.z);
}

inline Vector3 Matrix::GetTrans() const
{
	return Vector3(row0.w, row1.w, row2.w);
}

inline void Matrix::SetXAxis(const Vector3 xAxis)
{
	row0.x = xAxis.x;
	row1.x = xAxis.y;
	row2.x = xAxis.z;
	row3.x = 0.0f;
}

inline void Matrix::SetYAxis(const Vector3 yAxis)
{
	row0.y = yAxis.x;
	row1.y = yAxis.y;
	row2.y = yAxis.z;
	row3.y = 0.0f;
}

inline void Matrix::SetZAxis(const Vector3 zAxis)
{
	row0.z = zAxis.x;
	row1.z = zAxis.y;
	row2.z = zAxis.z;
	row3.z = 0.0f;
}

inline void Matrix::SetTrans(const Vector3 trans)
{
	row0.w = trans.x;
	row1.w = trans.y;
	row2.w = trans.z;
	row3.w = 1.0f;
}

inline void Matrix::SetScale(float scale)
{
	SetXAxis(GetXAxis() * scale);
	SetYAxis(GetYAxis() * scale);
	SetZAxis(GetZAxis() * scale);	
}

inline void Matrix::SetScale(const Vector3 scale)
{
	SetXAxis(GetXAxis() * scale.x);
	SetYAxis(GetYAxis() * scale.y);
	SetZAxis(GetZAxis() * scale.z);
}

inline const float& Matrix::Elem(int r, int c) const
{
	return ((float*)&row0)[r * 4 + c];
}

inline float& Matrix::Elem(int r, int c)
{
	return ((float*)&row0)[r * 4 + c];
}

inline Vector4 Matrix::GetCol(int c) const
{
	return Vector4(row0[c], row1[c], row2[c], row3[c]);
}

inline Vector3 Transform3x3(const Matrix& m, const Vector3& v)
{
	return Vector3(Dot3(m.row0, v), Dot3(m.row1, v), Dot3(m.row2, v));
}

inline Vector3 Transform3x4(const Matrix& m, const Vector3& v)
{
	return Vector3(Dot3(m.row0, v) + m.row0.w, Dot3(m.row1, v) + m.row1.w, Dot3(m.row2, v) + m.row2.w);
}

inline Vector4 Transform4x4(const Matrix& m, const Vector4& v)
{
	return Vector4(Dot(m.row0, v), Dot(m.row1, v), Dot(m.row2, v), Dot(m.row3, v));
}

inline Matrix Transpose(const Matrix& m)
{
	return Matrix(m.GetCol(0), m.GetCol(1), m.GetCol(2), m.GetCol(3));
}

inline Complex::Complex(float realIn, float imagIn) : real(realIn), imag(imagIn)
{

}

inline Complex Complex::conj() const
{
	return Complex(real, -imag);
}

inline float Complex::len() const
{
	return sqrt(dot(*this, *this));
}

inline Complex operator+(const Complex& a, const Complex& b)
{
	return Complex(a.real + b.real, a.imag + b.imag);
}

inline Complex operator-(const Complex& a, const Complex& b)
{
	return Complex(a.real - b.real, a.imag - b.imag);
}

inline Complex operator-(const Complex& a)
{
	return Complex(-a.real, -a.imag);
}

inline Complex operator+(const Complex& a)
{
	return Complex(a.real, a.imag);
}

inline Complex operator*(const Complex& a, const Complex& b)
{
	float aa = a.real;
	float bb = a.imag;
	float cc = b.real;
	float dd = b.imag;

	return Complex(aa * cc - (bb * dd), aa * dd + bb * cc);
}

inline Complex operator*(float f, const Complex& c)
{
	return Complex(f,0) * c;
}

inline Complex operator*(Complex& c, float f)
{
	return c * Complex(f,0);
}

inline Complex operator/(const Complex& a, const Complex& b)
{
	float aa = a.real;
	float bb = a.imag;
	float cc = b.real;
	float dd = b.imag;
	float denom = cc * cc + dd * dd;

	return Complex((aa * cc + bb * dd) / denom, (bb * cc - aa * dd) / denom);
}

inline Complex expi(float imag)
{
	return Complex(cos(imag), sin(imag));
}

inline float dot(const Complex& a, const Complex& b)
{
	return (a.real * b.real) + (b.imag * b.imag);
}

#endif


	
