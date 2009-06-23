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

inline float DegToRad(float deg)
{ 
	return deg * ((2 * PI) / 180.0f);
}

inline float RadToDeg(float rad)
{
	return rad * (180.0f / (2 * PI));
}

inline float Clamp(float value, float min, float max)
{
	float result = value;
	if (value > max)
		result = max;
	if (value < min)
		result = min;
	return result; 
}

// Limits angle between -PI & PI using modulus arithmetic
inline float LimitPi(float theta)
{
	float phi = Mod2Pi(theta);
	if (phi > PI)
		return phi - 2 * PI;
	else
		return phi;
}

// Limits angle between zero & PI using modulus arithmetic
inline float Mod2Pi(float theta)
{
	return fmod(theta, 2.0f * PI);
}

inline bool FuzzyEqual(float rhs, float lhs, float epsilon)
{
	return fabs(rhs - lhs) <= epsilon;
}

inline Vector2::Vector2(float xIn, float yIn) 
{ 
	x = xIn; 
	y = yIn; 
}

inline void Vector2::Set(float xIn, float yIn)
{
	x = xIn;
	y = yIn;
}

inline void Vector2::SetZero()
{
	x = 0;
	y = 0;
}

inline Vector2 Vector2::operator-() const
{
	return Vector2(-x, -y);
}

inline float Vector2::Dot(const Vector2& b) const
{
	return x * b.x + y * b.y;
}

inline Vector2 Vector2::Unit() const
{
	return *this / Len();
}

inline float Vector2::Len() const
{
	return sqrt((*this).Dot(*this));
}

inline float Vector2::LenSq() const
{
	return (*this).Dot(*this);
}

inline Vector2 operator-(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x - b.x, a.y - b.y);
}

inline Vector2 operator+(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x + b.x, a.y + b.y);
}

inline Vector2 operator*(const Vector2& v, float factor)
{
	return Vector2(factor * v.x, factor * v.y);
}

inline Vector2 operator*(float factor, const Vector2& v)
{
	return Vector2(factor * v.x, factor * v.y);
}

inline Vector2 operator*(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x * b.x, a.y * b.y);
}

inline Vector2 operator/(const Vector2& v, float denominator)
{
	return Vector2(v.x / denominator, v.y / denominator);
}

inline Vector2 operator/(float numerator, const Vector2& v)
{
	return Vector2(numerator / v.x, numerator / v.y);
}

inline Vector2 operator/(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x / b.x, a.y / b.y);
}

inline Vector3::Vector3(const Vector2& v, float zIn) 
{ 
	x = v.x; 
	y = v.y; 
	z = zIn; 
}

inline Vector3::Vector3(float xIn, float yIn, float zIn) 
{ 
	x = xIn; 
	y = yIn; 
	z = zIn; 
}

inline void Vector3::Set(float xIn, float yIn, float zIn)
{
	x = xIn;
	y = yIn;
	z = zIn;
}

inline void Vector3::SetZero()
{
	x = 0;
	y = 0;
	z = 0;
}

inline Vector3 operator-(const Vector3& v)
{
	return Vector3(-v.x, -v.y, -v.z);
}

inline Vector3 operator-(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Vector3 operator+(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline Vector3 operator*(const Vector3& v, float factor)
{
	return Vector3(v.x * factor, v.y * factor, v.z * factor);
}

inline Vector3 operator*(float factor, const Vector3& v)
{
	return Vector3(v.x * factor, v.y * factor, v.z * factor);
}

inline float operator*(const Vector3& a, const Vector3& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Vector3 CompMul(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x * b.x, a.y * b.y, a.z * b.z);
}

inline Vector3 operator/(const Vector3& v, float denominator)
{
	return Vector3(v.x / denominator, v.y / denominator, v.z / denominator);
}

inline Vector3 operator/(float numerator, const Vector3& v)
{
	return Vector3(numerator / v.x, numerator / v.y, numerator / v.z);
}

inline Vector3 CompDiv(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x / b.x, a.y / b.y, a.z / b.z);
}

inline float Len(const Vector3& v)
{
	return sqrt(v * v);
}

inline Vector3 UnitVec(const Vector3& v)
{
	float len = Len(v);
	return v / len;
}

inline Vector3 operator%(const Vector3& a, const Vector3& b)
{
	return Vector3((a.y * b.z) - (a.z * b.y), 
	               (a.z * b.x) - (a.x * b.z), 
	               (a.x * b.y) - (a.y * b.x));
}

inline Vector4::Vector4(const Vector3& v, float wIn) 
{ 
	x = v.x; 
	y = v.y; 
	z = v.z; 
    w = wIn; 
}

inline Vector4::Vector4(float xIn, float yIn, float zIn, float wIn) 
{ 
	x = xIn; 
	y = yIn; 
	z = zIn; 
	w = wIn; 
}

inline void Vector4::Set(float xIn, float yIn, float zIn, float wIn)
{
	x = xIn;
	y = yIn;
	z = zIn;
	w = wIn;
}

inline void Vector4::SetZero()
{
	x = 0;
	y = 0;
	z = 0;
	w = 0;
}

inline float Vector4::operator[](int i) const
{
	return *(&x + i);
}

inline float& Vector4::operator[](int i)
{
	return *(&x + i);
}

inline Vector4 operator-(const Vector4& v)
{
	return Vector4(-v.x, -v.y, -v.z, -v.w);
}

inline Vector4 operator-(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

inline Vector4 operator+(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

inline Vector4 operator*(const Vector4& v, float factor)
{
	return Vector4(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

inline Vector4 operator*(float factor, const Vector4& v)
{
	return Vector4(factor * v.x, factor * v.y, factor * v.z, factor * v.w);	
}

inline float operator*(const Vector4& a, const Vector4& b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

inline Vector4 CompMul(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

inline Vector4 operator/(const Vector4& v, float denominator)
{
	return Vector4(v.x / denominator, v.y / denominator, v.z / denominator, v.w / denominator);
}

inline Vector4 operator/(float numerator, const Vector4& v)
{
	return Vector4(numerator / v.x, numerator / v.y, numerator / v.z, numerator / v.w);
}

inline Vector4 CompDiv(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}

inline float Len(const Vector4& v)
{
	return sqrt(v * v);
}

inline Vector4 UnitVec(const Vector4& v)
{
	float len = Len(v);
	return v / len;
}

inline Quat::Quat(float xIn, float yIn, float zIn, float wIn) 
{ 
	x = xIn; 
	y = yIn; 
	z = zIn; 
	w = wIn; 
}

inline Quat::Quat(const Vector3& axis, float angle)
{
	Vector3 n = UnitVec(axis) * sin(angle/2.0f);
	x = n.x;
	y = n.y;
	z = n.z;
	w = cos(angle/2.0f);
}

inline Vector3 Rotate(const Quat& q, const Vector3& v)
{
	Quat r = q * Quat(v.x, v.y, v.z, 0.0f) * ~q;
	return Vector3(r.x, r.y, r.z);
}

inline Quat operator~(const Quat& v)
{
	return Quat(-v.x, -v.y, -v.z, v.w);
}

inline Quat operator-(const Quat& v)
{
	return Quat(-v.x, -v.y, -v.z, -v.w);
}

inline Quat operator-(const Quat& a, const Quat& b)
{
	return Quat(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

inline Quat operator+(const Quat& a, const Quat& b)
{
	return Quat(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

inline Quat operator*(const Quat& v, float factor)
{
	return Quat(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

inline Quat operator*(float factor, const Quat& v)
{
	return Quat(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

inline Quat operator*(const Quat& q1, const Quat& q2)
{
	return Quat( q1.x * q2.w + q1.y * q2.z - q1.z * q2.y + q1.w * q2.x,
				-q1.x * q2.z + q1.y * q2.w + q1.z * q2.x + q1.w * q2.y,
				 q1.x * q2.y - q1.y * q2.x + q1.z * q2.w + q1.w * q2.z,
				-q1.x * q2.x - q1.y * q2.y - q1.z * q2.z + q1.w * q2.w);
}

inline Quat CompMul(const Quat& a, const Quat& b)
{
	return Quat(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w);
}

inline Quat operator/(const Quat& v, float denominator)
{
	return Quat(v.x / denominator, v.y / denominator, v.z / denominator, v.w / denominator);
}

inline Quat operator/(float numerator, const Quat& v)
{
	return Quat(numerator / v.x, numerator / v.y, numerator / v.z, numerator / v.w);
}

inline Quat CompDiv(const Quat& a, const Quat& b)
{
	return Quat(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w);
}

inline float Len(const Quat& v)
{
	return sqrt(static_cast<const Vector4&>(v) * static_cast<const Vector4&>(v));	
}

inline Quat UnitVec(const Quat& v)
{
	float len = Len(v);
	return v / len;	
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
	Vector3 x = UnitVec(eye - target);
	Vector3 u = UnitVec(up);
	Vector3 z = x % u;
	Vector3 y = z % x;
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
	return Vector4(m.row0 * v, m.row1 * v, m.row2 * v, m.row3 * v);
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


	
