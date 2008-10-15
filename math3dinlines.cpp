#include "math3d.h"
#include <stdio.h>

#ifdef USE_SSE2

inline const float& VectorBase::operator[](int i) const 
{ 
	return reinterpret_cast<float*>(&data)[i]; 
}

inline float& VectorBase::operator[](int i) 
{ 
	return reinterpret_cast<float*>(&data)[i]; 
}

#else

inline const float& VectorBase::operator[](int i) const 
{ 
//	return data[i]; 
	return (&x)[i];
}

inline float& VectorBase::operator[](int i) 
{ 
//	return data[i]; 
	return (&x)[i];
}

#endif

inline const float& VectorBase::X() const 
{ 
	return (*this)[0];
}

inline const float& VectorBase::Y() const
{ 
	return (*this)[1]; 
}

inline const float& VectorBase::Z() const
{ 
	return (*this)[2]; 
}

inline const float& VectorBase::W() const 
{ 
	return (*this)[3]; 
}

inline float& VectorBase::X()
{ 
	return (*this)[0];
}

inline float& VectorBase::Y()
{ 
	return (*this)[1]; 
}

inline float& VectorBase::Z()
{ 
	return (*this)[2]; 
}

inline float& VectorBase::W()
{ 
	return (*this)[3]; 
}

inline float DegToRad(float deg)
{ 
	return deg * ((2 * PI) / 180.0f);
}

inline float RadToDeg(float rad)
{
	return rad * (180.0f / (2 * PI));
}

inline float Dot2(const VectorBase& a, const VectorBase& b)
{
	return (a.x * b.x) + (a.y * b.y);
}

inline float Dot3(const VectorBase& a, const VectorBase& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

inline float Dot4(const VectorBase& a, const VectorBase& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w);
}

inline Vector2::Vector2(float xIn, float yIn) 
{ 
	X() = xIn; 
	Y() = yIn; 
}

inline void Vector2::Set(float xIn, float yIn)
{
	x = xIn;
	y = yIn;
}

inline Vector2& Vector2::operator+=(const Vector2& v)
{
	x += v.x;
	y += v.y;
	return *this;
}

inline Vector2 operator-(const Vector2& v)
{
	return Vector2(-v.x, -v.y);
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

inline float operator*(const Vector2& a, const Vector2& b)
{
	return Dot2(a, b);
}

inline float Len(const Vector2& v)
{
	return sqrt(v * v);
}

inline Vector2 UnitVec(const Vector2& v)
{
	float len = Len(v);
	return v * (1.0f / len);
}

inline Vector3::Vector3(const Vector2& v, float zIn) 
{ 
	X() = v.X(); 
	Y() = v.Y(); 
	Z() = zIn; 
}

inline Vector3::Vector3(float xIn, float yIn, float zIn) 
{ 
	X() = xIn; 
	y = yIn; 
	z = zIn; 
}

inline void Vector3::Set(float xIn, float yIn, float zIn)
{
	x = xIn;
	y = yIn;
	z = zIn;
}

inline void Vector3::Normalize()
{
	float invLen = 1.0f / sqrt((x * x) + (y * y) + (z * z));
 	x *= invLen;
	y *= invLen;
	z *= invLen;
}

inline void Vector3::Mul3x4(const Matrix& m)
{
	Vector3 v(*this);
	x = Dot3(m.row0, v) + m.row0.w;
	y = Dot3(m.row1, v) + m.row1.w;
	z = Dot3(m.row2, v) + m.row2.w;
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
	return Dot3(a, b);
}

inline float Len(const Vector3& v)
{
	return sqrt(v * v);
}

inline Vector3 UnitVec(const Vector3& v)
{
	float len = Len(v);
	return v * (1.0f / len);
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
	return Dot4(a, b);
}

inline float Len(const Vector4& v)
{
	return sqrt(v * v);
}

inline Vector4 UnitVec(const Vector4& v)
{
	float len = Len(v);
	return v * (1.0f / len);
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

inline Vector3 Quat::Rotate(const Vector3& v) const
{
	Quat q = *this * Quat(v.x, v.y, v.z, 0.0f) * ~*this;
	return Vector3(q.x, q.y, q.z);
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
	return Quat(a.x + a.x, a.y + b.y, a.z + b.z, a.w + b.w);
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

inline float Len(const Quat& v)
{
	return sqrt(Dot4(v, v));	
}

inline Quat UnitVec(const Quat& v)
{
	float len = Len(v);
	return v * (1.0f / len);	
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

inline void Matrix::MakeProjection(float fovy, float aspect, float near, float far)
{
	float f = 1.0f / tan(fovy / 2.0f);
	row0.Set(f / aspect, 0, 0, 0);
	row1.Set(0, f, 0, 0);
	row2.Set(0, 0, (far + near) / (near - far), (2.0f * far * near) / (near - far));
	row3.Set(0, 0, -1, 0);
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

inline float Matrix::GetElem(int r, int c) const
{
	return ((float*)&row0)[r * 4 + c];
}

inline Vector4 Matrix::GetCol(int c) const
{
	return Vector4(row0[c], row1[c], row2[c], row3[c]);
}

inline Vector3 Mul3x3(const Matrix& m, const Vector3& v)
{
	return Vector3(Dot3(m.row0, v), Dot3(m.row1, v), Dot3(m.row2, v));
}

inline Vector3 Mul3x4(const Matrix& m, const Vector3& v)
{
	return Vector3(Dot3(m.row0, v) + m.row0.w, Dot3(m.row1, v) + m.row1.w, Dot3(m.row2, v) + m.row2.w);
}

inline Vector4 Mul4x4(const Matrix& m, const Vector4& v)
{
	return Vector4(m.row0 * v, m.row1 * v, m.row2 * v, m.row3 * v);
}

inline Matrix Transpose(const Matrix& m)
{
	return Matrix(m.GetCol(0), m.GetCol(1), m.GetCol(2), m.GetCol(3));
}
