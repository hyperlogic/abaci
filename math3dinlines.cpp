#include "math3d.h"
#include <stdio.h>

#ifdef USE_SSE

inline const float& VectorBase::operator[](int i) const 
{ 
	return reinterpret_cast<const float*>(&data)[i]; 
}

inline float& VectorBase::operator[](int i) 
{ 
	return reinterpret_cast<float*>(&data)[i]; 
}

#else

inline const float& VectorBase::operator[](int i) const 
{ 
	return data[i]; 
}

inline float& VectorBase::operator[](int i) 
{ 
	return data[i]; 
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
	return (a.X() * b.X()) + (a.Y() * b.Y());
}

inline float Dot3(const VectorBase& a, const VectorBase& b)
{
	return (a.X() * b.X()) + (a.Y() * b.Y()) + (a.Z() * b.Z());
}

inline float Dot4(const VectorBase& a, const VectorBase& b)
{
#ifdef USE_SSE
	__m128 r = _mm_mul_ps(a.data, b.data);
	r = _mm_hadd_ps(r, r);
	r = _mm_hadd_ps(r, r);
	return *reinterpret_cast<float*>(&r); 
#else
	return (a.X() * b.X()) + (a.Y() * b.Y()) + (a.Z() * b.Z()) + (a.W() * b.W());
#endif
}

inline Vector2::Vector2(float xIn, float yIn) 
{ 
	X() = xIn; 
	Y() = yIn; 
}

inline void Vector2::Set(float xIn, float yIn)
{
	X() = xIn;
	Y() = yIn;
}

inline Vector2& Vector2::operator+=(const Vector2& v)
{
	X() += v.X();
	Y() += v.Y();
	return *this;
}

inline Vector2 operator-(const Vector2& v)
{
	return Vector2(-v.X(), -v.Y());
}

inline Vector2 operator-(const Vector2& a, const Vector2& b)
{
	return Vector2(a.X() - b.X(), a.Y() - b.Y());
}

inline Vector2 operator+(const Vector2& a, const Vector2& b)
{
	return Vector2(a.X() + b.X(), a.Y() + b.Y());
}

inline Vector2 operator*(const Vector2& v, float factor)
{
	return Vector2(factor * v.X(), factor * v.Y());
}

inline Vector2 operator*(float factor, const Vector2& v)
{
	return Vector2(factor * v.X(), factor * v.Y());
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
	Y() = yIn; 
	Z() = zIn; 
}

inline void Vector3::Set(float xIn, float yIn, float zIn)
{
	X() = xIn;
	Y() = yIn;
	Z() = zIn;
}

inline void Vector3::Normalize()
{
	float invLen = 1.0f / sqrt((X() * X()) + (Y() * Y()) + (Z() * Z()));
 	X() *= invLen;
	Y() *= invLen;
	Z() *= invLen;
}

inline void Vector3::Mul3x4(const Matrix& m)
{
	Vector3 v(*this);
	X() = Dot3(m.row0, v) + m.row0.W();
	Y() = Dot3(m.row1, v) + m.row1.W();
	Z() = Dot3(m.row2, v) + m.row2.W();
}

inline Vector3 operator-(const Vector3& v)
{
	return Vector3(-v.X(), -v.Y(), -v.Z());
}

inline Vector3 operator-(const Vector3& a, const Vector3& b)
{
	return Vector3(a.X() - b.X(), a.Y() - b.Y(), a.Z() - b.Z());
}

inline Vector3 operator+(const Vector3& a, const Vector3& b)
{
#ifdef USE_SSE
	return Vector3(_mm_add_ps(a.data, b.data));
#else
	return Vector3(a.X() + b.X(), a.Y() + b.Y(), a.Z() + b.Z());
#endif
}

inline Vector3 operator*(const Vector3& v, float factor)
{
	return Vector3(v.X() * factor, v.Y() * factor, v.Z() * factor);
}

inline Vector3 operator*(float factor, const Vector3& v)
{
	return Vector3(v.X() * factor, v.Y() * factor, v.Z() * factor);
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
	return Vector3((a.Y() * b.Z()) - (a.Z() * b.Y()), 
	               (a.Z() * b.X()) - (a.X() * b.Z()), 
	               (a.X() * b.Y()) - (a.Y() * b.X()));
}

inline Vector4::Vector4(const Vector3& v, float wIn) 
{ 
	X() = v.X(); 
	Y() = v.Y(); 
	Z() = v.Z(); 
    W() = wIn; 
}

inline Vector4::Vector4(float xIn, float yIn, float zIn, float wIn) 
{ 
	X() = xIn; 
	Y() = yIn; 
	Z() = zIn; 
	W() = wIn; 
}

inline void Vector4::Set(float xIn, float yIn, float zIn, float wIn)
{
	X() = xIn;
	Y() = yIn;
	Z() = zIn;
	W() = wIn;
}

inline Vector4 operator-(const Vector4& v)
{
	return Vector4(-v.X(), -v.Y(), -v.Z(), -v.W());
}

inline Vector4 operator-(const Vector4& a, const Vector4& b)
{
	return Vector4(a.X() - b.X(), a.Y() - b.Y(), a.Z() - b.Z(), a.W() - b.W());
}

inline Vector4 operator+(const Vector4& a, const Vector4& b)
{
#ifdef USE_SSE
	return Vector4(_mm_add_ps(a.data, b.data));
#else
	return Vector4(a.X() + b.X(), a.Y() + b.Y(), a.Z() + b.Z(), a.W() + b.W());
#endif
}

inline Vector4 operator*(const Vector4& v, float factor)
{
	return Vector4(factor * v.X(), factor * v.Y(), factor * v.Z(), factor * v.W());
}

inline Vector4 operator*(float factor, const Vector4& v)
{
	return Vector4(factor * v.X(), factor * v.Y(), factor * v.Z(), factor * v.W());	
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
	X() = xIn; 
	Y() = yIn; 
	Z() = zIn; 
	W() = wIn; 
}

inline Quat::Quat(const Vector3& axis, float angle)
{
	Vector3 n = UnitVec(axis) * sin(angle/2.0f);
	X() = n.X();
	Y() = n.Y();
	Z() = n.Z();
	W() = cos(angle/2.0f);
}

inline Vector3 Quat::Rotate(const Vector3& v) const
{
	Quat q = *this * Quat(v.X(), v.Y(), v.Z(), 0.0f) * ~*this;
	return Vector3(q.X(), q.Y(), q.Z());
}

inline Quat operator~(const Quat& v)
{
	return Quat(-v.X(), -v.Y(), -v.Z(), v.W());
}

inline Quat operator-(const Quat& v)
{
	return Quat(-v.X(), -v.Y(), -v.Z(), -v.W());
}

inline Quat operator-(const Quat& a, const Quat& b)
{
	return Quat(a.X() - b.X(), a.Y() - b.Y(), a.Z() - b.Z(), a.W() - b.W());
}

inline Quat operator+(const Quat& a, const Quat& b)
{
	return Quat(a.X() + a.X(), a.Y() + b.Y(), a.Z() + b.Z(), a.W() + b.W());
}

inline Quat operator*(const Quat& v, float factor)
{
	return Quat(factor * v.X(), factor * v.Y(), factor * v.Z(), factor * v.W());
}

inline Quat operator*(float factor, const Quat& v)
{
	return Quat(factor * v.X(), factor * v.Y(), factor * v.Z(), factor * v.W());
}

inline Quat operator*(const Quat& q1, const Quat& q2)
{
	return Quat( q1.X() * q2.W() + q1.Y() * q2.Z() - q1.Z() * q2.Y() + q1.W() * q2.X(),
				-q1.X() * q2.Z() + q1.Y() * q2.W() + q1.Z() * q2.X() + q1.W() * q2.Y(),
				 q1.X() * q2.Y() - q1.Y() * q2.X() + q1.Z() * q2.W() + q1.W() * q2.Z(),
				-q1.X() * q2.X() - q1.Y() * q2.Y() - q1.Z() * q2.Z() + q1.W() * q2.W());
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
	return Vector3(row0.X(), row1.X(), row2.X());
}

inline Vector3 Matrix::GetYAxis() const
{
	return Vector3(row0.Y(), row1.Y(), row2.Y());
}

inline Vector3 Matrix::GetZAxis() const
{
	return Vector3(row0.Z(), row1.Z(), row2.Z());
}

inline Vector3 Matrix::GetTrans() const
{
	return Vector3(row0.W(), row1.W(), row2.W());
}

inline void Matrix::SetXAxis(const Vector3 xAxis)
{
	row0.X() = xAxis.X();
	row1.X() = xAxis.Y();
	row2.X() = xAxis.Z();
	row3.X() = 0.0f;
}

inline void Matrix::SetYAxis(const Vector3 yAxis)
{
	row0.Y() = yAxis.X();
	row1.Y() = yAxis.Y();
	row2.Y() = yAxis.Z();
	row3.Y() = 0.0f;
}

inline void Matrix::SetZAxis(const Vector3 zAxis)
{
	row0.Z() = zAxis.X();
	row1.Z() = zAxis.Y();
	row2.Z() = zAxis.Z();
	row3.Z() = 0.0f;
}

inline void Matrix::SetTrans(const Vector3 trans)
{
	row0.W() = trans.X();
	row1.W() = trans.Y();
	row2.W() = trans.Z();
	row3.W() = 1.0f;
}

inline void Matrix::SetScale(float scale)
{
	SetXAxis(GetXAxis() * scale);
	SetYAxis(GetYAxis() * scale);
	SetZAxis(GetZAxis() * scale);	
}

inline void Matrix::SetScale(const Vector3 scale)
{
	SetXAxis(GetXAxis() * scale.X());
	SetYAxis(GetYAxis() * scale.Y());
	SetZAxis(GetZAxis() * scale.Z());
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
	return Vector3(Dot3(m.row0, v) + m.row0.W(), Dot3(m.row1, v) + m.row1.W(), Dot3(m.row2, v) + m.row2.W());
}

inline Vector4 Mul4x4(const Matrix& m, const Vector4& v)
{
	return Vector4(m.row0 * v, m.row1 * v, m.row2 * v, m.row3 * v);
}

inline Matrix Transpose(const Matrix& m)
{
	return Matrix(m.GetCol(0), m.GetCol(1), m.GetCol(2), m.GetCol(3));
}
