#include "math3d.h"
#include <stdio.h>

float DegToRad(float deg)
{ 
	return deg * ((2 * PI) / 180.0f);
}

float RadToDeg(float rad)
{
	return rad * (180.0f / (2 * PI));
}

float Dot2(const Vector& a, const Vector& b)
{
	return (a.x * b.x) + (a.y * b.y);
}

float Dot3(const Vector& a, const Vector& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z);
}

float Dot4(const Vector& a, const Vector& b)
{
	return (a.x * b.x) + (a.y * b.y) + (a.z * b.z) + (a.w * b.w);
}

void Vector2::Set(float xIn, float yIn)
{
	x = xIn;
	y = yIn;
}

Vector2& Vector2::operator+=(const Vector2& v)
{
	x += v.x;
	y += v.y;
	return *this;
}

Vector2 operator-(const Vector2& v)
{
	return Vector2(-v.x, -v.y);
}

Vector2 operator-(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x - b.x, a.y - b.y);
}

Vector2 operator+(const Vector2& a, const Vector2& b)
{
	return Vector2(a.x + b.x, a.y + b.y);
}

Vector2 operator*(const Vector2& v, float factor)
{
	return Vector2(factor * v.x, factor * v.y);
}

Vector2 operator*(float factor, const Vector2& v)
{
	return Vector2(factor * v.x, factor * v.y);
}

float operator*(const Vector2& a, const Vector2& b)
{
	return Dot2(a, b);
}

float Len(const Vector2& v)
{
	return sqrt(v * v);
}

Vector2 UnitVec(const Vector2& v)
{
	float len = Len(v);
	return v * (1.0f / len);
}

void Vector3::Set(float xIn, float yIn, float zIn)
{
	x = xIn;
	y = yIn;
	z = zIn;
}

void Vector3::Normalize()
{
	float invLen = 1.0f / sqrt((x * x) + (y * y) + (z * z));
 	x *= invLen;
	y *= invLen;
	z *= invLen;
}

void Vector3::Mul3x4(const Matrix& m)
{
	Vector3 v(*this);
	x = Dot3(m.row0, v) + m.row0.w;
	y = Dot3(m.row1, v) + m.row1.w;
	z = Dot3(m.row2, v) + m.row2.w;
}

Vector3 operator-(const Vector3& v)
{
	return Vector3(-v.x, -v.y, -v.z);
}

Vector3 operator-(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x - b.x, a.y - b.y, a.z - b.z);
}

Vector3 operator+(const Vector3& a, const Vector3& b)
{
	return Vector3(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector3 operator*(const Vector3& v, float factor)
{
	return Vector3(v.x * factor, v.y * factor, v.z * factor);
}

Vector3 operator*(float factor, const Vector3& v)
{
	return Vector3(v.x * factor, v.y * factor, v.z * factor);
}

float operator*(const Vector3& a, const Vector3& b)
{
	return Dot3(a, b);
}

float Len(const Vector3& v)
{
	return sqrt(v * v);
}

Vector3 UnitVec(const Vector3& v)
{
	float len = Len(v);
	return v * (1.0f / len);
}

Vector3 operator%(const Vector3& a, const Vector3& b)
{
	return Vector3((a.y * b.z) - (a.z * b.y), 
	               (a.z * b.x) - (a.x * b.z), 
	               (a.x * b.y) - (a.y * b.x));
}

void Vector4::Set(float xIn, float yIn, float zIn, float wIn)
{
	x = xIn;
	y = yIn;
	z = zIn;
	w = wIn;
}

Vector4 operator-(const Vector4& v)
{
	return Vector4(-v.x, -v.y, -v.z, -v.w);
}

Vector4 operator-(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

Vector4 operator+(const Vector4& a, const Vector4& b)
{
	return Vector4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

Vector4 operator*(const Vector4& v, float factor)
{
	return Vector4(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

Vector4 operator*(float factor, const Vector4& v)
{
	return Vector4(factor * v.x, factor * v.y, factor * v.z, factor * v.w);	
}

float operator*(const Vector4& a, const Vector4& b)
{
	return Dot4(a, b);
}

float Len(const Vector4& v)
{
	return sqrt(v * v);
}

Vector4 UnitVec(const Vector4& v)
{
	float len = Len(v);
	return v * (1.0f / len);
}

Quat::Quat(const Vector3& axis, float angle)
{
	Vector3 n = UnitVec(axis) * sin(angle/2.0f);
	x = n.x;
	y = n.y;
	z = n.z;
	w = cos(angle/2.0f);
}

Vector3 Quat::Rotate(const Vector3& v) const
{
	Quat q = *this * Quat(v.x, v.y, v.z, 0.0f) * ~*this;
	return Vector3(q.x, q.y, q.z);
}

Quat operator~(const Quat& v)
{
	return Quat(-v.x, -v.y, -v.z, v.w);
}

Quat operator-(const Quat& v)
{
	return Quat(-v.x, -v.y, -v.z, -v.w);
}

Quat operator-(const Quat& a, const Quat& b)
{
	return Quat(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

Quat operator+(const Quat& a, const Quat& b)
{
	return Quat(a.x + a.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

Quat operator*(const Quat& v, float factor)
{
	return Quat(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

Quat operator*(float factor, const Quat& v)
{
	return Quat(factor * v.x, factor * v.y, factor * v.z, factor * v.w);
}

Quat operator*(const Quat& q1, const Quat& q2)
{
	return Quat( q1.x * q2.w + q1.y * q2.z - q1.z * q2.y + q1.w * q2.x,
				-q1.x * q2.z + q1.y * q2.w + q1.z * q2.x + q1.w * q2.y,
				 q1.x * q2.y - q1.y * q2.x + q1.z * q2.w + q1.w * q2.z,
				-q1.x * q2.x - q1.y * q2.y - q1.z * q2.z + q1.w * q2.w);
}

float Len(const Quat& v)
{
	return sqrt(Dot4(v, v));	
}

Quat UnitVec(const Quat& v)
{
	float len = Len(v);
	return v * (1.0f / len);	
}

Quat QuatExp(const Quat& q)
{
	float angle = Len(q);
	Vector3 n = UnitVec(Vector3(q.x, q.y, q.z)) * sin(angle/2.0f);
	return Quat(n.x, n.y, n.z, cos(angle/2.0f));
}

Quat QuatLog(const Quat& q)
{
	float cos_a = q.w;
	if (cos_a > 1.0f) cos_a = 1.0f;
	if (cos_a < -1.0f) cos_a = -1.0f;

    float sin_a = (float)sqrt(1.0f - cos_a * cos_a);

    if (fabs(sin_a) < 0.0005f)
		sin_a = 1.0f;
	else
		sin_a = 1.f/sin_a;

    float angle = 2.0f * (float)acos(cos_a);
	Quat log;
    log.x = q.x * sin_a * angle;
    log.y = q.y * sin_a * angle;
    log.z = q.z * sin_a * angle;
	log.w = 0.0f;
	return log;
}

Matrix::Matrix(const Vector3& xAxis, const Vector3& yAxis, const Vector3& zAxis, const Vector3& trans)
{
	SetXAxis(xAxis);
	SetYAxis(yAxis);
	SetZAxis(zAxis);
	SetTrans(trans);
}

Matrix::Matrix(const Vector4& row0In, const Vector4& row1In, const Vector4& row2In, const Vector4& row3In)
{
	row0 = row0In;
	row1 = row1In;
	row2 = row2In;
	row3 = row3In;
}

Matrix::Matrix(const Quat& q)
{
	SetXAxis(q.Rotate(Vector3(1.0f,0.0f,0.0f)));
	SetYAxis(q.Rotate(Vector3(0.0f,1.0f,0.0f)));
	SetZAxis(q.Rotate(Vector3(0.0f,0.0f,1.0f)));
	SetTrans(Vector3(0.0f,0.0f,0.0f));
}

Matrix::Matrix(const Quat& q, const Vector3& trans)
{
	SetXAxis(q.Rotate(Vector3(1.0f,0.0f,0.0f)));
	SetYAxis(q.Rotate(Vector3(0.0f,1.0f,0.0f)));
	SetZAxis(q.Rotate(Vector3(0.0f,0.0f,1.0f)));
	SetTrans(trans);	
}

Matrix::Matrix(const Vector3& scale, const Quat& q, const Vector3& trans)
{
	SetXAxis(q.Rotate(Vector3(scale.x,0.0f,0.0f)));
	SetYAxis(q.Rotate(Vector3(0.0f,scale.y,0.0f)));
	SetZAxis(q.Rotate(Vector3(0.0f,0.0f,scale.z)));
	SetTrans(trans);
}

Matrix::Matrix(const Vector3 axis, float angle)
{
	Quat q(axis, angle);
	SetXAxis(q.Rotate(Vector3(1.0f,0.0f,0.0f)));
	SetYAxis(q.Rotate(Vector3(0.0f,1.0f,0.0f)));
	SetZAxis(q.Rotate(Vector3(0.0f,0.0f,1.0f)));
	SetTrans(Vector3(0.0f,0.0f,0.0f));	
}

void Matrix::MakeIdent()
{
	row0.Set(1.0f, 0.0f, 0.0f, 0.0f);
	row1.Set(0.0f, 1.0f, 0.0f, 0.0f);
	row2.Set(0.0f, 0.0f, 1.0f, 0.0f);
	row3.Set(0.0f, 0.0f, 0.0f, 1.0f);
}

void Matrix::MakeProjection(float fovy, float aspect, float near, float far)
{
	float f = 1.0f / tan(fovy / 2.0f);
	row0.Set(f / aspect, 0, 0, 0);
	row1.Set(0, f, 0, 0);
	row2.Set(0, 0, (far + near) / (near - far), (2.0f * far * near) / (near - far));
	row3.Set(0, 0, -1, 0);
}

void Matrix::MakeLookAt(const Vector3& eye, const Vector3& target, const Vector3& up)
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

Vector3 Matrix::GetXAxis() const
{
	return Vector3(row0.x, row1.x, row2.x);
}

Vector3 Matrix::GetYAxis() const
{
	return Vector3(row0.y, row1.y, row2.y);
}

Vector3 Matrix::GetZAxis() const
{
	return Vector3(row0.z, row1.z, row2.z);
}

Vector3 Matrix::GetTrans() const
{
	return Vector3(row0.w, row1.w, row2.w);
}

void Matrix::SetXAxis(const Vector3 xAxis)
{
	row0.x = xAxis.x;
	row1.x = xAxis.y;
	row2.x = xAxis.z;
	row3.x = 0.0f;
}

void Matrix::SetYAxis(const Vector3 yAxis)
{
	row0.y = yAxis.x;
	row1.y = yAxis.y;
	row2.y = yAxis.z;
	row3.y = 0.0f;
}

void Matrix::SetZAxis(const Vector3 zAxis)
{
	row0.z = zAxis.x;
	row1.z = zAxis.y;
	row2.z = zAxis.z;
	row3.z = 0.0f;
}

void Matrix::SetTrans(const Vector3 trans)
{
	row0.w = trans.x;
	row1.w = trans.y;
	row2.w = trans.z;
	row3.w = 1.0f;
}

void Matrix::SetScale(float scale)
{
	SetXAxis(GetXAxis() * scale);
	SetYAxis(GetYAxis() * scale);
	SetZAxis(GetZAxis() * scale);	
}

void Matrix::SetScale(const Vector3 scale)
{
	SetXAxis(GetXAxis() * scale.x);
	SetYAxis(GetYAxis() * scale.y);
	SetZAxis(GetZAxis() * scale.z);
}

float Matrix::GetElem(int r, int c) const
{
	return ((float*)&row0)[r * 4 + c];
}

Vector4 Matrix::GetCol(int c) const
{
	return Vector4(row0[c], row1[c], row2[c], row3[c]);
}

Matrix operator*(const Matrix& a, const Matrix& b)
{
	Matrix bt = Transpose(b);
	return Matrix( Vector4(a.row0 * bt.row0, a.row0 * bt.row1, a.row0 * bt.row2, a.row0 * bt.row3), 
				   Vector4(a.row1 * bt.row0, a.row1 * bt.row1, a.row1 * bt.row2, a.row1 * bt.row3),
				   Vector4(a.row2 * bt.row0, a.row2 * bt.row1, a.row2 * bt.row2, a.row2 * bt.row3),
				   Vector4(a.row3 * bt.row0, a.row3 * bt.row1, a.row3 * bt.row2, a.row3 * bt.row3));
}

Vector3 Mul3x3(const Matrix& m, const Vector3& v)
{
	return Vector3(Dot3(m.row0, v), Dot3(m.row1, v), Dot3(m.row2, v));
}

Vector3 Mul3x4(const Matrix& m, const Vector3& v)
{
	return Vector3(Dot3(m.row0, v) + m.row0.w, Dot3(m.row1, v) + m.row1.w, Dot3(m.row2, v) + m.row2.w);
}

Vector4 Mul4x4(const Matrix& m, const Vector4& v)
{
	return Vector4(m.row0 * v, m.row1 * v, m.row2 * v, m.row3 * v);
}

Matrix Transpose(const Matrix& m)
{
	return Matrix(m.GetCol(0), m.GetCol(1), m.GetCol(2), m.GetCol(3));
}

Matrix OrthonormalInverse(const Matrix& m)
{
	Matrix r(m);
	r.SetTrans(Vector3(0, 0, 0));
	r = Transpose(r);
	r.SetTrans(-Mul3x4(r, m.GetTrans()));
	return r;
}
