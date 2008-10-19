#include "math3d.h"
#include <stdio.h>

Quat QuatExp(const Quat& q)
{
	float angle = Len(Vector3(q.X(), q.Y(), q.Z()));
	Vector3 n = UnitVec(Vector3(q.X(), q.Y(), q.Z())) * sin(angle/2.0f);
	return Quat(n.X(), n.Y(), n.Z(), cos(angle/2.0f));
}

Quat QuatLog(const Quat& q)
{
	float cos_a = q.W();
	if (cos_a > 1.0f) cos_a = 1.0f;
	if (cos_a < -1.0f) cos_a = -1.0f;

    float sin_a = (float)sqrt(1.0f - cos_a * cos_a);

    if (fabs(sin_a) < 0.0005f)
		sin_a = 1.0f;
	else
		sin_a = 1.f/sin_a;

    float angle = 2.0f * (float)acos(cos_a);
	Quat log;
    log.X() = q.X() * sin_a * angle;
    log.Y() = q.Y() * sin_a * angle;
    log.Z() = q.Z() * sin_a * angle;
	log.W() = 0.0f;
	return log;
}

Matrix::Matrix(const Quat& q)
{
	SetXAxis(Rotate(q, Vector3(1.0f,0.0f,0.0f)));
	SetYAxis(Rotate(q, Vector3(0.0f,1.0f,0.0f)));
	SetZAxis(Rotate(q, Vector3(0.0f,0.0f,1.0f)));
	SetTrans(Vector3(0.0f,0.0f,0.0f));
}

Matrix::Matrix(const Quat& q, const Vector3& trans)
{
	SetXAxis(Rotate(q, Vector3(1.0f,0.0f,0.0f)));
	SetYAxis(Rotate(q, Vector3(0.0f,1.0f,0.0f)));
	SetZAxis(Rotate(q, Vector3(0.0f,0.0f,1.0f)));
	SetTrans(trans);
}

Matrix::Matrix(const Vector3& scale, const Quat& q, const Vector3& trans)
{
	SetXAxis(Rotate(q, Vector3(scale.X(),0.0f,0.0f)));
	SetYAxis(Rotate(q, Vector3(0.0f,scale.Y(),0.0f)));
	SetZAxis(Rotate(q, Vector3(0.0f,0.0f,scale.Z())));
	SetTrans(trans);
}

Matrix::Matrix(const Vector3 axis, float angle)
{
	Quat q(axis, angle);
	SetXAxis(Rotate(q, Vector3(1.0f,0.0f,0.0f)));
	SetYAxis(Rotate(q, Vector3(0.0f,1.0f,0.0f)));
	SetZAxis(Rotate(q, Vector3(0.0f,0.0f,1.0f)));
	SetTrans(Vector3(0.0f,0.0f,0.0f));	
}

Matrix operator*(const Matrix& a, const Matrix& b)
{
	Matrix bt = Transpose(b);
	return Matrix( Vector4(a.row0 * bt.row0, a.row0 * bt.row1, a.row0 * bt.row2, a.row0 * bt.row3), 
				   Vector4(a.row1 * bt.row0, a.row1 * bt.row1, a.row1 * bt.row2, a.row1 * bt.row3),
				   Vector4(a.row2 * bt.row0, a.row2 * bt.row1, a.row2 * bt.row2, a.row2 * bt.row3),
				   Vector4(a.row3 * bt.row0, a.row3 * bt.row1, a.row3 * bt.row2, a.row3 * bt.row3));
}

Matrix OrthonormalInverse(const Matrix& m)
{
	Matrix r(m);
	r.SetTrans(Vector3(0, 0, 0));
	r = Transpose(r);
	r.SetTrans(-Transform3x4(r, m.GetTrans()));
	return r;
}
