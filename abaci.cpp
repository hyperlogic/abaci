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

#include "abaci.h"
#include <stdio.h>

#ifdef ABACI_NAMESPACE
namespace ABACI_NAMESPACE
{
#endif

Quat QuatExp(const Quat& q)
{
	float angle = Len(Vector3(q.X(), q.Y(), q.Z()));
	Vector3 n;
	if (angle > 0.0001f)
		n = UnitVec(Vector3(q.X(), q.Y(), q.Z())) * sin(angle/2.0f);
	else
		n.Set(0,0,0);

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

Matrix operator+(const Matrix& a, const Matrix& b)
{
	return Matrix(a.row0 + b.row0, a.row1 + b.row1, a.row2 + b.row2, a.row3 + b.row2);
}

Matrix operator-(const Matrix& a, const Matrix& b)
{
	return Matrix(a.row0 - b.row0, a.row1 - b.row1, a.row2 - b.row2, a.row3 - b.row2);
}

Matrix operator*(const Matrix& a, const Matrix& b)
{
	// TODO: slow.
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

Quat Matrix::GetQuat() const
{
	// TODO: add to unit test!

	int x = 0, y = 1, z = 2;
	// create a matrix with no scale
	Matrix m(UnitVec(GetXAxis()), UnitVec(GetYAxis()), UnitVec(GetZAxis()), Vector3(0,0,0));
	float trace = m.Elem(0,0) + m.Elem(1,1) + m.Elem(2,2);
	Quat q;
	if (trace > -1.0f)
	{
		int i = x, j = y, k = z;
		if (m.Elem(y,y) > m.Elem(x,x)) { i = y; j = z; k = x; }
		if (m.Elem(z,z) > m.Elem(i,i)) { i = z; j = x; k = y; }
		float r = sqrt(m.Elem(i,i) - m.Elem(j,j) - m.Elem(k,k) + 1.0f);
		q[i] = r / 2.0f;
		q[j] = (m.Elem(i,j) + m.Elem(j,i)) / (2.0f * r);
		q[k] = (m.Elem(k,i) + m.Elem(i,k)) / (2.0f * r);
		q[3] = (m.Elem(k,j) + m.Elem(j,k)) / (2.0f * r);
	}
    else
	{
		q.W() = sqrt(1.0f + trace) / 2.0f;
		q.X() = (row2[1] - row1[2]) / (4.0f * q.W());
		q.Y() = (row0[2] - row2[0]) / (4.0f * q.W());
		q.Z() = (row1[0] - row0[1]) / (4.0f * q.W());
	}
	return q;
}

#ifdef ABACI_NAMESPACE
}  // namespace
#endif
