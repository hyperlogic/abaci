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

// Quaternian exponential
Quat Quat::Exp() const
{
	float angle = Vector3(i, j, k).Len();
	Vector3 n;
	if (angle > 0.0001f)
		n = Vector3(i, j, k).Unit() * sin(angle / 2.0f);
	else
		n.Set(0,0,0);

	return Quat(n.x, n.y, n.z, cos(angle / 2.0f));
}

// Quaternian logarithm
Quat Quat::Log() const
{
	float cos_a = r;
	if (cos_a > 1.0f) cos_a = 1.0f;
	if (cos_a < -1.0f) cos_a = -1.0f;

    float sin_a = (float)sqrt(1.0f - cos_a * cos_a);

    if (fabs(sin_a) < 0.0005f)
		sin_a = 1.0f;
	else
		sin_a = 1.f/sin_a;

    float angle = 2.0f * (float)acos(cos_a);
	Quat log;
    log.i = i * sin_a * angle;
    log.j = j * sin_a * angle;
    log.k = k * sin_a * angle;
	log.r = 0.0f;
	return log;
}

// Create a Matrix from a Quat and a translation.
Matrix Matrix::QuatTrans(const Quat& q, const Vector3& trans)
{
	Matrix m;
	m.SetXAxis(q.Rotate(Vector3(1.0f,0.0f,0.0f)));
	m.SetYAxis(q.Rotate(Vector3(0.0f,1.0f,0.0f)));
	m.SetZAxis(q.Rotate(Vector3(0.0f,0.0f,1.0f)));
	m.SetTrans(Vector3(0.0f,0.0f,0.0f));
	return m;
}

// Create a Matrix from a Scale vector, a Quat and a translation.
Matrix Matrix::ScaleQuatTrans(const Vector3& scale, const Quat& q, const Vector3& trans)
{
	Matrix m;
	m.SetXAxis(q.Rotate(Vector3(scale.x,0.0f,0.0f)));
	m.SetYAxis(q.Rotate(Vector3(0.0f,scale.y,0.0f)));
	m.SetZAxis(q.Rotate(Vector3(0.0f,0.0f,scale.z)));
	m.SetTrans(trans);
	return m;
}

// Create a Matrix from a rotation represented by an axis and an angle.
Matrix Matrix::AxisAngle(const Vector3 axis, float angle)
{
	Matrix m;
	Quat q = Quat::AxisAngle(axis, angle);
	m.SetXAxis(q.Rotate(Vector3(1.0f,0.0f,0.0f)));
	m.SetYAxis(q.Rotate(Vector3(0.0f,1.0f,0.0f)));
	m.SetZAxis(q.Rotate(Vector3(0.0f,0.0f,1.0f)));
	m.SetTrans(Vector3(0.0f,0.0f,0.0f));
	return m;
}

// Returns the rotation component of this Matrix.
Quat Matrix::GetQuat() const
{
	// TODO: add to unit test!

	int x = 0, y = 1, z = 2;
	// create a matrix with no scale
	Matrix m = Axes(GetXAxis().Unit(), GetYAxis().Unit(), GetZAxis().Unit(), Vector3(0,0,0));
	float trace = m.Elem(0,0) + m.Elem(1,1) + m.Elem(2,2);
	Vector4 q;
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
		q.w = sqrt(1.0f + trace) / 2.0f;
		q.x = (row2[1] - row1[2]) / (4.0f * q.w);
		q.y = (row0[2] - row2[0]) / (4.0f * q.w);
		q.z = (row1[0] - row0[1]) / (4.0f * q.w);
	}
	return Quat(q.x, q.y, q.z, q.w);
}

// Matrix addition
Matrix operator+(const Matrix& a, const Matrix& b)
{
	return Matrix::Rows(a.row0 + b.row0, a.row1 + b.row1, a.row2 + b.row2, a.row3 + b.row2);
}

// Matrix subtraction
Matrix operator-(const Matrix& a, const Matrix& b)
{
	return Matrix::Rows(a.row0 - b.row0, a.row1 - b.row1, a.row2 - b.row2, a.row3 - b.row2);
}

// Matrix multiplication
Matrix operator*(const Matrix& a, const Matrix& b)
{
	Matrix bt = b.Transpose();
	return Matrix::Rows( Vector4(Dot(a.row0, bt.row0), Dot(a.row0, bt.row1), Dot(a.row0, bt.row2), Dot(a.row0, bt.row3)), 
						 Vector4(Dot(a.row1, bt.row0), Dot(a.row1, bt.row1), Dot(a.row1, bt.row2), Dot(a.row1, bt.row3)),
						 Vector4(Dot(a.row2, bt.row0), Dot(a.row2, bt.row1), Dot(a.row2, bt.row2), Dot(a.row2, bt.row3)),
						 Vector4(Dot(a.row3, bt.row0), Dot(a.row3, bt.row1), Dot(a.row3, bt.row2), Dot(a.row3, bt.row3)));
}

// If the 3x3 portion of this Matrix is Orthogonal (i.e. columns are orthogonal unit vectors)
// this will return the Inverse of that matrix.
Matrix Matrix::OrthoInverse() const
{
	Matrix r(*this);
	r.SetTrans(Vector3(0, 0, 0));
	r = r.Transpose();
	r.SetTrans(-r.Mul3x4(GetTrans()));
	return r;
}

// Full 4x4 Matrix Inverse, returns Identity if matrix has no inverse.
Matrix Matrix::FullInverse() const
{
	Matrix m;
	if (::FullInverse((*this), m))
		return m;
	else
		return Identity();
}

template <typename T>
void Swap(T& a, T& b)
{
	T temp = a;
	a = b;
	b = temp;
}

// Full 4x4 Matrix inverse, returns false if Matrix has no inverse.
bool FullInverse(const Matrix& m, Matrix& result)
{
    // Gaussian-Jordan Elimination, TODO: There are more numerically stable ways to do this...

	float temp[4][8];
	float* row[4];

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
	temp[0][4] = 1.0f; temp[0][5] = 0.0f; temp[0][6] = 0.0f; temp[0][7] = 0.0f;
	temp[1][4] = 0.0f; temp[1][5] = 1.0f; temp[1][6] = 0.0f; temp[1][7] = 0.0f;
	temp[2][4] = 0.0f; temp[2][5] = 0.0f; temp[2][6] = 1.0f; temp[2][7] = 0.0f;
	temp[3][4] = 0.0f; temp[3][5] = 0.0f; temp[3][6] = 0.0f; temp[3][7] = 1.0f;

	// bubble up row with largest leading number (partial pivot)
	if (fabs(row[0][0]) < fabs(row[1][0]))
		Swap(row[0], row[1]);
	if (fabs(row[0][0]) < fabs(row[2][0]))
		Swap(row[0], row[2]);
	if (fabs(row[0][0]) < fabs(row[3][0]))
		Swap(row[0], row[3]);

	if (fabs(row[0][0]) < 0.00001f)  // column is all zeros, there is no inverse.
		return false;

	// mult row[0] by 1/row[0][0].  To introduce a leading 1.
	float s = 1.0f / row[0][0];
	for (int c = 0; c < 8; ++c)	row[0][c] *= s;

	// add multiples of top row to lower rows so that all entries below leading 1 become zeros.
	for (int r = 1; r < 4; ++r)
	{
		float s = row[r][0];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[0][c];
	}

	// move row with largest leading number 
	if (fabs(row[1][1]) < fabs(row[2][1]))
		Swap(row[1], row[2]);
	if (fabs(row[1][1]) < fabs(row[3][1]))
		Swap(row[1], row[3]);

	if (fabs(row[1][1]) < 0.00001f)  // column is all zeros, there is no inverse.
		return false;

	// mult row[1] by 1/row[1][1].  To introduce a leading 1.
	s = 1.0f / row[1][1];
	for (int c = 0; c < 8; ++c)	row[1][c] *= s;
	
	// add multiples of top row to lower rows so that all entries below leading 1 become zeros.
	for (int r = 2; r < 4; ++r)
	{
		float s = row[r][1];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[1][c];
	}

	// move row with largest leading number 
	if (fabs(row[2][2]) < fabs(row[3][2]))
		Swap(row[2], row[3]);

	if (fabs(row[2][2]) < 0.00001f)  // column is all zeros, there is no inverse.
		return false;

	// mult row[2] by 1/row[2][2].  To introduce a leading 1.
	s = 1.0f / row[2][2];
	for (int c = 0; c < 8; ++c)	row[2][c] *= s;
	
	// add multiples of top row to lower rows so that all entries below leading 1 become zeros.
	for (int r = 3; r < 4; ++r)
	{
		float s = row[r][2];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[2][c];
	}

	// mult row[3] by 1/row[3][3]. To introduce a leading 1.
	s = 1.0f / row[3][3];
	for (int c = 0; c < 8; ++c)	row[3][c] *= s;

	// at this point row matrix should be in row-echelon form.

	// add multiples of row[3] to above rows, to zero out that column
	for (int r = 0; r < 3; ++r)
	{
		float s = row[r][3];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[3][c];
	}

	// add multiples of row[2] to above rows, to zero out that column
	for (int r = 0; r < 2; ++r)
	{
		float s = row[r][2];
		for (int c = 0; c < 8; ++c) row[r][c] -= s * row[2][c];
	}

	// add multiples of row[1] to above row, to zero out that column
	for (int r = 0; r < 1; ++r)
	{
		float s = row[r][1];
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
void PrintMatrix(const Matrix& m)
{
	printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.row0.x, m.row0.y, m.row0.z, m.row0.w);
	printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.row1.x, m.row1.y, m.row1.z, m.row1.w);
	printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.row2.x, m.row2.y, m.row2.z, m.row2.w);
	printf("| %15.5f, %15.5f, %15.5f, %15.5f |\n", m.row3.x, m.row3.y, m.row3.z, m.row3.w);
}

#ifdef ABACI_NAMESPACE
}  // namespace
#endif
