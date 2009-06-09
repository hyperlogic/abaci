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
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "unittest.h"

#ifdef ABACI_NAMESPACE
using namespace ABACI_NAMESPACE;
#endif

// static test containers
std::vector<float> s_floatVec;
std::vector<Vector2> s_vector2Vec;
std::vector<Vector3> s_vector3Vec;
std::vector<Vector4> s_vector4Vec;
std::vector<Quat> s_quatVec;
std::vector<Matrix> s_matrixVec;

bool FuzzyEqual(float rhs, float lhs)
{
	const float epsilon = 0.0001f;
	return fabs(rhs - lhs) <= epsilon;
}

bool FloatTest(float rhs, float lhs)
{
	// For testing purposes, two Nans are "equal" even if they have different bits.
	if (isnan(rhs) && isnan(lhs))
		return true;

	return rhs == lhs;
}

// 0..1
float RandomNormalizedFloat()
{
	return (float)rand() / RAND_MAX;
}

// -range to range
float RandomFloat(float range = 100.0f)
{
	return (((float)rand() / (RAND_MAX / 2)) - 1.0f) * range;
}

void InitTestData()
{
	srand((unsigned int)666);

	// add random values
	const int RANDOM_DATA_SIZE = 64;
	for (int i = 0; i < RANDOM_DATA_SIZE; ++i)
	{
		s_floatVec.push_back(RandomFloat());
		s_vector2Vec.push_back(Vector2(RandomFloat(), RandomFloat()));
		s_vector3Vec.push_back(Vector3(RandomFloat(), RandomFloat(), RandomFloat()));
		s_vector4Vec.push_back(Vector4(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat()));

		Vector4 row0(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
		Vector4 row1(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
		Vector4 row2(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
		Vector4 row3(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
		s_matrixVec.push_back(Matrix(row0, row1, row2, row3));
	}

	// zero
	s_floatVec.push_back(0);
	s_vector2Vec.push_back(Vector2(0, 0));
	s_vector3Vec.push_back(Vector3(0, 0, 0));
	s_vector4Vec.push_back(Vector4(0, 0, 0, 0));
	s_matrixVec.push_back(Matrix(Vector4(0, 0, 0, 0), Vector4(0, 0, 0, 0), Vector4(0, 0, 0, 0), Vector4(0, 0, 0, 0)));

	// ident-matrix
	s_matrixVec.push_back(Matrix(Vector4(1, 0, 0, 0), Vector4(0, 1, 0, 0), Vector4(0, 0, 1, 0), Vector4(0, 0, 0, 1)));

	// FLT_MAX
	s_floatVec.push_back(FLT_MAX);
	s_vector2Vec.push_back(Vector2(FLT_MAX, FLT_MAX));
	s_vector3Vec.push_back(Vector3(FLT_MAX, FLT_MAX, FLT_MAX));
	s_vector4Vec.push_back(Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX));
	s_matrixVec.push_back(Matrix(Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX), 
								 Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX), 
								 Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX), 
								 Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX)));

	s_floatVec.push_back(-FLT_MAX);
	s_vector2Vec.push_back(Vector2(-FLT_MAX, -FLT_MAX));
	s_vector3Vec.push_back(Vector3(-FLT_MAX, -FLT_MAX, -FLT_MAX));
	s_vector4Vec.push_back(Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX));
	s_matrixVec.push_back(Matrix(Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX), 
								 Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX), 
								 Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX), 
								 Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)));

	// FLT_MIN
	s_floatVec.push_back(FLT_MIN);
	s_vector2Vec.push_back(Vector2(FLT_MIN, FLT_MIN));
	s_vector3Vec.push_back(Vector3(FLT_MIN, FLT_MIN, FLT_MIN));
	s_vector4Vec.push_back(Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN));
	s_matrixVec.push_back(Matrix(Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN), 
								 Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN), 
								 Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN), 
								 Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN)));

	s_floatVec.push_back(-FLT_MIN);
	s_vector2Vec.push_back(Vector2(-FLT_MIN, -FLT_MIN));
	s_vector3Vec.push_back(Vector3(-FLT_MIN, -FLT_MIN, -FLT_MIN));
	s_vector4Vec.push_back(Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN));
	s_matrixVec.push_back(Matrix(Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN), 
								 Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN), 
								 Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN), 
								 Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN)));

	// inf
	float inf = FLT_MAX * FLT_MAX;
	s_floatVec.push_back(inf);
	s_vector2Vec.push_back(Vector2(inf, inf));
	s_vector3Vec.push_back(Vector3(inf, inf, inf));
	s_vector4Vec.push_back(Vector4(inf, inf, inf, inf));
	s_matrixVec.push_back(Matrix(Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf)));

	// -inf
	inf = -FLT_MAX * FLT_MAX;
	s_floatVec.push_back(inf);
	s_vector2Vec.push_back(Vector2(inf, inf));
	s_vector3Vec.push_back(Vector3(inf, inf, inf));
	s_vector4Vec.push_back(Vector4(inf, inf, inf, inf));
	s_matrixVec.push_back(Matrix(Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf)));

	// Generate the quats from the vec4s
	for(unsigned int i = 0; i < s_vector4Vec.size(); ++i)
	{
		Vector4 v = s_vector4Vec[i];
		v = UnitVec(v);
		s_quatVec.push_back(Quat(v.X(), v.Y(), v.Z(), v.W()));
	}

	// TODO: nans, q-nans & denormalized...
}

// 
// Vector2 tests
//

template <class UnaryOp>
class Vector2UnaryOpTest : public TestCase
{
public:
	Vector2UnaryOpTest() : TestCase(UnaryOp::GetName()) {}
	~Vector2UnaryOpTest() {}

	bool Test() const
	{
		UnaryOp op;
		for (unsigned int i = 0; i < s_vector2Vec.size(); ++i)
		{
			Vector2 a = s_vector2Vec[i];
			if (!op(a))
			{
				printf("a = (%.5f, %.5f)", a.X(), a.Y());
				return false;
			}
		}
		return true;
	}
};

class Vector2Negation
{
public:
	static const char* GetName() { return "Negation"; }
	bool operator() (const Vector2& a)
	{
		float ax = a.X(), ay = a.Y();
		float rx = -ax, ry = -ay;
		Vector2 r = -a;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y());
	}
};

class Vector2Length
{
public:
	static const char* GetName() { return "Length"; }
	bool operator() (const Vector2& a)
	{
		float ax = a.X(), ay = a.Y();
		float r = sqrt((ax * ax) + (ay * ay));
		float len = Len(a);
		return FloatTest(r, len);
	}
};

class Vector2UnitVec
{
public:
	static const char* GetName() { return "UnitVec"; }
	bool operator() (const Vector2& a)
	{
		float ax = a.X(), ay = a.Y();
		float len = sqrt((ax * ax) + (ay * ay));
		float rx = ax / len, ry = ay / len;
		Vector2 r = UnitVec(a);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y());
	}
};

class Vector2ScalarMultiplication
{
public:
	static const char* GetName() { return "Scalar Multiplication"; }
	bool operator() (const Vector2& a)
	{
		float ax = a.X(), ay = a.Y();
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector * scalar) and (scalar * vector)
			float s = s_floatVec[i];
			float r1x = ax * s, r1y = ay * s;
			float r2x = s * ax, r2y = s * ay;
			Vector2 r1 = a * s;
			Vector2 r2 = s * a;
			if (!(FloatTest(r1x, r1.X()) && FloatTest(r1y, r1.Y()) && FloatTest(r2x, r2.X()) && FloatTest(r2y, r2.Y()) && FloatTest(r1.X(), r2.X()) && FloatTest(r1.Y(), r2.Y())))
			{
				printf("s = %.5f", s);
				return false;
			}
		}
		return true;
	}
};

class Vector2ScalarDivision
{
public:
	static const char* GetName() { return "Scalar Division"; }
	bool operator() (const Vector2& a)
	{
		float ax = a.X(), ay = a.Y();
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector / scalar) and (scalar / vector)
			float s = s_floatVec[i];
			float r1x = ax / s, r1y = ay / s;
			float r2x = s / ax, r2y = s / ay;
			Vector2 r1 = a / s;
			Vector2 r2 = s / a;

			if (!(FloatTest(r1x, r1.X()) && FloatTest(r1y, r1.Y()) && FloatTest(r2x, r2.X()) && FloatTest(r2y, r2.Y())))
			{
				printf("s = %.5f\n", s);
				return false;
			}
		}
		return true;
	}
};

template <class BinaryOp>
class Vector2BinaryOpTest : public TestCase
{
public:
	Vector2BinaryOpTest() : TestCase(BinaryOp::GetName()) {}
	~Vector2BinaryOpTest() {}

	bool Test() const
	{
		BinaryOp op;
		for (unsigned int i = 0; i < s_vector2Vec.size(); ++i)
		{
			for (unsigned int j = 0; j < s_vector2Vec.size(); ++j)
			{
				Vector2 a = s_vector2Vec[i];
				Vector2 b = s_vector2Vec[j];
				if (!op(a, b))
				{
					printf("a = (%.5f, %.5f), b = (%.5f, %.5f)", a.X(), a.Y(), b.X(), b.Y());
					return false;
				}
			}
		}
		return true;
	}
};

class Vector2Addition
{
public:
	static const char* GetName() { return "Addition"; }
	bool operator()(const Vector2& a, const Vector2& b) const
	{
		float ax = a.X(), ay = a.Y();
		float bx = b.X(), by = b.Y();
		float rx = ax + bx, ry = ay + by;
		Vector2 r = a + b;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y());
	}
};

class Vector2Subtraction
{
public:
	static const char* GetName() { return "Subtraction"; }
	bool operator()(const Vector2& a, const Vector2& b) const
	{
		float ax = a.X(), ay = a.Y();
		float bx = b.X(), by = b.Y();
		float rx = ax - bx, ry = ay - by;
		Vector2 r = a - b;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y());
	}
};

class Vector2DotProduct
{
public:
	static const char* GetName() { return "Dot Product"; }
	bool operator()(const Vector2& a, const Vector2& b) const
	{
		float ax = a.X(), ay = a.Y();
		float bx = b.X(), by = b.Y();
		float r = (ax * bx) + (ay * by);
		float dot = a * b;
		return FloatTest(dot, r);
	}
};

class Vector2CompMul
{
public:
	static const char* GetName() { return "Comp Mul"; }
	bool operator()(const Vector2& a, const Vector2& b) const
	{
		float ax = a.X(), ay = a.Y();
		float bx = b.X(), by = b.Y();
		float rx = (ax * bx), ry = (ay * by);
		Vector2 r = CompMul(a, b);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y());
	}
};

class Vector2CompDiv
{
public:
	static const char* GetName() { return "Comp Div"; }
	bool operator()(const Vector2& a, const Vector2& b) const
	{
		float ax = a.X(), ay = a.Y();
		float bx = b.X(), by = b.Y();
		float rx = (ax / bx), ry = (ay / by);
		Vector2 r = CompDiv(a, b);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y());
	}
};

//
// Vector3 tests
//

template <class UnaryOp>
class Vector3UnaryOpTest : public TestCase
{
public:
	Vector3UnaryOpTest() : TestCase(UnaryOp::GetName()) {}
	~Vector3UnaryOpTest() {}

	bool Test() const
	{
		UnaryOp op;
		for (unsigned int i = 0; i < s_vector3Vec.size(); ++i)
		{
			Vector3 a = s_vector3Vec[i];
			if (!op(a))
			{
				printf("a = (%.5f, %.5f, %.5f)", a.X(), a.Y(), a.Z());
				return false;
			}
		}
		return true;
	}
};

class Vector3Negation
{
public:
	static const char* GetName() { return "Negation"; }
	bool operator() (const Vector3& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		float rx = -ax, ry = -ay, rz = -az;
		Vector3 r = -a;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z());
	}
};

class Vector3Length
{
public:
	static const char* GetName() { return "Length"; }
	bool operator() (const Vector3& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		float r = sqrt((ax * ax) + (ay * ay) + (az * az));
		float len = Len(a);
		return FloatTest(r, len);
	}
};

class Vector3UnitVec
{
public:
	static const char* GetName() { return "UnitVec"; }
	bool operator() (const Vector3& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		float len = sqrt((ax * ax) + (ay * ay) + (az * az));
		float rx = ax / len, ry = ay / len, rz = az / len;
		Vector3 r = UnitVec(a);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz ,r.Z());
	}
};

class Vector3ScalarMultiplication
{
public:
	static const char* GetName() { return "Scalar Multiplication"; }
	bool operator() (const Vector3& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector * scalar) and (scalar * vector)
			float s = s_floatVec[i];
			float r1x = ax * s, r1y = ay * s, r1z = az * s;
			float r2x = s * ax, r2y = s * ay, r2z = s * az;
			Vector3 r1 = a * s;
			Vector3 r2 = s * a;
			if (!(FloatTest(r1x, r1.X()) && FloatTest(r1y, r1.Y()) && FloatTest(r1z, r1.Z()) && 
				  FloatTest(r2x, r2.X()) && FloatTest(r2y, r2.Y()) && FloatTest(r2z, r2.Z()) &&
				  FloatTest(r1.X(), r2.X()) && FloatTest(r1.Y(), r2.Y()) && FloatTest(r1.Z(), r2.Z())))
			{
				printf("s = %.5f", s);
				return false;
			}
		}
		return true;
	}
};

class Vector3ScalarDivision
{
public:
	static const char* GetName() { return "Scalar Division"; }
	bool operator() (const Vector3& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector / scalar) and (scalar / vector)
			float s = s_floatVec[i];
			float r1x = ax / s, r1y = ay / s, r1z = az / s;
			float r2x = s / ax, r2y = s / ay, r2z = s / az;
			Vector3 r1 = a / s;
			Vector3 r2 = s / a;

			if (!(FloatTest(r1x, r1.X()) && FloatTest(r1y, r1.Y()) && FloatTest(r1z, r1.Z()) && 
				  FloatTest(r2x, r2.X()) && FloatTest(r2y, r2.Y()) && FloatTest(r2z, r2.Z())))
			{
				printf("s = %.5f\n", s);
				return false;
			}
		}
		return true;
	}
};

template <class BinaryOp>
class Vector3BinaryOpTest : public TestCase
{
public:
	Vector3BinaryOpTest() : TestCase(BinaryOp::GetName()) {}
	~Vector3BinaryOpTest() {}

	bool Test() const
	{
		BinaryOp op;
		for (unsigned int i = 0; i < s_vector3Vec.size(); ++i)
		{
			for (unsigned int j = 0; j < s_vector3Vec.size(); ++j)
			{
				Vector3 a = s_vector3Vec[i];
				Vector3 b = s_vector3Vec[j];
				if (!op(a, b))
				{
					printf("a = (%.5f, %.5f, %.5f), b = (%.5f, %.5f, %.5f)", a.X(), a.Y(), a.Z(), b.X(), b.Y(), b.Z());
					return false;
				}
			}
		}
		return true;
	}
};

class Vector3Addition
{
public:
	static const char* GetName() { return "Addition"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		float bx = b.X(), by = b.Y(), bz = b.Z();
		float rx = ax + bx, ry = ay + by, rz = az + bz;
		Vector3 r = a + b;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z());
	}
};

class Vector3Subtraction
{
public:
	static const char* GetName() { return "Subtraction"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		float bx = b.X(), by = b.Y(), bz = b.Z();
		float rx = ax - bx, ry = ay - by, rz = az - bz;
		Vector3 r = a - b;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z());
	}
};

class Vector3DotProduct
{
public:
	static const char* GetName() { return "Dot Product"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		float bx = b.X(), by = b.Y(), bz = b.Z();
		float r = (ax * bx) + (ay * by) + (az * bz);
		float dot = a * b;
		return FloatTest(dot, r);
	}
};

class Vector3CompMul
{
public:
	static const char* GetName() { return "Comp Mul"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		float bx = b.X(), by = b.Y(), bz = b.Z();
		float rx = (ax * bx), ry = (ay * by), rz = (az * bz);
		Vector3 r = CompMul(a, b);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z());
	}
};

class Vector3CompDiv
{
public:
	static const char* GetName() { return "Comp Div"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		float bx = b.X(), by = b.Y(), bz = b.Z();
		float rx = (ax / bx), ry = (ay / by), rz = (az / bz);
		Vector3 r = CompDiv(a, b);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z());
	}
};

class Vector3CrossProduct
{
public:
	static const char* GetName() { return "Cross Product"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z();
		float bx = b.X(), by = b.Y(), bz = b.Z();

		float rx = (ay * bz) - (az * by);
		float ry = (az * bx) - (ax * bz);
		float rz = (ax * by) - (ay * bx);

		Vector3 r = a % b;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z());
	}
};

//
// Vector4 tests
//

template <class UnaryOp>
class Vector4UnaryOpTest : public TestCase
{
public:
	Vector4UnaryOpTest() : TestCase(UnaryOp::GetName()) {}
	~Vector4UnaryOpTest() {}

	bool Test() const
	{
		UnaryOp op;
		for (unsigned int i = 0; i < s_vector4Vec.size(); ++i)
		{
			Vector4 a = s_vector4Vec[i];
			if (!op(a))
			{
				printf("a = (%.5f, %.5f, %.5f, %.5f)", a.X(), a.Y(), a.Z(), a.W());
				return false;
			}
		}
		return true;
	}
};

class Vector4Negation
{
public:
	static const char* GetName() { return "Negation"; }
	bool operator() (const Vector4& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float rx = -ax, ry = -ay, rz = -az, rw = -aw;
		Vector4 r = -a;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class Vector4Length
{
public:
	static const char* GetName() { return "Length"; }
	bool operator() (const Vector4& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float r = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
		float len = Len(a);
		return FloatTest(r, len);
	}
};

class Vector4UnitVec
{
public:
	static const char* GetName() { return "UnitVec"; }
	bool operator() (const Vector4& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float len = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
		float rx = ax / len, ry = ay / len, rz = az / len, rw = aw / len;
		Vector4 r = UnitVec(a);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class Vector4ScalarMultiplication
{
public:
	static const char* GetName() { return "Scalar Multiplication"; }
	bool operator() (const Vector4& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector * scalar) and (scalar * vector)
			float s = s_floatVec[i];
			float r1x = ax * s, r1y = ay * s, r1z = az * s, r1w = aw * s;
			float r2x = s * ax, r2y = s * ay, r2z = s * az, r2w = s * aw;
			Vector4 r1 = a * s;
			Vector4 r2 = s * a;
			if (!(FloatTest(r1x, r1.X()) && FloatTest(r1y, r1.Y()) && FloatTest(r1z, r1.Z()) && FloatTest(r1w, r1.W()) &&
				  FloatTest(r2x, r2.X()) && FloatTest(r2y, r2.Y()) && FloatTest(r2z, r2.Z()) && FloatTest(r2w, r2.W()) &&
				  FloatTest(r1.X(), r2.X()) && FloatTest(r1.Y(), r2.Y()) && FloatTest(r1.Z(), r2.Z()) && FloatTest(r1.W(), r2.W())))
			{
				printf("s = %.5f", s);
				return false;
			}
		}
		return true;
	}
};

class Vector4ScalarDivision
{
public:
	static const char* GetName() { return "Scalar Division"; }
	bool operator() (const Vector4& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector / scalar) and (scalar / vector)
			float s = s_floatVec[i];
			float r1x = ax / s, r1y = ay / s, r1z = az / s, r1w = aw / s;
			float r2x = s / ax, r2y = s / ay, r2z = s / az, r2w = s / aw;
			Vector4 r1 = a / s;
			Vector4 r2 = s / a;

			if (!(FloatTest(r1x, r1.X()) && FloatTest(r1y, r1.Y()) && FloatTest(r1z, r1.Z()) && FloatTest(r1w, r1.W()) &&
				  FloatTest(r2x, r2.X()) && FloatTest(r2y, r2.Y()) && FloatTest(r2z, r2.Z()) && FloatTest(r2w, r2.W())))
			{
				printf("s = %.5f\n", s);
				return false;
			}
		}
		return true;
	}
};

template <class BinaryOp>
class Vector4BinaryOpTest : public TestCase
{
public:
	Vector4BinaryOpTest() : TestCase(BinaryOp::GetName()) {}
	~Vector4BinaryOpTest() {}

	bool Test() const
	{
		BinaryOp op;
		for (unsigned int i = 0; i < s_vector4Vec.size(); ++i)
		{
			for (unsigned int j = 0; j < s_vector4Vec.size(); ++j)
			{
				Vector4 a = s_vector4Vec[i];
				Vector4 b = s_vector4Vec[j];
				if (!op(a, b))
				{
					printf("a = (%.5f, %.5f, %.5f, %.5f), b = (%.5f, %.5f, %.5f, %.5f)", a.X(), a.Y(), a.Z(), a.W(), b.X(), b.Y(), b.Z(), b.W());
					return false;
				}
			}
		}
		return true;
	}
};

class Vector4Addition
{
public:
	static const char* GetName() { return "Addition"; }
	bool operator()(const Vector4& a, const Vector4& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float bx = b.X(), by = b.Y(), bz = b.Z(), bw = b.W();
		float rx = ax + bx, ry = ay + by, rz = az + bz, rw = aw + bw;
		Vector4 r = a + b;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class Vector4Subtraction
{
public:
	static const char* GetName() { return "Subtraction"; }
	bool operator()(const Vector4& a, const Vector4& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float bx = b.X(), by = b.Y(), bz = b.Z(), bw = b.W();
		float rx = ax - bx, ry = ay - by, rz = az - bz, rw = aw - bw;
		Vector4 r = a - b;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class Vector4DotProduct
{
public:
	static const char* GetName() { return "Dot Product"; }
	bool operator()(const Vector4& a, const Vector4& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float bx = b.X(), by = b.Y(), bz = b.Z(), bw = b.W();
		float r = (ax * bx) + (ay * by) + (az * bz) + (aw * bw);
		float dot = a * b;
		return FloatTest(dot, r);
	}
};

class Vector4CompMul
{
public:
	static const char* GetName() { return "Comp Mul"; }
	bool operator()(const Vector4& a, const Vector4& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float bx = b.X(), by = b.Y(), bz = b.Z(), bw = b.W();
		float rx = (ax * bx), ry = (ay * by), rz = (az * bz), rw = (aw * bw);
		Vector4 r = CompMul(a, b);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class Vector4CompDiv
{
public:
	static const char* GetName() { return "Comp Div"; }
	bool operator()(const Vector4& a, const Vector4& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float bx = b.X(), by = b.Y(), bz = b.Z(), bw = b.W();
		float rx = (ax / bx), ry = (ay / by), rz = (az / bz), rw = (aw / bw);
		Vector4 r = CompDiv(a, b);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

//
// Quat tests
//

template <class UnaryOp>
class QuatUnaryOpTest : public TestCase
{
public:
	QuatUnaryOpTest() : TestCase(UnaryOp::GetName()) {}
	~QuatUnaryOpTest() {}

	bool Test() const
	{
		UnaryOp op;
		for (unsigned int i = 0; i < s_quatVec.size(); ++i)
		{
			Quat a = s_quatVec[i];
			if (!op(a))
			{
				printf("a = (%.5f, %.5f, %.5f, %.5f)", a.X(), a.Y(), a.Z(), a.W());
				return false;
			}
		}
		return true;
	}
};

void QuatMul(float* rx, float* ry, float* rz, float* rw, float ax, float ay, float az, float aw, float bx, float by, float bz, float bw)
{
    *rx =  ax * bw + ay * bz - az * by + aw * bx;
    *ry = -ax * bz + ay * bw + az * bx + aw * by;
    *rz =  ax * by - ay * bx + az * bw + aw * bz;
    *rw = -ax * bx - ay * by - az * bz + aw * bw;
}

class QuatRotate
{
public:
	static const char* GetName() { return "Rotate"; }
	bool operator() (const Quat& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float rx, ry, rz, rw;

		for (unsigned int i = 0; i < s_vector3Vec.size(); ++i)
		{
			Vector3 v = s_vector3Vec[i];
			float bx = v.X(), by = v.Y(), bz = v.Z(), bw = 0.0f;

			QuatMul(&rx, &ry, &rz, &rw, ax, ay, az, aw, bx, by, bz, bw);
			QuatMul(&rx, &ry, &rz, &rw, rx, ry, rz, rw, -ax, -ay, -az, aw);

			Vector3 r = Rotate(a, v);
			if (!(FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z())))
			{
				printf("v = (%.5f, %.5f, %.5f)", v.X(), v.Y(), v.Z());
				return false;
			}
		}
		return true;
	}
};

class QuatConjugate
{
public:
	static const char* GetName() { return "Conjugate"; }
	bool operator() (const Quat& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float rx = -ax, ry = -ay, rz = -az, rw = aw;
		Quat r = ~a;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class QuatNegation
{
public:
	static const char* GetName() { return "Negation"; }
	bool operator() (const Quat& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float rx = -ax, ry = -ay, rz = -az, rw = -aw;
		Quat r = -a;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class QuatScalarMultiplication
{
public:
	static const char* GetName() { return "Scalar Multiplication"; }
	bool operator() (const Quat& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector * scalar) and (scalar * vector)
			float s = s_floatVec[i];
			float r1x = ax * s, r1y = ay * s, r1z = az * s, r1w = aw * s;
			float r2x = s * ax, r2y = s * ay, r2z = s * az, r2w = s * aw;
			Quat r1 = a * s;
			Quat r2 = s * a;
			if (!(FloatTest(r1x, r1.X()) && FloatTest(r1y, r1.Y()) && FloatTest(r1z, r1.Z()) && FloatTest(r1w, r1.W()) &&
				  FloatTest(r2x, r2.X()) && FloatTest(r2y, r2.Y()) && FloatTest(r2z, r2.Z()) && FloatTest(r2w, r2.W()) &&
				  FloatTest(r1.X(), r2.X()) && FloatTest(r1.Y(), r2.Y()) && FloatTest(r1.Z(), r2.Z()) && FloatTest(r1.W(), r2.W())))
			{
				printf("s = %.5f", s);
				return false;
			}
		}
		return true;
	}
};

class QuatScalarDivision
{
public:
	static const char* GetName() { return "Scalar Division"; }
	bool operator() (const Quat& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector / scalar) and (scalar / vector)
			float s = s_floatVec[i];
			float r1x = ax / s, r1y = ay / s, r1z = az / s, r1w = aw / s;
			float r2x = s / ax, r2y = s / ay, r2z = s / az, r2w = s / aw;
			Quat r1 = a / s;
			Quat r2 = s / a;

			if (!(FloatTest(r1x, r1.X()) && FloatTest(r1y, r1.Y()) && FloatTest(r1z, r1.Z()) && FloatTest(r1w, r1.W()) &&
				  FloatTest(r2x, r2.X()) && FloatTest(r2y, r2.Y()) && FloatTest(r2z, r2.Z()) && FloatTest(r2w, r2.W())))
			{
				printf("s = %.5f\n", s);
				return false;
			}
		}
		return true;
	}
};

class QuatLength
{
public:
	static const char* GetName() { return "Length"; }
	bool operator() (const Quat& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float r = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
		float len = Len(a);
		return FloatTest(r, len);
	}
};

class QuatUnitVec
{
public:
	static const char* GetName() { return "UnitVec"; }
	bool operator() (const Quat& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float len = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
		float rx = ax / len, ry = ay / len, rz = az / len, rw = aw / len;
		Quat r = UnitVec(a);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

void QuatExp(float* rx, float* ry, float* rz, float* rw, float qx, float qy, float qz, float qw)
{
	float angle = sqrt((qx * qx) + (qy * qy) + (qz * qz));
	float sin_a = sin(angle / 2.0f);

	*rx = (qx / angle) * sin_a;
	*ry = (qy / angle) * sin_a;
	*rz = (qz / angle) * sin_a;
	*rw = cos(angle / 2.0f);
}

class QuatExponential
{
public:
	static const char* GetName() { return "Quat Exponential"; }
	bool operator() (const Quat& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float rx, ry, rz, rw;
		QuatExp(&rx, &ry, &rz, &rw, ax, ay, az, aw);
		Quat r = QuatExp(a);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

void QuatLog(float* rx, float* ry, float* rz, float* rw, float qx, float qy, float qz, float qw)
{
	float cos_a = qw;
	if (cos_a > 1.0f) cos_a = 1.0f;
	if (cos_a < -1.0f) cos_a = -1.0f;

    float sin_a = (float)sqrt(1.0f - cos_a * cos_a);

    if (fabs(sin_a) < 0.0005f)
		sin_a = 1.0f;
	else
		sin_a = 1.f/sin_a;

    float angle = 2.0f * (float)acos(cos_a);

    *rx = qx * sin_a * angle;
    *ry = qy * sin_a * angle;
    *rz = qz * sin_a * angle;
	*rw = 0.0f;
}

class QuatLogarithm
{
public:
	static const char* GetName() { return "Quat Logarithm"; }
	bool operator() (const Quat& a)
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float rx, ry, rz, rw;
		QuatLog(&rx, &ry, &rz, &rw, ax, ay, az, aw);
		Quat r = QuatLog(a);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

template <class BinaryOp>
class QuatBinaryOpTest : public TestCase
{
public:
	QuatBinaryOpTest() : TestCase(BinaryOp::GetName()) {}
	~QuatBinaryOpTest() {}

	bool Test() const
	{
		BinaryOp op;
		for (unsigned int i = 0; i < s_quatVec.size(); ++i)
		{
			for (unsigned int j = 0; j < s_quatVec.size(); ++j)
			{
				Quat a = s_quatVec[i];
				Quat b = s_quatVec[j];
				if (!op(a, b))
				{
					printf("a = (%.5f, %.5f, %.5f, %.5f), b = (%.5f, %.5f, %.5f, %.5f)", a.X(), a.Y(), a.Z(), a.W(), b.X(), b.Y(), b.Z(), b.W());
					return false;
				}
			}
		}
		return true;
	}
};

class QuatAddition
{
public:
	static const char* GetName() { return "Addition"; }
	bool operator()(const Quat& a, const Quat& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float bx = b.X(), by = b.Y(), bz = b.Z(), bw = b.W();
		float rx = ax + bx, ry = ay + by, rz = az + bz, rw = aw + bw;
		Quat r = a + b;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class QuatSubtraction
{
public:
	static const char* GetName() { return "Subtraction"; }
	bool operator()(const Quat& a, const Quat& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float bx = b.X(), by = b.Y(), bz = b.Z(), bw = b.W();
		float rx = ax - bx, ry = ay - by, rz = az - bz, rw = aw - bw;
		Quat r = a - b;
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class QuatCompMul
{
public:
	static const char* GetName() { return "Comp Mul"; }
	bool operator()(const Quat& a, const Quat& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float bx = b.X(), by = b.Y(), bz = b.Z(), bw = b.W();
		float rx = (ax * bx), ry = (ay * by), rz = (az * bz), rw = (aw * bw);
		Quat r = CompMul(a, b);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

class QuatCompDiv
{
public:
	static const char* GetName() { return "Comp Div"; }
	bool operator()(const Quat& a, const Quat& b) const
	{
		float ax = a.X(), ay = a.Y(), az = a.Z(), aw = a.W();
		float bx = b.X(), by = b.Y(), bz = b.Z(), bw = b.W();
		float rx = (ax / bx), ry = (ay / by), rz = (az / bz), rw = (aw / bw);
		Quat r = CompDiv(a, b);
		return FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W());
	}
};

//
// Matrix tests
//

template <class UnaryOp>
class MatrixUnaryOpTest : public TestCase
{
public:
	MatrixUnaryOpTest() : TestCase(UnaryOp::GetName()) {}
	~MatrixUnaryOpTest() {}

	bool Test() const
	{
		UnaryOp op;
		for (unsigned int i = 0; i < s_matrixVec.size(); ++i)
		{
			Matrix m = s_matrixVec[i];
			if (!op(m))
			{
				Vector4 row0 = m.row0;
				Vector4 row1 = m.row1;
				Vector4 row2 = m.row2;
				Vector4 row3 = m.row3;
				printf("m = (%.5f, %.5f, %.5f, %.5f)", row0.X(), row0.Y(), row0.Z(), row0.W());
				printf(", (%.5f, %.5f, %.5f, %.5f)", row1.X(), row1.Y(), row1.Z(), row1.W());
				printf(", (%.5f, %.5f, %.5f, %.5f)", row2.X(), row2.Y(), row2.Z(), row2.W());
				printf(", (%.5f, %.5f, %.5f, %.5f)", row3.X(), row3.Y(), row3.Z(), row3.W());
				return false;
			}
		}
		return true;
	}
};

void MatrixToFloatVec(float *fv, const Matrix& m)
{
	for(int r = 0; r < 4; ++r)
	{
		for(int c = 0; c < 4; ++c)
		{
			fv[4 * r + c] = m.Elem(r,c);
		}
	}
}

class MatrixTransform3x3
{
public:
	static const char* GetName() { return "Transform3x3"; }
	bool operator() (const Matrix& a)
	{
		float fv[16];
		MatrixToFloatVec(fv, a);

		for (unsigned int i = 0; i < s_vector3Vec.size(); ++i)
		{
			Vector3 v = s_vector3Vec[i];
			float bx = v.X(), by = v.Y(), bz = v.Z();

			float rx = fv[0] * bx + fv[1] * by + fv[2] * bz;
			float ry = fv[4] * bx + fv[5] * by + fv[6] * bz;
			float rz = fv[8] * bx + fv[9] * by + fv[10] * bz;

			Vector3 r = Transform3x3(a, v);
			if (!(FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z())))
			{
				printf("v = (%.5f, %.5f, %.5f)", v.X(), v.Y(), v.Z());
				return false;
			}
		}
		return true;
	}
};

class MatrixTransform3x4
{
public:
	static const char* GetName() { return "Transform3x4"; }
	bool operator() (const Matrix& a)
	{
		float fv[16];
		MatrixToFloatVec(fv, a);

		for (unsigned int i = 0; i < s_vector3Vec.size(); ++i)
		{
			Vector3 v = s_vector3Vec[i];
			float bx = v.X(), by = v.Y(), bz = v.Z();

			float rx = fv[0] * bx + fv[1] * by + fv[2] * bz + fv[3];
			float ry = fv[4] * bx + fv[5] * by + fv[6] * bz + fv[7];
			float rz = fv[8] * bx + fv[9] * by + fv[10] * bz + fv[11];

			Vector3 r = Transform3x4(a, v);
			if (!(FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z())))
			{
				printf("v = (%.5f, %.5f, %.5f)", v.X(), v.Y(), v.Z());
				return false;
			}
		}
		return true;
	}
};

class MatrixTransform4x4
{
public:
	static const char* GetName() { return "Transform4x4"; }
	bool operator() (const Matrix& a)
	{
		float fv[16];
		MatrixToFloatVec(fv, a);

		for (unsigned int i = 0; i < s_vector4Vec.size(); ++i)
		{
			Vector4 v = s_vector4Vec[i];
			float bx = v.X(), by = v.Y(), bz = v.Z(), bw = v.W();

			float rx = fv[0] * bx + fv[1] * by + fv[2] * bz + fv[3] * bw;
			float ry = fv[4] * bx + fv[5] * by + fv[6] * bz + fv[7] * bw;
			float rz = fv[8] * bx + fv[9] * by + fv[10] * bz + fv[11] * bw;
			float rw = fv[12] * bx + fv[13] * by + fv[14] * bz + fv[15] * bw;

			Vector4 r = Transform4x4(a, v);
			if (!(FloatTest(rx, r.X()) && FloatTest(ry, r.Y()) && FloatTest(rz, r.Z()) && FloatTest(rw, r.W())))
			{
				printf("v = (%.5f, %.5f, %.5f, %.5f)", v.X(), v.Y(), v.Z(), v.W());
				return false;
			}
		}
		return true;
	}
};

template <class BinaryOp>
class MatrixBinaryOpTest : public TestCase
{
public:
	MatrixBinaryOpTest() : TestCase(BinaryOp::GetName()) {}
	~MatrixBinaryOpTest() {}

	bool Test() const
	{
		BinaryOp op;
		for (unsigned int i = 0; i < s_matrixVec.size(); ++i)
		{
			for (unsigned int j = 0; j < s_matrixVec.size(); ++j)
			{
				Matrix a = s_matrixVec[i];
				Matrix b = s_matrixVec[j];
				if (!op(a, b))
				{
					Vector4 row0 = a.row0;
					Vector4 row1 = a.row1;
					Vector4 row2 = a.row2;
					Vector4 row3 = a.row3;
					printf("a = (%.5f, %.5f, %.5f, %.5f)", row0.X(), row0.Y(), row0.Z(), row0.W());
					printf(", (%.5f, %.5f, %.5f, %.5f)", row1.X(), row1.Y(), row1.Z(), row1.W());
					printf(", (%.5f, %.5f, %.5f, %.5f)", row2.X(), row2.Y(), row2.Z(), row2.W());
					printf(", (%.5f, %.5f, %.5f, %.5f)", row3.X(), row3.Y(), row3.Z(), row3.W());

					row0 = b.row0;
					row1 = b.row1;
					row2 = b.row2;
					row3 = b.row3;
					printf(" b = (%.5f, %.5f, %.5f, %.5f)", row0.X(), row0.Y(), row0.Z(), row0.W());
					printf(", (%.5f, %.5f, %.5f, %.5f)", row1.X(), row1.Y(), row1.Z(), row1.W());
					printf(", (%.5f, %.5f, %.5f, %.5f)", row2.X(), row2.Y(), row2.Z(), row2.W());
					printf(", (%.5f, %.5f, %.5f, %.5f)", row3.X(), row3.Y(), row3.Z(), row3.W());
					return false;
				}
			}
		}
		return true;
	}
};



class MatrixMultiplication
{
public:
	static const char* GetName() { return "Multiplication"; }
	bool operator()(const Matrix& a, const Matrix& b) const
	{
		float afv[16], bfv[16];
		MatrixToFloatVec(afv, a);
		MatrixToFloatVec(bfv, b);

		float rfv[16];
		static int ri[4][4] = {{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}};
		static int ci[4][4] = {{0, 4, 8, 12}, {1, 5, 9, 13}, {2, 6, 10, 14}, {3, 7, 11, 15}};
		for (int r = 0; r < 4; ++r)
		{
			for (int c = 0; c < 4; ++c)
			{
				float terms[4];				
				for (int i = 0; i < 4; ++i)
					terms[i] = afv[ri[r][i]] * bfv[ci[c][i]];
				rfv[r * 4 + c] = terms[0] + terms[1] + terms[2] + terms[3];
			}
		}

		Matrix r = a * b;

		return FloatTest(rfv[0], r.Elem(0,0)) && FloatTest(rfv[1], r.Elem(0,1)) && FloatTest(rfv[2], r.Elem(0,2)) && FloatTest(rfv[3], r.Elem(0,3)) &&
		    FloatTest(rfv[4], r.Elem(1,0)) && FloatTest(rfv[5], r.Elem(1,1)) && FloatTest(rfv[6], r.Elem(1,2)) && FloatTest(rfv[7], r.Elem(1,3)) &&
			FloatTest(rfv[8], r.Elem(2,0)) && FloatTest(rfv[9], r.Elem(2,1)) && FloatTest(rfv[10], r.Elem(2,2)) && FloatTest(rfv[11], r.Elem(2,3)) &&
			FloatTest(rfv[12], r.Elem(3,0)) && FloatTest(rfv[13], r.Elem(3,1)) && FloatTest(rfv[14], r.Elem(3,2)) && FloatTest(rfv[15], r.Elem(3,3));
	}
};

int main(int argc, char* argv[])
{
	InitTestData();

	TestSuite vector2Suite("Vector2");
	vector2Suite.AddTest(new Vector2UnaryOpTest<Vector2Negation>());
	vector2Suite.AddTest(new Vector2UnaryOpTest<Vector2Length>());
	vector2Suite.AddTest(new Vector2UnaryOpTest<Vector2UnitVec>());
	vector2Suite.AddTest(new Vector2UnaryOpTest<Vector2ScalarMultiplication>());
	vector2Suite.AddTest(new Vector2UnaryOpTest<Vector2ScalarDivision>());
	vector2Suite.AddTest(new Vector2BinaryOpTest<Vector2Addition>());
	vector2Suite.AddTest(new Vector2BinaryOpTest<Vector2Subtraction>());
	vector2Suite.AddTest(new Vector2BinaryOpTest<Vector2DotProduct>());
	vector2Suite.AddTest(new Vector2BinaryOpTest<Vector2CompMul>());
	vector2Suite.AddTest(new Vector2BinaryOpTest<Vector2CompDiv>());
	vector2Suite.RunTests();

	TestSuite vector3Suite("Vector3");
	vector3Suite.AddTest(new Vector3UnaryOpTest<Vector3Negation>());
	vector3Suite.AddTest(new Vector3UnaryOpTest<Vector3Length>());
	vector3Suite.AddTest(new Vector3UnaryOpTest<Vector3UnitVec>());
	vector3Suite.AddTest(new Vector3UnaryOpTest<Vector3ScalarMultiplication>());
	vector3Suite.AddTest(new Vector3UnaryOpTest<Vector3ScalarDivision>());
	vector3Suite.AddTest(new Vector3BinaryOpTest<Vector3Addition>());
	vector3Suite.AddTest(new Vector3BinaryOpTest<Vector3Subtraction>());
	vector3Suite.AddTest(new Vector3BinaryOpTest<Vector3DotProduct>());
	vector3Suite.AddTest(new Vector3BinaryOpTest<Vector3CompMul>());
	vector3Suite.AddTest(new Vector3BinaryOpTest<Vector3CompDiv>());
	vector3Suite.AddTest(new Vector3BinaryOpTest<Vector3CrossProduct>());
	vector3Suite.RunTests();

	TestSuite vector4Suite("Vector4");
	vector4Suite.AddTest(new Vector4UnaryOpTest<Vector4Negation>());
	vector4Suite.AddTest(new Vector4UnaryOpTest<Vector4Length>());
	vector4Suite.AddTest(new Vector4UnaryOpTest<Vector4UnitVec>());
	vector4Suite.AddTest(new Vector4UnaryOpTest<Vector4ScalarMultiplication>());
	vector4Suite.AddTest(new Vector4UnaryOpTest<Vector4ScalarDivision>());
	vector4Suite.AddTest(new Vector4BinaryOpTest<Vector4Addition>());
	vector4Suite.AddTest(new Vector4BinaryOpTest<Vector4Subtraction>());
	vector4Suite.AddTest(new Vector4BinaryOpTest<Vector4DotProduct>());
	vector4Suite.AddTest(new Vector4BinaryOpTest<Vector4CompMul>());
	vector4Suite.AddTest(new Vector4BinaryOpTest<Vector4CompDiv>());
	vector4Suite.RunTests();

	TestSuite quatSuite("Quat");
	quatSuite.AddTest(new QuatUnaryOpTest<QuatRotate>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatConjugate>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatNegation>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatScalarMultiplication>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatScalarDivision>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatLength>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatUnitVec>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatExponential>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatLogarithm>());
	quatSuite.AddTest(new QuatBinaryOpTest<QuatAddition>());
	quatSuite.AddTest(new QuatBinaryOpTest<QuatSubtraction>());
	quatSuite.AddTest(new QuatBinaryOpTest<QuatCompMul>());
	quatSuite.AddTest(new QuatBinaryOpTest<QuatCompDiv>());
	quatSuite.RunTests();

	TestSuite matrixSuite("Matrix");
	matrixSuite.AddTest(new MatrixUnaryOpTest<MatrixTransform3x3>());
	matrixSuite.AddTest(new MatrixUnaryOpTest<MatrixTransform3x4>());
	matrixSuite.AddTest(new MatrixUnaryOpTest<MatrixTransform4x4>());
	matrixSuite.AddTest(new MatrixBinaryOpTest<MatrixMultiplication>());
	
	// TODO: matrix from quat
	// TODO: matrix from quat & trans
	// TODO: matrix from quat, trans & scale
	// TODO: matrix from axis angle.
	// TODO: matrix set scale.
	// TODO: make projection
	// TODO: make look-at
	// TODO: make ident
	// TODO: matrix addition
	// TODO: matrix subtraction
	// TODO: Inverse
	// TODO: OrthonormalInverse
//	matrixSuite.RunTests();

	// TODO: Rad2Deg
	// TODO: Deg2Rad

	// quick matrix test
	Matrix m(Quat(Vector3(0.3f,1.5f,0.7f), -PI/3.0f), Vector3(1.0, 2.0f, 3.0f));	
	Matrix m_inv;
	if (!Inverse(m_inv, m))
		printf("Could not find inverse!\n");
	else
		printf("Found inverse!\n");

	// Check if product is identity.
	Matrix i = m * m_inv;

	if (FuzzyEqual(i.row0[0], 1.0f) && FuzzyEqual(i.row0[1], 0.0f) && FuzzyEqual(i.row0[2], 0.0f) && FuzzyEqual(i.row0[3], 0.0f) &&
		FuzzyEqual(i.row1[0], 0.0f) && FuzzyEqual(i.row1[1], 1.0f) && FuzzyEqual(i.row1[2], 0.0f) && FuzzyEqual(i.row1[3], 0.0f) &&
		FuzzyEqual(i.row2[0], 0.0f) && FuzzyEqual(i.row2[1], 0.0f) && FuzzyEqual(i.row2[2], 1.0f) && FuzzyEqual(i.row2[3], 0.0f) &&
		FuzzyEqual(i.row3[0], 0.0f) && FuzzyEqual(i.row3[1], 0.0f) && FuzzyEqual(i.row3[2], 0.0f) && FuzzyEqual(i.row3[3], 1.0f))
		printf("Product test passed!\n");
	else
		printf("Product test FAIL!\n");

	return 0;
}
