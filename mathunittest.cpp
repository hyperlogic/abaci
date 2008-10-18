#include "math3d.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "unittest.h"

#include <CoreServices/CoreServices.h>
#include <mach/mach.h>
#include <mach/mach_time.h>

// static test containers
std::vector<float> s_floatVec;
std::vector<Vector2> s_vector2Vec;
std::vector<Vector3> s_vector3Vec;
std::vector<Vector4> s_vector4Vec;
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
	srand((unsigned int)mach_absolute_time());

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

	return 0;
}
