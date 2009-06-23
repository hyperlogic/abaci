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
#include <math.h>
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

bool FuzzyFloatTest(float rhs, float lhs, float epsilon=0.0001f)
{
	// NaN's are equal, and so are Infs
	if ((isnan(rhs) && isnan(lhs)) ||
		(isinf(rhs) && isinf(lhs)))
		return true;
	else
		return FuzzyEqual(rhs, lhs, epsilon);
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
		v = v.Unit();
		s_quatVec.push_back(Quat(v.x, v.y, v.z, v.w));
	}

	// TODO: nans, q-nans & denormalized...
}

//
// Float tests
//

template <class UnaryOp>
class FloatUnaryOpTest : public TestCase
{
public:
	FloatUnaryOpTest() : TestCase(UnaryOp::GetName()) {}
	~FloatUnaryOpTest() {}

	bool Test() const
	{
		UnaryOp op;
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			float a = s_floatVec[i];
			if (!op(a))
			{
				printf("a = %.5f", a);
				return false;
			}
		}
		return true;
	}
};

class FloatDegToRad
{
public:
	static const char* GetName() { return "DegToRad"; }
	bool operator() (float a)
	{
		float r = a * (PI / 180.0f);
		bool retval = FuzzyFloatTest(DegToRad(a), r);
		if (!retval)
			printf("r = %.5f\n", r);
		return retval;
	}
};

class FloatRadToDeg
{
public:
	static const char* GetName() { return "RadToDeg"; }
	bool operator() (float a)
	{
		float r = a * (180.0f / PI);
		return FuzzyFloatTest(RadToDeg(a), r);
	}
};

class FloatLimitPi
{
public:
	static const char* GetName() { return "LimitPi"; }
	bool operator() (float a)
	{
		float r = a;
		if ((a > PI) || (a < -PI))
			r = fmod(a + PI, 2.0f * PI) - PI;
		return FuzzyFloatTest(LimitPi(a), r);
	}
};

class FloatMod2Pi
{
public:
	static const char* GetName() { return "Mod2Pi"; }
	bool operator() (float a)
	{
		float r = a;
		if ((a > 0) || (a < 2.0f * PI))
			r = fmod(a, 2.0f * PI);
		return FuzzyFloatTest(Mod2Pi(a), r);
	}
};

template <class BinaryOp>
class FloatBinaryOpTest : public TestCase
{
public:
	FloatBinaryOpTest() : TestCase(BinaryOp::GetName()) {}
	~FloatBinaryOpTest() {}

	bool Test() const
	{
		BinaryOp op;
		for (unsigned int i = 0; i < s_vector2Vec.size(); ++i)
		{
			for (unsigned int j = 0; j < s_vector2Vec.size(); ++j)
			{
				float a = s_floatVec[i];
				float b = s_floatVec[j];
				if (!op(a, b))
				{
					printf("a = %.5f b = %.5f", a, b);
					return false;
				}
			}
		}
		return true;
	}
};

class FloatFuzzyEqual
{
public:
	static const char* GetName() { return "FuzzyEqual"; }
	bool operator() (float a, float b)
	{
		const float epsilon = 0.001f;
		bool r = fabs(a - b) <= epsilon;
		return !(r ^ FuzzyEqual(a, b, epsilon));
	}
};


template <class TernaryOp>
class FloatTernaryOpTest : public TestCase
{
public:
	FloatTernaryOpTest() : TestCase(TernaryOp::GetName()) {}
	~FloatTernaryOpTest() {}

	bool Test() const
	{
		TernaryOp op;
		for (unsigned int i = 0; i < s_vector2Vec.size(); ++i)
		{
			for (unsigned int j = 0; j < s_vector2Vec.size(); ++j)
			{
				for (unsigned int k = 0; k < s_vector2Vec.size(); ++k)
				{
					float a = s_floatVec[i];
					float b = s_floatVec[j];
					float c = s_floatVec[k];
					if (!op(a, b, c))
					{
						printf("a = %.5f b = %.5f, c = %.5f", a, b, c);
						return false;
					}
				}
			}
		}
		return true;
	}
};

class FloatClamp
{
public:
	static const char* GetName() { return "Clamp"; }
	bool operator() (float a, float min, float max)
	{
		float r = a;

		// don't check degenerate values. I'm assuming the user does the right thing.
		if (min > max)
			return true;

		if (a < min)
			r = min;
		else if (a > max)
			r = max;

		return FloatTest(Clamp(a, min, max), r);
	}
};

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
				printf("a = (%.5f, %.5f)", a.x, a.y);
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
		float ax = a.x, ay = a.y;
		float rx = -ax, ry = -ay;
		Vector2 r = -a;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y);
	}
};

class Vector2Length
{
public:
	static const char* GetName() { return "Length"; }
	bool operator() (const Vector2& a)
	{
		float ax = a.x, ay = a.y;
		float r = sqrt((ax * ax) + (ay * ay));
		float len = a.Len();
		return FloatTest(r, len);
	}
};

class Vector2UnitVec
{
public:
	static const char* GetName() { return "UnitVec"; }
	bool operator() (const Vector2& a)
	{
		float ax = a.x, ay = a.y;
		float len = sqrt((ax * ax) + (ay * ay));
		float rx = ax / len, ry = ay / len;
		Vector2 r = a.Unit();
		return FloatTest(rx, r.x) && FloatTest(ry, r.y);
	}
};

class Vector2ScalarMultiplication
{
public:
	static const char* GetName() { return "Scalar Multiplication"; }
	bool operator() (const Vector2& a)
	{
		float ax = a.x, ay = a.y;
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector * scalar) and (scalar * vector)
			float s = s_floatVec[i];
			float r1x = ax * s, r1y = ay * s;
			float r2x = s * ax, r2y = s * ay;
			Vector2 r1 = a * s;
			Vector2 r2 = s * a;
			if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r1.x, r2.x) && FloatTest(r1.y, r2.y)))
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
		float ax = a.x, ay = a.y;
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector / scalar) and (scalar / vector)
			float s = s_floatVec[i];
			float r1x = ax / s, r1y = ay / s;
			float r2x = s / ax, r2y = s / ay;
			Vector2 r1 = a / s;
			Vector2 r2 = s / a;

			if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y)))
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
					printf("a = (%.5f, %.5f), b = (%.5f, %.5f)", a.x, a.y, b.x, b.y);
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
		float ax = a.x, ay = a.y;
		float bx = b.x, by = b.y;
		float rx = ax + bx, ry = ay + by;
		Vector2 r = a + b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y);
	}
};

class Vector2Subtraction
{
public:
	static const char* GetName() { return "Subtraction"; }
	bool operator()(const Vector2& a, const Vector2& b) const
	{
		float ax = a.x, ay = a.y;
		float bx = b.x, by = b.y;
		float rx = ax - bx, ry = ay - by;
		Vector2 r = a - b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y);
	}
};

class Vector2DotProduct
{
public:
	static const char* GetName() { return "Dot Product"; }
	bool operator()(const Vector2& a, const Vector2& b) const
	{
		float ax = a.x, ay = a.y;
		float bx = b.x, by = b.y;
		float r = (ax * bx) + (ay * by);
		float dot = Dot(a, b);
		return FloatTest(dot, r);
	}
};

class Vector2CompMul
{
public:
	static const char* GetName() { return "Comp Mul"; }
	bool operator()(const Vector2& a, const Vector2& b) const
	{
		float ax = a.x, ay = a.y;
		float bx = b.x, by = b.y;
		float rx = (ax * bx), ry = (ay * by);
		Vector2 r = a * b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y);
	}
};

class Vector2CompDiv
{
public:
	static const char* GetName() { return "Comp Div"; }
	bool operator()(const Vector2& a, const Vector2& b) const
	{
		float ax = a.x, ay = a.y;
		float bx = b.x, by = b.y;
		float rx = (ax / bx), ry = (ay / by);
		Vector2 r = a / b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y);
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
				printf("a = (%.5f, %.5f, %.5f)", a.x, a.y, a.z);
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
		float ax = a.x, ay = a.y, az = a.z;
		float rx = -ax, ry = -ay, rz = -az;
		Vector3 r = -a;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
	}
};

class Vector3Length
{
public:
	static const char* GetName() { return "Length"; }
	bool operator() (const Vector3& a)
	{
		float ax = a.x, ay = a.y, az = a.z;
		float r = sqrt((ax * ax) + (ay * ay) + (az * az));
		float len = a.Len();
		return FloatTest(r, len);
	}
};

class Vector3UnitVec
{
public:
	static const char* GetName() { return "UnitVec"; }
	bool operator() (const Vector3& a)
	{
		float ax = a.x, ay = a.y, az = a.z;
		float len = sqrt((ax * ax) + (ay * ay) + (az * az));
		float rx = ax / len, ry = ay / len, rz = az / len;
		Vector3 r = a.Unit();
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz ,r.z);
	}
};

class Vector3ScalarMultiplication
{
public:
	static const char* GetName() { return "Scalar Multiplication"; }
	bool operator() (const Vector3& a)
	{
		float ax = a.x, ay = a.y, az = a.z;
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector * scalar) and (scalar * vector)
			float s = s_floatVec[i];
			float r1x = ax * s, r1y = ay * s, r1z = az * s;
			float r2x = s * ax, r2y = s * ay, r2z = s * az;
			Vector3 r1 = a * s;
			Vector3 r2 = s * a;
			if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) && 
				  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z) &&
				  FloatTest(r1.x, r2.x) && FloatTest(r1.y, r2.y) && FloatTest(r1.z, r2.z)))
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
		float ax = a.x, ay = a.y, az = a.z;
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector / scalar) and (scalar / vector)
			float s = s_floatVec[i];
			float r1x = ax / s, r1y = ay / s, r1z = az / s;
			float r2x = s / ax, r2y = s / ay, r2z = s / az;
			Vector3 r1 = a / s;
			Vector3 r2 = s / a;

			if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) && 
				  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z)))
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
					printf("a = (%.5f, %.5f, %.5f), b = (%.5f, %.5f, %.5f)", a.x, a.y, a.z, b.x, b.y, b.z);
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
		float ax = a.x, ay = a.y, az = a.z;
		float bx = b.x, by = b.y, bz = b.z;
		float rx = ax + bx, ry = ay + by, rz = az + bz;
		Vector3 r = a + b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
	}
};

class Vector3Subtraction
{
public:
	static const char* GetName() { return "Subtraction"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.x, ay = a.y, az = a.z;
		float bx = b.x, by = b.y, bz = b.z;
		float rx = ax - bx, ry = ay - by, rz = az - bz;
		Vector3 r = a - b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
	}
};

class Vector3DotProduct
{
public:
	static const char* GetName() { return "Dot Product"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.x, ay = a.y, az = a.z;
		float bx = b.x, by = b.y, bz = b.z;
		float r = (ax * bx) + (ay * by) + (az * bz);
		float dot = Dot(a, b);
		return FloatTest(dot, r);
	}
};

class Vector3CompMul
{
public:
	static const char* GetName() { return "Comp Mul"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.x, ay = a.y, az = a.z;
		float bx = b.x, by = b.y, bz = b.z;
		float rx = (ax * bx), ry = (ay * by), rz = (az * bz);
		Vector3 r = a * b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
	}
};

class Vector3CompDiv
{
public:
	static const char* GetName() { return "Comp Div"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.x, ay = a.y, az = a.z;
		float bx = b.x, by = b.y, bz = b.z;
		float rx = (ax / bx), ry = (ay / by), rz = (az / bz);
		Vector3 r = a / b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
	}
};

class Vector3CrossProduct
{
public:
	static const char* GetName() { return "Cross Product"; }
	bool operator()(const Vector3& a, const Vector3& b) const
	{
		float ax = a.x, ay = a.y, az = a.z;
		float bx = b.x, by = b.y, bz = b.z;

		float rx = (ay * bz) - (az * by);
		float ry = (az * bx) - (ax * bz);
		float rz = (ax * by) - (ay * bx);

		Vector3 r = Cross(a, b);
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z);
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
				printf("a = (%.5f, %.5f, %.5f, %.5f)", a.x, a.y, a.z, a.w);
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
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float rx = -ax, ry = -ay, rz = -az, rw = -aw;
		Vector4 r = -a;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
	}
};

class Vector4Length
{
public:
	static const char* GetName() { return "Length"; }
	bool operator() (const Vector4& a)
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float r = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
		float len = a.Len();
		return FloatTest(r, len);
	}
};

class Vector4UnitVec
{
public:
	static const char* GetName() { return "UnitVec"; }
	bool operator() (const Vector4& a)
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float len = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
		float rx = ax / len, ry = ay / len, rz = az / len, rw = aw / len;
		Vector4 r = a.Unit();
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
	}
};

class Vector4ScalarMultiplication
{
public:
	static const char* GetName() { return "Scalar Multiplication"; }
	bool operator() (const Vector4& a)
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector * scalar) and (scalar * vector)
			float s = s_floatVec[i];
			float r1x = ax * s, r1y = ay * s, r1z = az * s, r1w = aw * s;
			float r2x = s * ax, r2y = s * ay, r2z = s * az, r2w = s * aw;
			Vector4 r1 = a * s;
			Vector4 r2 = s * a;
			if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) && FloatTest(r1w, r1.w) &&
				  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z) && FloatTest(r2w, r2.w) &&
				  FloatTest(r1.x, r2.x) && FloatTest(r1.y, r2.y) && FloatTest(r1.z, r2.z) && FloatTest(r1.w, r2.w)))
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
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector / scalar) and (scalar / vector)
			float s = s_floatVec[i];
			float r1x = ax / s, r1y = ay / s, r1z = az / s, r1w = aw / s;
			float r2x = s / ax, r2y = s / ay, r2z = s / az, r2w = s / aw;
			Vector4 r1 = a / s;
			Vector4 r2 = s / a;

			if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) && FloatTest(r1w, r1.w) &&
				  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z) && FloatTest(r2w, r2.w)))
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
					printf("a = (%.5f, %.5f, %.5f, %.5f), b = (%.5f, %.5f, %.5f, %.5f)", a.x, a.y, a.z, a.w, b.x, b.y, b.z, b.w);
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
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float bx = b.x, by = b.y, bz = b.z, bw = b.w;
		float rx = ax + bx, ry = ay + by, rz = az + bz, rw = aw + bw;
		Vector4 r = a + b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
	}
};

class Vector4Subtraction
{
public:
	static const char* GetName() { return "Subtraction"; }
	bool operator()(const Vector4& a, const Vector4& b) const
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float bx = b.x, by = b.y, bz = b.z, bw = b.w;
		float rx = ax - bx, ry = ay - by, rz = az - bz, rw = aw - bw;
		Vector4 r = a - b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
	}
};

class Vector4DotProduct
{
public:
	static const char* GetName() { return "Dot Product"; }
	bool operator()(const Vector4& a, const Vector4& b) const
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float bx = b.x, by = b.y, bz = b.z, bw = b.w;
		float r = (ax * bx) + (ay * by) + (az * bz) + (aw * bw);
		float dot = Dot(a, b);
		return FloatTest(dot, r);
	}
};

class Vector4CompMul
{
public:
	static const char* GetName() { return "Comp Mul"; }
	bool operator()(const Vector4& a, const Vector4& b) const
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float bx = b.x, by = b.y, bz = b.z, bw = b.w;
		float rx = (ax * bx), ry = (ay * by), rz = (az * bz), rw = (aw * bw);
		Vector4 r = a * b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
	}
};

class Vector4CompDiv
{
public:
	static const char* GetName() { return "Comp Div"; }
	bool operator()(const Vector4& a, const Vector4& b) const
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float bx = b.x, by = b.y, bz = b.z, bw = b.w;
		float rx = (ax / bx), ry = (ay / by), rz = (az / bz), rw = (aw / bw);
		Vector4 r = a / b;
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
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
				printf("a = (%.5f, %.5f, %.5f, %.5f)", a.i, a.j, a.k, a.r);
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
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		float rx, ry, rz, rw;

		for (unsigned int i = 0; i < s_vector3Vec.size(); ++i)
		{
			Vector3 v = s_vector3Vec[i];
			float bx = v.x, by = v.y, bz = v.z, bw = 0.0f;

			QuatMul(&rx, &ry, &rz, &rw, ax, ay, az, aw, bx, by, bz, bw);
			QuatMul(&rx, &ry, &rz, &rw, rx, ry, rz, rw, -ax, -ay, -az, aw);

			Vector3 r = a.Rotate(v);
			if (!(FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z)))
			{
				printf("v = (%.5f, %.5f, %.5f)", v.x, v.y, v.z);
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
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		float rx = -ax, ry = -ay, rz = -az, rw = aw;
		Quat r = ~a;
		return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
	}
};

class QuatNegation
{
public:
	static const char* GetName() { return "Negation"; }
	bool operator() (const Quat& a)
	{
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		float rx = -ax, ry = -ay, rz = -az, rw = -aw;
		Quat r = -a;
		return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
	}
};

/*
class QuatScalarMultiplication
{
public:
	static const char* GetName() { return "Scalar Multiplication"; }
	bool operator() (const Quat& a)
	{
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector * scalar) and (scalar * vector)
			float s = s_floatVec[i];
			float r1x = ax * s, r1y = ay * s, r1z = az * s, r1w = aw * s;
			float r2x = s * ax, r2y = s * ay, r2z = s * az, r2w = s * aw;
			Quat r1 = a * s;
			Quat r2 = s * a;
			if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) && FloatTest(r1w, r1.w) &&
				  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z) && FloatTest(r2w, r2.w) &&
				  FloatTest(r1.x, r2.x) && FloatTest(r1.y, r2.y) && FloatTest(r1.z, r2.z) && FloatTest(r1.w, r2.w)))
			{
				printf("s = %.5f", s);
				return false;
			}
		}
		return true;
	}
};
*/

/*
class QuatScalarDivision
{
public:
	static const char* GetName() { return "Scalar Division"; }
	bool operator() (const Quat& a)
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		for (unsigned int i = 0; i < s_floatVec.size(); ++i)
		{
			// test (vector / scalar) and (scalar / vector)
			float s = s_floatVec[i];
			float r1x = ax / s, r1y = ay / s, r1z = az / s, r1w = aw / s;
			float r2x = s / ax, r2y = s / ay, r2z = s / az, r2w = s / aw;
			Quat r1 = a / s;
			Quat r2 = s / a;

			if (!(FloatTest(r1x, r1.x) && FloatTest(r1y, r1.y) && FloatTest(r1z, r1.z) && FloatTest(r1w, r1.w) &&
				  FloatTest(r2x, r2.x) && FloatTest(r2y, r2.y) && FloatTest(r2z, r2.z) && FloatTest(r2w, r2.w)))
			{
				printf("s = %.5f\n", s);
				return false;
			}
		}
		return true;
	}
};
*/

class QuatLength
{
public:
	static const char* GetName() { return "Length"; }
	bool operator() (const Quat& a)
	{
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		float r = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
		float len = a.Len();
		return FloatTest(r, len);
	}
};

class QuatUnitVec
{
public:
	static const char* GetName() { return "UnitVec"; }
	bool operator() (const Quat& a)
	{
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		float len = sqrt((ax * ax) + (ay * ay) + (az * az) + (aw * aw));
		float rx = ax / len, ry = ay / len, rz = az / len, rw = aw / len;
		Quat r = a.Unit();
		return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
	}
};

void QuatExp(float* rx, float* ry, float* rz, float* rw, float qx, float qy, float qz, float qw)
{
	float angle = sqrt((qx * qx) + (qy * qy) + (qz * qz));
	float sin_a = sin(angle / 2.0f);
	// prevent div by zero
	if (angle > 0.0001f)
	{
		*rx = (qx / angle) * sin_a;
		*ry = (qy / angle) * sin_a;
		*rz = (qz / angle) * sin_a;
	}
	else
		*rx = *ry = *rz = 0.0f;

	*rw = cos(angle / 2.0f);
}

class QuatExponential
{
public:
	static const char* GetName() { return "Quat Exponential"; }
	bool operator() (const Quat& a)
	{
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		float rx, ry, rz, rw;
		QuatExp(&rx, &ry, &rz, &rw, ax, ay, az, aw);
		Quat r = a.Exp();
		bool rval = FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
		if (!rval)
		{
			printf("rx = %.5f, ry = %.5f, rz = %.5f, rw = %.5f\n", rx, ry, rz, rw);
			printf("r = %.5f, %.5f, %.5f, %.5f\n", r.i, r.j, r.k, r.r);
		}
		return rval;
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
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		float rx, ry, rz, rw;
		QuatLog(&rx, &ry, &rz, &rw, ax, ay, az, aw);
		Quat r = a.Log();
		return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
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
					printf("a = (%.5f, %.5f, %.5f, %.5f), b = (%.5f, %.5f, %.5f, %.5f)", a.i, a.j, a.k, a.r, b.i, b.j, b.k, b.r);
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
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		float bx = b.i, by = b.j, bz = b.k, bw = b.r;
		float rx = ax + bx, ry = ay + by, rz = az + bz, rw = aw + bw;
		Quat r = a + b;
		return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
	}
};

class QuatSubtraction
{
public:
	static const char* GetName() { return "Subtraction"; }
	bool operator()(const Quat& a, const Quat& b) const
	{
		float ax = a.i, ay = a.j, az = a.k, aw = a.r;
		float bx = b.i, by = b.j, bz = b.k, bw = b.r;
		float rx = ax - bx, ry = ay - by, rz = az - bz, rw = aw - bw;
		Quat r = a - b;
		return FloatTest(rx, r.i) && FloatTest(ry, r.j) && FloatTest(rz, r.k) && FloatTest(rw, r.r);
	}
};

/*
class QuatCompMul
{
public:
	static const char* GetName() { return "Comp Mul"; }
	bool operator()(const Quat& a, const Quat& b) const
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float bx = b.x, by = b.y, bz = b.z, bw = b.w;
		float rx = (ax * bx), ry = (ay * by), rz = (az * bz), rw = (aw * bw);
		Quat r = CompMul(a, b);
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
	}
};
*/
/*
class QuatCompDiv
{
public:
	static const char* GetName() { return "Comp Div"; }
	bool operator()(const Quat& a, const Quat& b) const
	{
		float ax = a.x, ay = a.y, az = a.z, aw = a.w;
		float bx = b.x, by = b.y, bz = b.z, bw = b.w;
		float rx = (ax / bx), ry = (ay / by), rz = (az / bz), rw = (aw / bw);
		Quat r = CompDiv(a, b);
		return FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w);
	}
};
*/

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
				printf("m = \n");
				PrintMatrix(m);
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
			float bx = v.x, by = v.y, bz = v.z;

			float rx = fv[0] * bx + fv[1] * by + fv[2] * bz;
			float ry = fv[4] * bx + fv[5] * by + fv[6] * bz;
			float rz = fv[8] * bx + fv[9] * by + fv[10] * bz;

			Vector3 r = Transform3x3(a, v);
			if (!(FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z)))
			{
				printf("v = (%.5f, %.5f, %.5f)", v.x, v.y, v.z);
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
			float bx = v.x, by = v.y, bz = v.z;

			float rx = fv[0] * bx + fv[1] * by + fv[2] * bz + fv[3];
			float ry = fv[4] * bx + fv[5] * by + fv[6] * bz + fv[7];
			float rz = fv[8] * bx + fv[9] * by + fv[10] * bz + fv[11];

			Vector3 r = Transform3x4(a, v);
			if (!(FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z)))
			{
				printf("v = (%.5f, %.5f, %.5f)", v.x, v.y, v.z);
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
			float bx = v.x, by = v.y, bz = v.z, bw = v.w;

			float rx = fv[0] * bx + fv[1] * by + fv[2] * bz + fv[3] * bw;
			float ry = fv[4] * bx + fv[5] * by + fv[6] * bz + fv[7] * bw;
			float rz = fv[8] * bx + fv[9] * by + fv[10] * bz + fv[11] * bw;
			float rw = fv[12] * bx + fv[13] * by + fv[14] * bz + fv[15] * bw;

			Vector4 r = Transform4x4(a, v);
			if (!(FloatTest(rx, r.x) && FloatTest(ry, r.y) && FloatTest(rz, r.z) && FloatTest(rw, r.w)))
			{
				printf("v = (%.5f, %.5f, %.5f, %.5f)", v.x, v.y, v.z, v.w);
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
					printf("a = \n");
					PrintMatrix(a);

					printf("b = \n");
					PrintMatrix(b);
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


/*
 * Compute inverse of 4x4 transformation matrix.
 * Code contributed by Jacques Leroy jle@star.be
 * Return GL_TRUE for success, GL_FALSE for failure (singular matrix)
 */
static bool
invert_matrix(const float * m, float * out)
{
/* NB. OpenGL Matrices are COLUMN major. */
#define SWAP_ROWS(a, b) { float *_tmp = a; (a)=(b); (b)=_tmp; }
#define MAT(m,r,c) (m)[(c)*4+(r)]

   float wtmp[4][8];
   float m0, m1, m2, m3, s;
   float *r0, *r1, *r2, *r3;

   r0 = wtmp[0], r1 = wtmp[1], r2 = wtmp[2], r3 = wtmp[3];

   r0[0] = MAT(m, 0, 0), r0[1] = MAT(m, 0, 1),
      r0[2] = MAT(m, 0, 2), r0[3] = MAT(m, 0, 3),
      r0[4] = 1.0, r0[5] = r0[6] = r0[7] = 0.0,
      r1[0] = MAT(m, 1, 0), r1[1] = MAT(m, 1, 1),
      r1[2] = MAT(m, 1, 2), r1[3] = MAT(m, 1, 3),
      r1[5] = 1.0, r1[4] = r1[6] = r1[7] = 0.0,
      r2[0] = MAT(m, 2, 0), r2[1] = MAT(m, 2, 1),
      r2[2] = MAT(m, 2, 2), r2[3] = MAT(m, 2, 3),
      r2[6] = 1.0, r2[4] = r2[5] = r2[7] = 0.0,
      r3[0] = MAT(m, 3, 0), r3[1] = MAT(m, 3, 1),
      r3[2] = MAT(m, 3, 2), r3[3] = MAT(m, 3, 3),
      r3[7] = 1.0, r3[4] = r3[5] = r3[6] = 0.0;

   /* choose pivot - or die */
   if (fabs(r3[0]) > fabs(r2[0]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[0]) > fabs(r1[0]))
      SWAP_ROWS(r2, r1);
   if (fabs(r1[0]) > fabs(r0[0]))
      SWAP_ROWS(r1, r0);
   if (0.0 == r0[0])
      return false;

   /* eliminate first variable     */
   m1 = r1[0] / r0[0];
   m2 = r2[0] / r0[0];
   m3 = r3[0] / r0[0];
   s = r0[1];
   r1[1] -= m1 * s;
   r2[1] -= m2 * s;
   r3[1] -= m3 * s;
   s = r0[2];
   r1[2] -= m1 * s;
   r2[2] -= m2 * s;
   r3[2] -= m3 * s;
   s = r0[3];
   r1[3] -= m1 * s;
   r2[3] -= m2 * s;
   r3[3] -= m3 * s;
   s = r0[4];
   if (s != 0.0) {
      r1[4] -= m1 * s;
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r0[5];
   if (s != 0.0) {
      r1[5] -= m1 * s;
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r0[6];
   if (s != 0.0) {
      r1[6] -= m1 * s;
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r0[7];
   if (s != 0.0) {
      r1[7] -= m1 * s;
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[1]) > fabs(r2[1]))
      SWAP_ROWS(r3, r2);
   if (fabs(r2[1]) > fabs(r1[1]))
      SWAP_ROWS(r2, r1);
   if (0.0 == r1[1])
	   return false;

   /* eliminate second variable */
   m2 = r2[1] / r1[1];
   m3 = r3[1] / r1[1];
   r2[2] -= m2 * r1[2];
   r3[2] -= m3 * r1[2];
   r2[3] -= m2 * r1[3];
   r3[3] -= m3 * r1[3];
   s = r1[4];
   if (0.0 != s) {
      r2[4] -= m2 * s;
      r3[4] -= m3 * s;
   }
   s = r1[5];
   if (0.0 != s) {
      r2[5] -= m2 * s;
      r3[5] -= m3 * s;
   }
   s = r1[6];
   if (0.0 != s) {
      r2[6] -= m2 * s;
      r3[6] -= m3 * s;
   }
   s = r1[7];
   if (0.0 != s) {
      r2[7] -= m2 * s;
      r3[7] -= m3 * s;
   }

   /* choose pivot - or die */
   if (fabs(r3[2]) > fabs(r2[2]))
      SWAP_ROWS(r3, r2);
   if (0.0 == r2[2])
      return false;

   /* eliminate third variable */
   m3 = r3[2] / r2[2];
   r3[3] -= m3 * r2[3], r3[4] -= m3 * r2[4],
      r3[5] -= m3 * r2[5], r3[6] -= m3 * r2[6], r3[7] -= m3 * r2[7];

   /* last check */
   if (0.0 == r3[3])
      return false;

   s = 1.0 / r3[3];		/* now back substitute row 3 */
   r3[4] *= s;
   r3[5] *= s;
   r3[6] *= s;
   r3[7] *= s;

   m2 = r2[3];			/* now back substitute row 2 */
   s = 1.0 / r2[2];
   r2[4] = s * (r2[4] - r3[4] * m2), r2[5] = s * (r2[5] - r3[5] * m2),
      r2[6] = s * (r2[6] - r3[6] * m2), r2[7] = s * (r2[7] - r3[7] * m2);
   m1 = r1[3];
   r1[4] -= r3[4] * m1, r1[5] -= r3[5] * m1,
      r1[6] -= r3[6] * m1, r1[7] -= r3[7] * m1;
   m0 = r0[3];
   r0[4] -= r3[4] * m0, r0[5] -= r3[5] * m0,
      r0[6] -= r3[6] * m0, r0[7] -= r3[7] * m0;

   m1 = r1[2];			/* now back substitute row 1 */
   s = 1.0 / r1[1];
   r1[4] = s * (r1[4] - r2[4] * m1), r1[5] = s * (r1[5] - r2[5] * m1),
      r1[6] = s * (r1[6] - r2[6] * m1), r1[7] = s * (r1[7] - r2[7] * m1);
   m0 = r0[2];
   r0[4] -= r2[4] * m0, r0[5] -= r2[5] * m0,
      r0[6] -= r2[6] * m0, r0[7] -= r2[7] * m0;

   m0 = r0[1];			/* now back substitute row 0 */
   s = 1.0 / r0[0];
   r0[4] = s * (r0[4] - r1[4] * m0), r0[5] = s * (r0[5] - r1[5] * m0),
      r0[6] = s * (r0[6] - r1[6] * m0), r0[7] = s * (r0[7] - r1[7] * m0);

   MAT(out, 0, 0) = r0[4];
   MAT(out, 0, 1) = r0[5], MAT(out, 0, 2) = r0[6];
   MAT(out, 0, 3) = r0[7], MAT(out, 1, 0) = r1[4];
   MAT(out, 1, 1) = r1[5], MAT(out, 1, 2) = r1[6];
   MAT(out, 1, 3) = r1[7], MAT(out, 2, 0) = r2[4];
   MAT(out, 2, 1) = r2[5], MAT(out, 2, 2) = r2[6];
   MAT(out, 2, 3) = r2[7], MAT(out, 3, 0) = r3[4];
   MAT(out, 3, 1) = r3[5], MAT(out, 3, 2) = r3[6];
   MAT(out, 3, 3) = r3[7];

   return true;

#undef MAT
#undef SWAP_ROWS
}

/*
 * Compute the inverse of a 4x4 matrix.
 *
 * From an algorithm by V. Strassen, 1969, _Numerishe Mathematik_, vol. 13,
 * pp. 354-356.
 * 60 multiplies, 24 additions, 10 subtractions, 8 negations, 2 divisions,
 * 48 assignments, _0_ branches
 *
 * This implementation by Scott McCaskill
 */

typedef float Mat2[2][2];

enum {
    M00 = 0, M01 = 4, M02 = 8, M03 = 12,
    M10 = 1, M11 = 5, M12 = 9, M13 = 13,
    M20 = 2, M21 = 6, M22 = 10,M23 = 14,
    M30 = 3, M31 = 7, M32 = 11,M33 = 15
};

static void invert_matrix_general( const float *m, float *out )
{
   Mat2 r1, r2, r3, r4, r5, r6, r7;
   const float * A = m;
   float *       C = out;
   float one_over_det;

   /*
    * A is the 4x4 source matrix (to be inverted).
    * C is the 4x4 destination matrix
    * a11 is the 2x2 matrix in the upper left quadrant of A
    * a12 is the 2x2 matrix in the upper right quadrant of A
    * a21 is the 2x2 matrix in the lower left quadrant of A
    * a22 is the 2x2 matrix in the lower right quadrant of A
    * similarly, cXX are the 2x2 quadrants of the destination matrix
    */

   /* R1 = inverse( a11 ) */
   one_over_det = 1.0f / ( ( A[M00] * A[M11] ) - ( A[M10] * A[M01] ) );
   r1[0][0] = one_over_det * A[M11];
   r1[0][1] = one_over_det * -A[M01];
   r1[1][0] = one_over_det * -A[M10];
   r1[1][1] = one_over_det * A[M00];

   /* R2 = a21 x R1 */
   r2[0][0] = A[M20] * r1[0][0] + A[M21] * r1[1][0];
   r2[0][1] = A[M20] * r1[0][1] + A[M21] * r1[1][1];
   r2[1][0] = A[M30] * r1[0][0] + A[M31] * r1[1][0];
   r2[1][1] = A[M30] * r1[0][1] + A[M31] * r1[1][1];

   /* R3 = R1 x a12 */
   r3[0][0] = r1[0][0] * A[M02] + r1[0][1] * A[M12];
   r3[0][1] = r1[0][0] * A[M03] + r1[0][1] * A[M13];
   r3[1][0] = r1[1][0] * A[M02] + r1[1][1] * A[M12];
   r3[1][1] = r1[1][0] * A[M03] + r1[1][1] * A[M13];

   /* R4 = a21 x R3 */
   r4[0][0] = A[M20] * r3[0][0] + A[M21] * r3[1][0];
   r4[0][1] = A[M20] * r3[0][1] + A[M21] * r3[1][1];
   r4[1][0] = A[M30] * r3[0][0] + A[M31] * r3[1][0];
   r4[1][1] = A[M30] * r3[0][1] + A[M31] * r3[1][1];

   /* R5 = R4 - a22 */
   r5[0][0] = r4[0][0] - A[M22];
   r5[0][1] = r4[0][1] - A[M23];
   r5[1][0] = r4[1][0] - A[M32];
   r5[1][1] = r4[1][1] - A[M33];

   /* R6 = inverse( R5 ) */
   one_over_det = 1.0f / ( ( r5[0][0] * r5[1][1] ) - ( r5[1][0] * r5[0][1] ) );
   r6[0][0] = one_over_det * r5[1][1];
   r6[0][1] = one_over_det * -r5[0][1];
   r6[1][0] = one_over_det * -r5[1][0];
   r6[1][1] = one_over_det * r5[0][0];

   /* c12 = R3 x R6 */
   C[M02] = r3[0][0] * r6[0][0] + r3[0][1] * r6[1][0];
   C[M03] = r3[0][0] * r6[0][1] + r3[0][1] * r6[1][1];
   C[M12] = r3[1][0] * r6[0][0] + r3[1][1] * r6[1][0];
   C[M13] = r3[1][0] * r6[0][1] + r3[1][1] * r6[1][1];

   /* c21 = R6 x R2 */
   C[M20] = r6[0][0] * r2[0][0] + r6[0][1] * r2[1][0];
   C[M21] = r6[0][0] * r2[0][1] + r6[0][1] * r2[1][1];
   C[M30] = r6[1][0] * r2[0][0] + r6[1][1] * r2[1][0];
   C[M31] = r6[1][0] * r2[0][1] + r6[1][1] * r2[1][1];

   /* R7 = R3 x c21 */
   r7[0][0] = r3[0][0] * C[M20] + r3[0][1] * C[M30];
   r7[0][1] = r3[0][0] * C[M21] + r3[0][1] * C[M31];
   r7[1][0] = r3[1][0] * C[M20] + r3[1][1] * C[M30];
   r7[1][1] = r3[1][0] * C[M21] + r3[1][1] * C[M31];

   /* c11 = R1 - R7 */
   C[M00] = r1[0][0] - r7[0][0];
   C[M01] = r1[0][1] - r7[0][1];
   C[M10] = r1[1][0] - r7[1][0];
   C[M11] = r1[1][1] - r7[1][1];

   /* c22 = -R6 */
   C[M22] = -r6[0][0];
   C[M23] = -r6[0][1];
   C[M32] = -r6[1][0];
   C[M33] = -r6[1][1];
}




class MatrixInverse
{
public:
	static const char* GetName() { return "Inverse"; }
	bool operator()(const Matrix& a) const
	{
//		Matrix m(Quat(Vector3(0.3f,1.5f,0.7f), -PI/2.0f), Vector3(1.0, 2.0f, 3.0f));
//		m.SetScale(Vector3(1,2,1));

		float fv[16];
		
		// convert the transpose of a into a float vec. (why? because invert_matrix expects an opengl style matrix)
		MatrixToFloatVec(fv, Transpose(a));
		float rfv[16];
		bool invertable = invert_matrix(fv, rfv);

		float rfv2[16];
		invert_matrix_general(fv, rfv2);
		
		Matrix r;
		bool matrix_invertable = Inverse(a, r);

		if (invertable != matrix_invertable)
		{
			printf("Failure because invertable mismatch\n");
			return false;  // test failed
		}

		if (invertable)
		{
			bool rval;
			// compare r with the transpose of rfv.  (why? becaaue invert_matrix returns an opengl style matrix)
			const float epsilon = 0.05f;
			rval = FuzzyFloatTest(rfv[0], r.Elem(0,0), epsilon) && 
				FuzzyFloatTest(rfv[4], r.Elem(0,1), epsilon) && 
				FuzzyFloatTest(rfv[8], r.Elem(0,2), epsilon) && 
				FuzzyFloatTest(rfv[12], r.Elem(0,3), epsilon) &&
				FuzzyFloatTest(rfv[1], r.Elem(1,0), epsilon) && 
				FuzzyFloatTest(rfv[5], r.Elem(1,1), epsilon) && 
				FuzzyFloatTest(rfv[9], r.Elem(1,2), epsilon) && 
				FuzzyFloatTest(rfv[13], r.Elem(1,3), epsilon) &&
				FuzzyFloatTest(rfv[2], r.Elem(2,0), epsilon) && 
				FuzzyFloatTest(rfv[6], r.Elem(2,1), epsilon) && 
				FuzzyFloatTest(rfv[10], r.Elem(2,2), epsilon) && 
				FuzzyFloatTest(rfv[14], r.Elem(2,3), epsilon) &&
				FuzzyFloatTest(rfv[3], r.Elem(3,0), epsilon) && 
				FuzzyFloatTest(rfv[7], r.Elem(3,1), epsilon) && 
				FuzzyFloatTest(rfv[11], r.Elem(3,2), epsilon) && 
				FuzzyFloatTest(rfv[15], r.Elem(3,3), epsilon);

			if (!rval)
			{
				printf("Failure because matrices are not equal\n");

				printf("rfv = \n");
				printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv[0], rfv[4], rfv[8], rfv[12]);
				printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv[1], rfv[5], rfv[9], rfv[13]);
				printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv[2], rfv[6], rfv[10], rfv[14]);
				printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv[3], rfv[7], rfv[11], rfv[15]);

				printf("rfv2 = \n");
				printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv2[0], rfv2[4], rfv2[8], rfv2[12]);
				printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv2[1], rfv2[5], rfv2[9], rfv2[13]);
				printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv2[2], rfv2[6], rfv2[10], rfv2[14]);
				printf("| %10.5f, %10.5f, %10.5f, %10.5f |\n", rfv2[3], rfv2[7], rfv2[11], rfv2[15]);


				printf("r = \n");
				PrintMatrix(r);
			}

			return rval;
		}
		else
		{
			return true;  // both matrices were determined to be non-invertable.
		}
	}
};

int main(int argc, char* argv[])
{
	InitTestData();

	TestSuite floatSuite("float");
	floatSuite.AddTest(new FloatUnaryOpTest<FloatDegToRad>());
	floatSuite.AddTest(new FloatUnaryOpTest<FloatRadToDeg>());
	floatSuite.AddTest(new FloatTernaryOpTest<FloatClamp>());
	floatSuite.AddTest(new FloatUnaryOpTest<FloatLimitPi>());
	floatSuite.AddTest(new FloatUnaryOpTest<FloatMod2Pi>());
	floatSuite.AddTest(new FloatBinaryOpTest<FloatFuzzyEqual>());
	floatSuite.RunTests();

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
//	quatSuite.AddTest(new QuatUnaryOpTest<QuatScalarMultiplication>());
//	quatSuite.AddTest(new QuatUnaryOpTest<QuatScalarDivision>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatLength>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatUnitVec>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatExponential>());
	quatSuite.AddTest(new QuatUnaryOpTest<QuatLogarithm>());
	quatSuite.AddTest(new QuatBinaryOpTest<QuatAddition>());
	quatSuite.AddTest(new QuatBinaryOpTest<QuatSubtraction>());
//	quatSuite.AddTest(new QuatBinaryOpTest<QuatCompMul>());
//	quatSuite.AddTest(new QuatBinaryOpTest<QuatCompDiv>());
	quatSuite.RunTests();

	TestSuite matrixSuite("Matrix");
	matrixSuite.AddTest(new MatrixUnaryOpTest<MatrixTransform3x3>());
	matrixSuite.AddTest(new MatrixUnaryOpTest<MatrixTransform3x4>());
	matrixSuite.AddTest(new MatrixUnaryOpTest<MatrixTransform4x4>());
	matrixSuite.AddTest(new MatrixBinaryOpTest<MatrixMultiplication>());
	matrixSuite.AddTest(new MatrixUnaryOpTest<MatrixInverse>());
	
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
	matrixSuite.RunTests();

	// TODO: Rad2Deg
	// TODO: Deg2Rad

	return 0;
}
