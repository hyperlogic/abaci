#include "math3d.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "unittest.h"

#include <CoreServices/CoreServices.h>
#include <mach/mach.h>
#include <mach/mach_time.h>

// static test containers
std::vector<Vector2> s_vector2Vec;
std::vector<Vector3> s_vector3Vec;
std::vector<Vector4> s_vector4Vec;
std::vector<Matrix> s_matrixVec;

float s_floatResult;
Vector4 s_vector4Result;

// 0..1
float RandomNormalizedFloat()
{
	return (float)rand() / RAND_MAX;
}

// -range to range
float RandomFloat(float range = 10000.0f)
{
	return (((float)rand() / (RAND_MAX / 2)) - 1.0f) * range;
}

void InitData()
{
	srand((unsigned int)mach_absolute_time());

	// add random values
	const int RANDOM_DATA_SIZE = 64;
	for (int i = 0; i < RANDOM_DATA_SIZE; ++i)
	{
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
	s_vector2Vec.push_back(Vector2(0, 0));
	s_vector3Vec.push_back(Vector3(0, 0, 0));
	s_vector4Vec.push_back(Vector4(0, 0, 0, 0));
	s_matrixVec.push_back(Matrix(Vector4(0, 0, 0, 0), Vector4(0, 0, 0, 0), Vector4(0, 0, 0, 0), Vector4(0, 0, 0, 0)));

	// ident
	s_matrixVec.push_back(Matrix(Vector4(1, 0, 0, 0), Vector4(0, 1, 0, 0), Vector4(0, 0, 1, 0), Vector4(0, 0, 0, 1)));

	// FLT_MAX
	s_vector2Vec.push_back(Vector2(FLT_MAX, FLT_MAX));
	s_vector3Vec.push_back(Vector3(FLT_MAX, FLT_MAX, FLT_MAX));
	s_vector4Vec.push_back(Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX));
	s_matrixVec.push_back(Matrix(Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX), 
								 Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX), 
								 Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX), 
								 Vector4(FLT_MAX, FLT_MAX, FLT_MAX, FLT_MAX)));

	s_vector2Vec.push_back(Vector2(-FLT_MAX, -FLT_MAX));
	s_vector3Vec.push_back(Vector3(-FLT_MAX, -FLT_MAX, -FLT_MAX));
	s_vector4Vec.push_back(Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX));
	s_matrixVec.push_back(Matrix(Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX), 
								 Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX), 
								 Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX), 
								 Vector4(-FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX)));

	// FLT_MIN
	s_vector2Vec.push_back(Vector2(FLT_MIN, FLT_MIN));
	s_vector3Vec.push_back(Vector3(FLT_MIN, FLT_MIN, FLT_MIN));
	s_vector4Vec.push_back(Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN));
	s_matrixVec.push_back(Matrix(Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN), 
								 Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN), 
								 Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN), 
								 Vector4(FLT_MIN, FLT_MIN, FLT_MIN, FLT_MIN)));

	s_vector2Vec.push_back(Vector2(-FLT_MIN, -FLT_MIN));
	s_vector3Vec.push_back(Vector3(-FLT_MIN, -FLT_MIN, -FLT_MIN));
	s_vector4Vec.push_back(Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN));
	s_matrixVec.push_back(Matrix(Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN), 
								 Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN), 
								 Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN), 
								 Vector4(-FLT_MIN, -FLT_MIN, -FLT_MIN, -FLT_MIN)));

	// inf
	float inf = FLT_MAX * FLT_MAX;

	s_vector2Vec.push_back(Vector2(inf, inf));
	s_vector3Vec.push_back(Vector3(inf, inf, inf));
	s_vector4Vec.push_back(Vector4(inf, inf, inf, inf));
	s_matrixVec.push_back(Matrix(Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf)));

	// -inf
	inf = -FLT_MAX * FLT_MAX;
	s_vector2Vec.push_back(Vector2(inf, inf));
	s_vector3Vec.push_back(Vector3(inf, inf, inf));
	s_vector4Vec.push_back(Vector4(inf, inf, inf, inf));
	s_matrixVec.push_back(Matrix(Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf), 
								 Vector4(inf, inf, inf, inf)));

	// TODO: nans, q-nans & denormalized...
}

class Vector2AddTest : public TestCase
{
public:
	Vector2AddTest() : TestCase("Add") {}
	~Vector2AddTest() {}

	bool Test() const;
};

bool FuzzyEqual(float rhs, float lhs)
{
	const float epsilon = 0.0001f;
	return fabs(rhs - lhs) <= epsilon;
}

bool ftest(float rhs, float lhs)
{
	// For testing purposes, two Nans are "equal" even if they have different bits.
	if (isnan(rhs) && isnan(lhs))
		return true;

	return rhs == lhs;
}

bool Vector2AddTest::Test() const
{
	for (unsigned int i = 0; i < s_vector2Vec.size(); ++i)
	{
		for (unsigned int j = 0; j < s_vector2Vec.size(); ++j)
		{
			Vector2 a = s_vector2Vec[i];
			Vector2 b = s_vector2Vec[j];

			float ax = a.X();
			float ay = a.Y();
			float bx = b.X();
			float by = b.Y();
			float rx = ax + bx;
			float ry = ay + by;

			Vector2 r = a + b;

			if (!ftest(r.X(), rx) || !ftest(r.Y(), ry))
			{
				printf("a = (%.5f, %.5f), b = (%.5f, %.5f)\n", a.X(), a.Y(), b.X(), b.Y());
				printf("r = (%.5f, %.5f), rx = %.5f, ry = %.5f\n", r.X(), r.Y(), rx, ry);
				return false;
			}
		}
	}
	return true;
}

int main(int argc, char* argv[])
{
	InitData();

	TestSuite vector2Suite("Vector2");
	vector2Suite.AddTest(new Vector2AddTest());
	vector2Suite.RunTests();

	return 0;
}
