#include "math3d.h"
#include <stdio.h>
#include <stdlib.h>
#include "unittest.h"

#include <CoreServices/CoreServices.h>
#include <mach/mach.h>
#include <mach/mach_time.h>

static const int RANDOM_DATA_SIZE = 1024 * 1024;

Vector3 s_randomVector3[RANDOM_DATA_SIZE];
Vector4 s_randomVector4[RANDOM_DATA_SIZE];
Matrix s_randomMatrix[RANDOM_DATA_SIZE];

float s_floatResult;
Vector4 s_vector4Result;

// 0..1
float RandomFloat()
{
	return (float)rand() / RAND_MAX;
}

void InitRandomData()
{
	srand((unsigned int)mach_absolute_time());
	for (int i = 0; i < RANDOM_DATA_SIZE; ++i)
	{
		s_randomVector3[i].Set(RandomFloat(), RandomFloat(), RandomFloat());
		s_randomVector4[i].Set(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());

		Vector4 row0(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
		Vector4 row1(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
		Vector4 row2(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
		Vector4 row3(RandomFloat(), RandomFloat(), RandomFloat(), RandomFloat());
		s_randomMatrix[i] = Matrix(row0, row1, row2, row3);
	}
}

// returns microseconds
float TimeDiff(uint64_t start, uint64_t end)
{
    uint64_t absTime = end - start;
    Nanoseconds nanosec = AbsoluteToNanoseconds( *(AbsoluteTime*)&absTime);
    return (float) UnsignedWideToUInt64(nanosec) / 1000.0;
}

class BlockTimer
{
public:
	BlockTimer(const char* descIn) 
	{
		desc = descIn; 
		start = mach_absolute_time();
	}

	~BlockTimer()
	{
		uint64_t end = mach_absolute_time();
		printf("%s took %.1f usec\n", desc, TimeDiff(start, end));
	}

	uint64_t start;
	const char* desc;
};

void TimeVector4Add()
{
	Vector4 r;
	{
		BlockTimer timer("Vector4 addition");

		r = s_randomVector4[0];
		for (int i = 0; i < RANDOM_DATA_SIZE - 3; ++i)
		{
			r = s_randomVector4[i] + s_randomVector4[i+1] + s_randomVector4[i+2] + s_randomVector4[i+3] + r;
		}
	}
	s_vector4Result = r;
}

void TimeVector4DotProduct()
{
	float a = 1.0f;
	{
		BlockTimer timer("Vector4 dot");

		for (int i = 0; i < RANDOM_DATA_SIZE - 1; ++i)
		{
			a += s_randomVector4[i] * s_randomVector4[i+1];
		}
	}
	s_floatResult = a;
}

void TimeMul4x4()
{
	Vector4 a;
	{
		BlockTimer timer("Mul4x4");
		
		for (int i = 0; i < RANDOM_DATA_SIZE; ++i)
		{
			a = a + Mul4x4(s_randomMatrix[i], s_randomVector4[i]);
		}
	}
	s_vector4Result = a;
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

bool Vector2AddTest::Test() const
{
	Vector2 a(1, 1);
	Vector2 b(2, 2);
	Vector2 c = a + b;

	if (FuzzyEqual(c.X(), 3.0f) && FuzzyEqual(c.Y(), 3.0f))
		return true;
	else
		return false;
}

int main(int argc, char* argv[])
{
	InitRandomData();

	TimeVector4Add();
	TimeVector4DotProduct();
	TimeMul4x4();

	TestSuite vector2Suite("Vector2");
	vector2Suite.AddTest(new Vector2AddTest());
	vector2Suite.RunTests();

	return 0;
}
