#include "math3d.h"
#include <stdio.h>
#include <stdlib.h>
#include <CoreServices/CoreServices.h>
#include <mach/mach.h>
#include <mach/mach_time.h>

float u[3];
float v[3];
Vector3 x;
Vector3 y;
Vector3 r;

// returns microseconds
float TimeDiff(uint64_t start, uint64_t end)
{
    uint64_t absTime = end - start;
    Nanoseconds nanosec = AbsoluteToNanoseconds( *(AbsoluteTime*)&absTime);
    return (float) UnsignedWideToUInt64(nanosec) / 1000.0;
}

int main(int argc, char* argv[])
{
	u[0] = atof(argv[1]);
	u[1] = atof(argv[2]);
	u[2] = atof(argv[3]);

	v[0] = atof(argv[4]);
	v[1] = atof(argv[5]);
	v[2] = atof(argv[6]);

	x.Set(u[0], u[1], u[2]);
	y.Set(v[0], v[1], v[2]);

	uint64_t start = mach_absolute_time();
	
	r = x;
	for (int i = 0; i < 100000; ++i)
	{
		r = y + x + x + y + r;
	}

	uint64_t end = mach_absolute_time();

	printf("r = %.5f, %.5f, %.5f\n", r.X(), r.Y(), r.Z());

	printf("took %.5f usec\n", TimeDiff(start, end));

	return 0;
}
