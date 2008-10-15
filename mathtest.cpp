#include "math3d.h"
#include <stdio.h>
#include <stdlib.h>

float u[3];
float v[3];
Vector3 x;
Vector3 y;
Vector3 r;

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

	r = x + y;

	printf("r = %.5f, %.5f, %.5f\n", r.x, r.y, r.z);

	return 0;
}
