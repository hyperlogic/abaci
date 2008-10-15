#include "math3d.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[])
{
	float u[3];
	u[0] = atof(argv[1]);
	u[1] = atof(argv[2]);
	u[2] = atof(argv[3]);

	float v[3];
	v[0] = atof(argv[4]);
	v[1] = atof(argv[5]);
	v[2] = atof(argv[6]);

	Vector3 x(u[0], u[1], u[2]);
	Vector3 y(v[0], v[1], v[2]);

	Vector3 r = x + y;

	printf("r = %.5f, %.5f, %.5f\n", r.x, r.y, r.z);

	return 0;
}
