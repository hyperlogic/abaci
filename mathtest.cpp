#include "math3d.h"
#include "stdio.h"

int main()
{
	Vector3 x(1.0f, 1.0f, 1.0f);
	Vector3 y(1.0f, 1.0f, 1.0f);

	Vector3 r = x + y;

	printf("r = %.5f, %.5f, %.5f\n", r.x, r.y, r.z);

	return 0;
}
