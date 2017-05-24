#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, const char *argv[])
{
	for (int i = atof(argv[1]); i <= atof(argv[3]); i++) {
		printf("%g \n", ceil(pow(10,atof(argv[2])*i)));
	}
	return 0;
}
