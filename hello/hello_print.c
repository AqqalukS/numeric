#include <stdio.h>
#include <stdlib.h>


void print_hello () {
	printf("Hello, ");
}

void print_me () {
	
	printf("%g", getlogin_r ());
}
