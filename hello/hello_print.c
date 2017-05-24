#include <stdio.h>
#include <stdlib.h>


void print_hello () {
	printf("Hello, ");
}

void print_me () {
	char *user = getenv("USER");
	printf("%s", user);
}
