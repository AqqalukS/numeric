CFLAGS += -I../matrix -I../vector -Wall -O
LDLIBS += -lm #-lgc `gsl-config --libs`

outA.txt log : mainA
	./$< 1> $@ 2> log

mainA : mainA.o ode.o vector.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o mainA *.txt
