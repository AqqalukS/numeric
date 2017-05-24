CFLAGS += -I../matrix -I../vector -Wall -O
LDLIBS += -lm 
CC = gcc --std=gnu99

data : main 
	./$< > $@

main : main.o lineq.o matrix.o vector.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o main
