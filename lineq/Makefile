CFLAGS += -I../matrix -I../vector -Wall -O
LDLIBS += -lm 
CC = gcc-4.9 -std=c11

all : data
	cat ./$<

data : mainA mainB mainC
	./$< > $@
	./$(word 2, $^) >> $@
	./$(word 3, $^) >> $@

mainC : mainC.o lineq.o matrix.o vector.o

mainB : mainB.o lineq.o matrix.o vector.o

mainA : mainA.o lineq.o matrix.o vector.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o mainA mainB mainC data
