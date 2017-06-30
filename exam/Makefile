CFLAGS += -I../matrix -I../vector -Wall -O
LDLIBS += -lm #-lgc `gsl-config --libs`
N = 5

all : outA.txt outB.txt
	cat ./$^

outA.txt log : mainA
	./$< $(N) 1> $@ 2> log

outB.txt log : mainB
	./$< $(N) 1> $@ 2>> log

mainA : mainA.o lanczos.o gs_ortho.o  vector.o matrix.o

mainB : mainB.o arnoldi.o lanczos.o gs_ortho.o vector.o matrix.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o mainA mainB log out* 
