CFLAGS += -I../matrix -I../vector -I../lineq -Wall -O
LDLIBS += -lm 
CC = gcc-4.9 -std=c11

all : print plot.png

plot.png : plot.gpi log mainB
	gnuplot ./$<

print : data
	cat ./$<

data log : mainA mainB
	./$< 1> $@ 2> log 
	./$(word 2, $^) 1>> $@ 2>> log 

mainB : mainB.o funcs.o optim.o matrix.o vector.o lineq.o

mainA : mainA.o funcs.o optim.o matrix.o vector.o lineq.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@
 
vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

lineq.o : ../lineq/lineq.c ../lineq/lineq.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o main data log *.png
