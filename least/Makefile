CFLAGS += -I../matrix -I../vector -I../lineq -Wall -O
LDLIBS += -lm 
CC = gcc --std=gnu99

all : plot.png

plot.png : plot.gpi data log
	gnuplot ./$<

data log : main
	./$< 1> $@ 2> log

main : least_square.o matrix.o vector.o lineq.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@
 
vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

lineq.o : ../lineq/lineq.c ../lineq/lineq.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o main data log *.png
