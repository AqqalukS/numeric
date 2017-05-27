CFLAGS += -I../matrix -I../vector -I../lineq -Wall -O
LDLIBS += -lm 
CC = gcc-4.9 -std=c11
PLOTTER = gnuplot
N = 20 #number of points

plot.png : plot.gpi data log
	$(PLOTTER) $<

data log : main Makefile
	./$<  $(N) 1> data 2> log

main : main.o interp.o vector.o 

main.o interp.o : interp.h

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o *.png main data log
