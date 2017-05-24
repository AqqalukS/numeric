CFLAGS = -Wall -std=gnu99
PLOTTER = gnuplot

plot.png : plot.gpi data log
	$(PLOTTER) $<

data log : main
	./$< 1> data 2> log

main : main.o interp.o vector.o

main.o interp.o vector.o : interp.h vector.h

clean :
	$(RM) *.o *.png main data log
