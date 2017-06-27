CFLAGS += -I../vector -I../lineq -Wall -O
LDLIBS += -lm 
PLOTTER = gnuplot
N = 20 #number of points

plot.png : plot.gpi interpolation.txt datapoints.txt
	$(PLOTTER) $<

interpolation.txt datapoints.txt : main Makefile
	./$< $(N) 1> $@ 2> datapoints.txt

main : main.o interp.o vector.o 

main.o interp.o : interp.h

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o *.png main *.txt
