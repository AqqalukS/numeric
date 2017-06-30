CFLAGS += -I../vector -I../lineq -Wall -O
LDLIBS += -lm 
PLOTTER = gnuplot
N = 20 #number of points

all : plotA.png plotB.png plotC.png

plotA.png  plotB.png plotC.png : plot.gpi outA.txt dataA.txt outB.txt dataB.txt outC.txt dataC.txt
	$(PLOTTER) $<

outA.txt dataA.txt : mainA Makefile
	./$< $(N) 1> $@ 2> dataA.txt

outB.txt dataB.txt : mainB Makefile
	./$< $(N) 1> $@ 2> dataB.txt

outC.txt dataC.txt : mainC Makefile
	./$< $(N) 1> $@ 2> dataC.txt

mainA : mainA.o interp.o vector.o 

mainB : mainB.o interp.o vector.o 

mainC : mainC.o interp.o vector.o 

mainA.o mainB.o mainC.o interp.o : interp.h

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o *.png mainA mainB mainC *.txt
