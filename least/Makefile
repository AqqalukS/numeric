CFLAGS += -I../matrix -I../vector -I../lineq -I../eigen -Wall -O
LDLIBS += -lm 

all : plot.png

plot.png : plot.gpi data_AB.txt data_C.txt 
	gnuplot ./$<

data_AB.txt : mainAB
	./$< > $@

data_C.txt  : mainC
	./$< > $@ 

mainAB : funs.o least_square.o matrix.o vector.o lineq.o jacobi.o

mainC : funs.o least_singular.o  matrix.o vector.o lineq.o jacobi.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@
 
vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

lineq.o : ../lineq/lineq.c ../lineq/lineq.h
	$(CC) -c $(CFLAGS) $< -o $@

jacobi.o : ../eigen/jacobi.c ../eigen/jacobi.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o mainAB mainC  *.txt *.png
