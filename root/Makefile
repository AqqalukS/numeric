CFLAGS += -I../matrix -I../vector -I../lineq -Wall -O
LDLIBS += -lm 

all : dataA.txt dataB.txt dataC.txt
	cat $^

dataA.txt : mainA
	./$< > $@

dataB.txt : mainB
	./$< > $@

dataC.txt : mainC
	./$< > $@

mainA : mainA.o root.o funs.o matrix.o vector.o lineq.o

mainB : mainB.o root.o funs.o matrix.o vector.o lineq.o

mainC : mainC.o root.o funs.o matrix.o vector.o lineq.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@
 
vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

lineq.o : ../lineq/lineq.c ../lineq/lineq.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o mainA mainB mainC *.txt *.png
