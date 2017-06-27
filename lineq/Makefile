CFLAGS += -I../matrix -I../vector -Wall -O
LDLIBS += -lm 

all : outputA.txt outputB.txt outputC.txt
	cat ./$^

outputA.txt : mainA 
	./$< > $@

outputB.txt : mainB 
	./$< > $@

outputC.txt : mainC 
	./$< > $@

mainC : mainC.o lineq.o matrix.o vector.o

mainB : mainB.o lineq.o matrix.o vector.o

mainA : mainA.o lineq.o matrix.o vector.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o mainA mainB mainC *.txt
