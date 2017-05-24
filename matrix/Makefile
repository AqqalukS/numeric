CFLAGS += -I../vector/ -Wall -O

all : data
	cat $<

data log : matrix_test
	./$< 1> data 2> log

matrix_test : matrix_test.o matrix.o vector.o

matrix.o : matrix.h 

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o matrix_test data 
