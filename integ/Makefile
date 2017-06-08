CFLAGS += -Wall -O `gsl-config --cflags`
LDLIBS += `gsl-config --libs` 
LDLIBS += -lm 

all : data
	cat ./$<
data log : main
	./$< 1> $@ 2> log

main : main.o integ.o integ_gsl.o 

clean :
	$(RM) *.o main data log *.png
