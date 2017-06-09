CFLAGS += -I../matrix -I../vector -Wall -O
LDLIBS += -lm -lgc `gsl-config --libs`
#CC = gcc --std=gnu99
NSTART = 15 # 10^(NSTEP*(NSTART:NEND)
NSTEP = 0.1
NEND = 25

.PHONEY : jacobi

all : plot.png 

plot.png : plot.gpi times.txt
	gnuplot ./$<

times.txt rotations.txt : main_cyclic main_by_row main_classic plot_seq.txt Makefile
	cat /dev/null > $@
	cat /dev/null > rotations.txt
	for N in `cat $(word 4, $^)`; do \
	echo "N=$$N"; time -f "$$N %U " -ao $@ ./$< $$N >> rotations.txt; \
	done
	echo " " >> $@
	echo " " >> $@
	echo " " >> rotations.txt
	echo " " >> rotations.txt
	for N in `cat $(word 4, $^)`; do \
	echo "N=$$N"; time -f "$$N %U " -ao $@ ./$(word 2, $^) $$N >> rotations.txt; \
	done
	echo " " >> $@
	echo " " >> $@
	echo " " >> rotations.txt
	echo " " >> rotations.txt
	for N in `cat $(word 4, $^)`; do \
	echo "N=$$N"; time -f "$$N %U " -ao $@ ./$(word 3, $^) $$N >> rotations.txt; \
	#echo "N=$$N"; \time -ao $@ -f "$$N %U" ./$(word 3, $^) $$N >> rotations.txt; \
	done

plot_seq.txt : plot_seq Makefile
	./$< $(NSTART) $(NSTEP) $(NEND) > $@

main_cyclic : main_cyclic.o  jacobi.o matrix.o vector.o 

main_by_row : main_by_row.o  jacobi.o matrix.o vector.o 

main_classic : main_classic.o  jacobi.o matrix.o vector.o 

jacobi : jacobi_test.txt
	cat ./$<

jacobi_test.txt : jacobi_test 
	./$< 5 > $@

jacobi_test : jacobi_test.o jacobi.o matrix.o vector.o

matrix.o : ../matrix/matrix.c ../matrix/matrix.h
	$(CC) -c $(CFLAGS) $< -o $@

vector.o : ../vector/vector.c ../vector/vector.h
	$(CC) -c $(CFLAGS) $< -o $@

clean :
	$(RM) *.o jacobi_test main_cyclic main_by_row main_classic plot_seq *.txt *.png
