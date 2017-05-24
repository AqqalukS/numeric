CFLAGS += -I../matrix -I../vector -Wall -O
LDLIBS += -lm -lgc `gsl-config --libs`
CC = gcc --std=gnu99
NSTART = 10 # 10^(NSTEP*(NSTART:NEND)
NSTEP = 0.1
NEND = 30

.PHONEY : jacobi

all : plot.times.png plot.rotations.png

plot.rotations.png : plot.rotations.gpi rotations2.txt
	gnuplot ./$<

plot.times.png : plot.times.gpi times2.txt
	gnuplot ./$<

#times.txt rotations.txt : main_cyclic main_by_row plot_seq.txt Makefile
########cat /dev/null > $@
########cat /dev/null > rotations.txt
########for N in `cat $(word 3, $^)`; do \
########echo "N=$$N"; \gtime -ao $@ -f "$$N %U" ./$< $$N >> rotations.txt; \
########done
########echo " " >> $@
########echo " " >> $@
########echo " " >> rotations.txt$
########for N in `cat $(word 3, $^)`; do \
########echo "N=$$N"; \gtime -ao $@ -f "$$N %U" ./$(word 2, $^) $$N >> rotations.txt; \
########done

plot_seq.txt : plot_seq Makefile
	./$< $(NSTART) $(NSTEP) $(NEND) > $@

main_cyclic : main_cyclic.o  jacobi.o matrix.o vector.o 

main_by_row : main_by_row.o  jacobi.o matrix.o vector.o 

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
	$(RM) *.o jacobi_test main_cyclic main_by_row plot_seq 
