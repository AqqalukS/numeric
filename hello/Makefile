.PHONEY : A

A : out.A.txt
	cat $<

out.A.txt : main-A
	./$< > $@

.PHONEY : clean

clean : 
	$(RM) main-A *.txt
