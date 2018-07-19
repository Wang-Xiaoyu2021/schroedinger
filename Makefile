FC=gfortran

qm1d.exe: qm1d.o
	$(FC) $^ -o $@

%.o: %.f
	$(FC) -c $<

clean:
	rm -rf *.o
