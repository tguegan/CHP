all: cleanall
	mpicc main.c init.c solver.c schwarz.c -o main -lm -O3 -Wno-unused-result
cleanall:
	rm -rf main *~ *# *.plt
clean:
	rm -rf *.plt *~ *#
