all: clean
	mpicc.openmpi main.c init.c solver.c schwarz.c -o main -lm
clean:
	rm -rf *~ main *# *.plt
