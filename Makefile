FLAGS=-O3 -Wall
RM=rm -f
CC=gcc

clean:
	$(RM) xhtcg xopenmp xopenacc

serial:
	$(CC) $(FLAGS) htcg.c -o xhtcg -lm -fopenmp

openmp:
	$(CC) $(FLAGS) openmp.c -o xopenmp -lm -fopenmp

openacc:
	nvc -fast -Minfo=acc -acc openacc.c -o xopenacc

xserial:
	@echo "EXEC=xhtcg" > make.sh

xopenmp:
	@echo "EXEC=xopenmp" > make.sh

xopenacc:
	@echo "EXEC=xopenacc" > make.sh

big:
	@echo "INPUTFILE=input_big.dat" >> make.sh

small:
	@echo "INPUTFILE=input_small.dat" >> make.sh

input:
	@echo "INPUTFILE=$(FILE)" >> make.sh

run:
	@. ./make.sh; ./$$EXEC < ./inputs/$$INPUTFILE > ./outputs/out_$$EXEC.txt
