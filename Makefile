FLAGS=-O3 -Wall
RM=rm -f
CC=gcc

EXEC=acc2

all: $(EXEC)

$(EXEC):
	$(CC) $(FLAGS) $(EXEC).c -o x$(EXEC) -lm -fopenmp

clean:
	$(RM) $(EXEC).o $(EXEC)


teste: 
	nvc -fast -Minfo=acc -acc -gpu=cc86 acc2.c -o xacc2
	nvc -fast -Minfo=acc -acc -gp./xacc2 < input > gpu.log
	nvc -fast -Minfo=acc -ta=multicore 