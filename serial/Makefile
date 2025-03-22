FLAGS=-O3 -Wall
RM=rm -f
CC=gcc

EXEC=acc

all: $(EXEC)

$(EXEC):
	$(CC) $(FLAGS) $(EXEC).c -o x$(EXEC) -lm -fopenmp -pg

clean:
	$(RM) $(EXEC).o $(EXEC)
