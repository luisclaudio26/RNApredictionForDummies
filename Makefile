CC=gcc
CFLAGS=-std=c99 -g -DVERBOSE

main: main.o nussanov.o
	$(CC) $(CFLAGS) -o RNApredict main.o nussanov.o

main.o: main.c
	$(CC) $(CFLAGS) -c main.c -o main.o

nussanov.o: nussanov.c
	$(CC) $(CFLAGS) -c nussanov.c -o nussanov.o

clean:
	rm *.o
