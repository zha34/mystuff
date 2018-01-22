all: rowColBlock

CC = gcc
CCFLAGS = -Wall  -O3
LDFLAGS = -Wall  -O3
LFLAGS = -lm -lpthread

rowColBlock: rowColBlock.o
	$(CC) $(LDFLAGS) -o rowColBlock rowColBlock.o  $(LFLAGS)

rowColBlock.o: rowColBlock.c
	$(CC) $(CCFLAGS) -c rowColBlock.c -o rowColBlock.o

clean:
	-rm *.o *~ rowColBlock
