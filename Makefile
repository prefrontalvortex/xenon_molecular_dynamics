CC=gcc
CFLAGS=-Wall -pg
# Library flag is important!
LDFLAGS=-pthread -pg -lm
OFLAGS=-O3
TARGET=mol.run
SOURCES=main.c aux.c

OBJS=$(SOURCES:.c=.o)

all: $(TARGET)

original:
	g++ -O3 original.cpp -o mol_orig.run

fast: $(OBJS)
	$(CC) $(OFLAGS) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $(TARGET)

fastc:
	gcc main.c aux.c $(OFLAGS) $(CFLAGS) $(LDFLAGS) -o mol_c.run

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(CFLAGS) $(OBJS) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $(TARGET)

clean:
	rm -f *.o core $(TARGET)
