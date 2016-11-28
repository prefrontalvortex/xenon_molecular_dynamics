CC=g++
CFLAGS=-Wall -pg
# Library flag is important!
LDFLAGS=-pthread -pg
OFLAGS=-O3
TARGET=mol
SOURCES=main.cpp

OBJS=$(SOURCES:.c=.o)

all: $(TARGET)

fast: $(OBJS)
	$(CC) -o $(TARGET) $(OFLAGS) $(CFLAGS) $(OBJS) $(LDFLAGS)


$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(CFLAGS) $(OBJS) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $(TARGET)

clean:
	rm -f *.o core $(TARGET)