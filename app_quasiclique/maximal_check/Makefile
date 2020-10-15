CC = g++
CFLAGS = -O3 -Wall
OBJECTS = main.o graph.o RmNonMax.o data.o
TARGET = quasiCliques

$(TARGET) : $(OBJECTS) Makefile
	$(CC) $(CFLAGS) $(OBJECTS) -o $(TARGET)

main.o: main.cpp graph.h RmNonMax.h
	$(CC) -c $(CFLAGS) main.cpp

graph.o: graph.cpp graph.h  data.h
	$(CC) -c $(CFLAGS) graph.cpp

data.o: data.cpp data.h
	$(CC) -c $(CFLAGS) data.cpp

RmNonMax.o: RmNonMax.cpp RmNonMax.h graph.h
	$(CC) -c $(CFLAGS) RmNonMax.cpp

clean :
	rm *.o

