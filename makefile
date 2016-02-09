# Makefile for REGCAP
# 2/2/16 LIR
CC=g++

OBJECTS=main.o functions.o
EXE=rc

regcap: $(OBJECTS) functions.h
	$(CC) $(OBJECTS) -o $(EXE)

main.o: main.cpp functions.h
	$(CC) -c main.cpp

functions.o: functions.cpp functions.h
	$(CC) -c functions.cpp

clean:
	rm $(OBJECTS) $(EXE)

