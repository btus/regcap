# Makefile for REGCAP
# 2/2/16 LIR
CC=g++

OBJECTS=main.o functions.o config.o log.o
EXE=rc

regcap: $(OBJECTS) functions.h config/config.h
	$(CC) $(OBJECTS) -o $(EXE)

main.o: main.cpp functions.h
	$(CC) -c main.cpp

functions.o: functions.cpp functions.h
	$(CC) -c functions.cpp

config.o: config/config.cpp config/config.h
	$(CC) -c config/config.cpp

log.o: config/log.cpp config/log.h
	$(CC) -c config/log.cpp

clean:
	rm $(OBJECTS) $(EXE)

