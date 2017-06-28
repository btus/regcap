# Makefile for REGCAP
# 3/16/16 LIR
CC=g++

OBJECTS=main.o functions.o config.o log.o weather.o psychro.o equip.o gauss.o moisture.o
EXE=rc

regcap: $(OBJECTS) functions.h config/config.h
	$(CC) $(OBJECTS) -o $(EXE)

main.o: main.cpp functions.h weather.h 
	$(CC) -c main.cpp

functions.o: functions.cpp functions.h constants.h gauss.h
	$(CC) -c functions.cpp

gauss.o: gauss.cpp gauss.h
	$(CC) -c gauss.cpp

weather.o: weather.cpp weather.h constants.h
	$(CC) -c weather.cpp

equip.o: equip.cpp equip.h constants.h psychro.h
	$(CC) -c equip.cpp

moisture.o: moisture.cpp moisture.h constants.h psychro.h gauss.h
	$(CC) -c moisture.cpp

psychro.o: psychro.cpp psychro.h constants.h
	$(CC) -c psychro.cpp

config.o: config/config.cpp config/config.h
	$(CC) -c config/config.cpp

log.o: config/log.cpp config/log.h
	$(CC) -c config/log.cpp

clean:
	rm $(OBJECTS) $(EXE)

