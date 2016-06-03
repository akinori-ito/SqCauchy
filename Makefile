CCC=g++
CCFLAGS=-O3 -Wall
#CCFLAGS=-g -Wall

.cc.o:
	$(CCC) $(CCFLAGS) -c $<

mmrecog: main.o mixture.o addlog.o
	$(CCC) $(CCFLAGS) -o mmrecog main.o mixture.o addlog.o -lm

main.o: main.cc mixture.h distribution.h
mixture.o: mixture.cc mixture.h addlog.h
addlog.o: addlog.cc addlog.h
