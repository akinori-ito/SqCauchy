CCC=g++
CCFLAGS=-O3 -Wall
#CCFLAGS=-g -Wall

.cc.o:
	$(CCC) $(CCFLAGS) -c $<

sqcupdate: SQCupdate.o SQC_train.o
	$(CCC) $(CCFLAGS) -o sqcupdate SQCupdate.o SQC_train.o -lm

mmrecog: main.o mixture.o addlog.o
	$(CCC) $(CCFLAGS) -o mmrecog main.o mixture.o addlog.o -lm

main.o: main.cc mixture.h distribution.h
mixture.o: mixture.cc mixture.h addlog.h
addlog.o: addlog.cc addlog.h

SQCupdate.o: SQCupdate.cc SQC_train.h Adam.h
SQC_train.o: SQC_train.cc SQC_train.h
