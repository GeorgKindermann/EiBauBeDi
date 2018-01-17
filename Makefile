CC = g++ -Wall -O3 -funroll-loops -march=native -mtune=native

all:	eiBauBeDi.o testEiBauBeDi

eiBauBeDi.o:	eiBauBeDi.cc eiBauBeDi.h
	$(CC) eiBauBeDi.cc -c

testEiBauBeDi:	testEiBauBeDi.cc eiBauBeDi.o
	$(CC) testEiBauBeDi.cc -otestEiBauBeDi eiBauBeDi.o
	strip testEiBauBeDi
