CC = g++ -Wall -O3 -march=native -mtune=native -pipe

all:	eiBauBeDi.o testEiBauBeDi

eiBauBeDi.o:	eiBauBeDi.cc eiBauBeDi.h
	$(CC) eiBauBeDi.cc -c

testEiBauBeDi:	testEiBauBeDi.cc eiBauBeDi.o
	$(CC) testEiBauBeDi.cc -otestEiBauBeDi eiBauBeDi.o
	strip testEiBauBeDi
