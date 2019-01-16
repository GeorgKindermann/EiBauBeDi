CC = g++ -Wall -Wextra -O3 -march=native -std=c++17

all:	eiBauBeDi.o testEiBauBeDi

eiBauBeDi.o:	eiBauBeDi.cc eiBauBeDi.h
	$(CC) eiBauBeDi.cc -c

testEiBauBeDi:	testEiBauBeDi.cc eiBauBeDi.o
	$(CC) testEiBauBeDi.cc -otestEiBauBeDi eiBauBeDi.o
	strip testEiBauBeDi
