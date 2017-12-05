CC = g++ -Wall -O3 -funroll-loops -march=native -mtune=native

all:	eiBauBeDi

eiBauBeDi:	eiBauBeDi.cc
	$(CC) eiBauBeDi.cc -oeiBauBeDi
	strip eiBauBeDi
