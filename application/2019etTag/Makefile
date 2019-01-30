CC = g++ -Wall -Wextra -O3 -march=native -std=c++17

all:	2antOhneXy 3antMitXy

2antOhneXy:	2antOhneXy.cc
	$(CC) 2antOhneXy.cc -o2antOhneXy
	strip 2antOhneXy

3antMitXy:	3antMitXy.cc ../../eiBauBeDi.o
	$(CC) 3antMitXy.cc -o3antMitXy ../../eiBauBeDi.o
	strip 3antMitXy
