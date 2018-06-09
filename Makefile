CC = g++
CFLAG = -std=c++11 -m64 -Ofast -march=native -flto -Wall
CLIB = -lntl -lgmp -fopenmp -lgsl -lgslcblas
INCLUDE = .

all:
	$(CC) $(CFLAG) main.cpp $(CLIB) -I $(INCLUDE)

clean:
	@rm a.out
