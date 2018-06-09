CC	= g++
SRC	= main.cpp
BOOST_DIR = ./boost_1_66_0
CFLAGS	= -std=c++11 -march=native -Ofast -funroll-loops -flto -fomit-frame-pointer -fprefetch-loop-arrays -msse4 -fpermissive
LDFLAGS	= -lgmp -lmpfr -lgsl -lntl -fopenmp
INCL	= -I. -I${BOOST_DIR}

all:
	${CC} ${CFLAGS} ${SRC} ${INCL} ${LDFLAGS}


