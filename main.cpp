//comment out when publishing
#define __develop
//#define __ide

int debug_output=0;
//#define debug_display_info

#include <lattice/pbkz.hpp>
#include "lattice/bkztest.cpp"


int main(int argc, char** argv) {

    bkztest();  exit(0);

}

//working compiling command: 
//g++ main.cpp -I. -I../../boost_1_64_0 -fopenmp -lntl -lgsl -lgmp -lmpfr -fpermissive -std=c++11 
//Faster options
//g++ main.cpp -I. -I../../boost_1_64_0 -fopenmp -lntl -lgsl -lgmp -lmpfr -fpermissive -std=c++11 -march=native -msse4 -msse4.1 -mavx --param max-inline-insns-single=5000 --fast-math -flto -fomit-frame-pointer -fprefetch-loop-arrays

//g++ main.cpp -I. -I../../boost_1_64_0 -lntl -lgmp -fopenmp -lgsl -std=c++11 -lntl -O6 -fopenmp -funroll-loops -march=native --param max-inline-insns-single=5000 -ffast-math -flto -fomit-frame-pointer -fprefetch-loop-arrays -msse4 -lmpfr

