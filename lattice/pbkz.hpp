#ifndef _inc_pbkz_hpp
#define _inc_pbkz_hpp

//configures
#define _include_svpchallenge_generator
#define _include_idealchallenge_generator

//system includes
#include <cstdlib>
#include <sys/times.h>
#include <sys/stat.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>     
#include <algorithm>
#include <functional>        
using namespace std;

//NTL
#include <NTL/ZZ.h>
#include <NTL/matrix.h>
#include <NTL/vec_ZZ.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_RR.h>
#include <NTL/LLL.h>
NTL_CLIENT

#include <omp.h>

//BKZ options (old interface)
#define OPT_PREPROCESS 0x01
#define OPT_MULTITHREAD 0x02
#define OPT_FIRSTINDEX 0x04
#define OPT_PREPROCESS2 0x08
#define OPT_TIMELOG 0x10
#define OPT_EXTENDBLOCKSIZE 0x20
#define OPT_GHBASEDSKIP 0x40
#define OPT_FIND_SHORT 0x80
#define OPT_OPTIMIZE_PRUNING_FUNCTION 0x100

        
//To define verbose level (for readability in function call)
#define VL0 0
#define VL1 1
#define VL2 2
#define VL3 3
#define VL4 4
        
#define optlevel0 0
#define optlevel1 1
#define optlevel2 2
        
//Define the maximum dimention of lattices
//(used to allocate the shared memories)
#define latticemaxdim 1000

//Max. # vectors returned by the enumeration subroutine
#define maxfoundvec 1000
        
        
//Fundemental tools
#include "timetools.cpp"        
#include "filetools.cpp"        
#include "latticetools.cpp"        
#include "vectorset.cpp"        
        
//LLL subroutines are used
#include "../external/LLL_RR.c"

//The following defines are adopted to NTL 9.7.0
#define RowTransform RowTransformQP
#define RowTransform2 RowTransform2QP
#define red_fudge red_fudgeQP
#define init_red_fudge init_red_fudgeQP
#define inc_red_fudge inc_red_fudgeQP
#define log_red log_redQP
#define verbose verboseQP
#define NumSwaps NumSwapsQP
#define StartTime StartTimeQP
#define LastTime LastTimeQP
#define LLLStatus LLLStatusQP
#define BKZConstant BKZConstantQP
#define BKZThresh BKZThreshQP
#define BKZStatus BKZStatusQP
#define ComputeBKZConstant ComputeBKZConstantQP
#include "../external/LLL_QP.c"
/*
#undef RowTransform
#undef RowTransform2
#undef red_fudge
#undef log_red
#undef verbose
   */      
//Tools for lattice vector enumeration
#include "pruningfunc.cpp"
#include "vectorenumeration.cpp"
#include "vectorenumeration_close.cpp"

        //progressive BKZ routines
#include "pbkzsharemem.cpp"
#include "pbkzproperty.cpp"
#include "pbkzsupport.cpp"
#include "pbkzmain.cpp"
#include "pbkzinterface.cpp"


#ifdef _include_svpchallenge_generator
    #include "../external/tools.h"
#endif

#ifdef _include_idealchallenge_generator
    #include <NTL/ZZX.h>
    #include "../external/ideal.h"
#endif
        
        
#include <string>
#include <vector>
#include <map>

#include <boost/algorithm/string.hpp>

typedef std::map<std::string,std::string> stringmap;
typedef std::vector<std::string> bkzstrategy;
typedef std::vector<double> pruningfunction;





#endif
