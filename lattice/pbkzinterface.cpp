#ifndef _inc_pbkzinterface_cpp
#define _inc_pbkzinterface_cpp

#include "pbkzproperty.cpp"
#include "pruningfunc.cpp"


namespace progressive_bkz {
    
    void fixedparamBKZ(mat_ZZ& L,std::string opts,BKZproperty& BP) {
        
        extractoptions(BP,opts);
        //
        //Todo: detecting ENUM speeds
        //
        sharememalloc();
        
        //reading pruning coefficients file
        pruning_func::init_pruning_func();
        
        BKZmain(L,0,BP);

    }
    
    void fixedparamBKZ(mat_ZZ& L,std::string opts) {
        BKZproperty BP;
        fixedparamBKZ(L,opts,BP);
    }


    int fixedparamBKZ(mat_ZZ &L,int index,int beta,double prob, double alpha, int tourlim,int vl,int opt) {

        sharememalloc();
        pruning_func::init_pruning_func();

        BKZproperty BP;
        BP.enumlim[0] = 1000000;    // limit of #processed nodes (/10^8) 
        BP.enumlim[1] = 1; 

        BP.beta[0] = beta;  //blocksize
        BP.beta[1] = beta*0.6;    //blocksize in preprocess strategy1
        BP.tourlim = tourlim;
        BP.breakindex = -1;


        //multithread?
        if (opt & OPT_MULTITHREAD) {
            BP.multithread = true;
            BP.numthreads = omp_get_max_threads();
            BP.MTlimit = BP.numthreads*10000000;
            cout << "setting # threads = " << BP.numthreads << endl;
        } else {
            BP.multithread = false;
            BP.MTlimit = 10000000;
        }

        //optimize pruning function
        if (opt & OPT_OPTIMIZE_PRUNING_FUNCTION) {
            BP.optimizepf = true;
            BP.optimizepf_at_least = 1000;
            lattice_enum::enum_speed_bench(BP.numthreads);
        }

        //do preprocess?
        if (opt & (OPT_PREPROCESS | OPT_PREPROCESS2)) {
            BP.preprocess = true;
            BP.preprocess_at_least = 1;  
            BP.preprocess_strategy = 1;
            if (opt & OPT_PREPROCESS2) BP.preprocess_strategy = 2;
            lattice_enum::enum_speed_bench(BP.numthreads);
        } else {
            BP.preprocess = false;
        }

        if (opt & OPT_EXTENDBLOCKSIZE) {
            BP.extend_blocksize=true;
            BP.extend_blocksizemult = 2;   //extend blocksize while expected cost <= max cost of the current tour
        }


        //break at a specific index
        if (opt & OPT_FIRSTINDEX) {
            BP.breakindex = 1;
        }

        //break at a specific index
        if (opt & OPT_GHBASEDSKIP) {
            BP.ghskip = true;
            BP.ghskipcoeff = 1.025;  //if |b*i| < a*GH(L), skip the block
        }

        //output log file
        if (opt & OPT_TIMELOG) {
             BP.tlname="bkzlog.txt"; 
        }

        BP.verboselevel = vl;
        BP.pruning_prob = prob;

        BP.init_radius = alpha;
        BP.init_mode = 'G';  //R=alpha*GH(L)

        BP.holdvecs = 16;

        //A heuristic strategy for finding short vectors 
        if (opt & OPT_FIND_SHORT) {
            BP.process_entire_basis = true;
            BP.ec_boundalpha = 1.05;
        }    

        cputime = 0;    //CPU time in second
        start = clock();

        return BKZmain(L,0,BP);
    }

    void progressiveBKZ(mat_ZZ &L,int startbeta,int endbeta) {
        int beta;
        double r,abar,alpha,prob;
        std::stringstream opts;

        int n = L.NumRows();
        quad_float* c;
        c = new quad_float[n+1];
        double logfec;
        
        for (beta=startbeta;beta<=min(endbeta,n);beta++) {
            r= compute_target_r(beta);
            abar = 1.0 / lattice_tools::gaussian_heuristic_constant[beta] * exp(-0.25*(beta-1.0)*log(r));
            alpha = abar * ((1.0+beta)/beta);
            if (alpha<1.00) alpha = 1.0;
            prob = 2.0 * exp(-1.0*beta*log(alpha));

            BKZGSsim(c,n,beta-1);  //outputs GS length with det=1    
            logfec = to_double(log(lattice_tools::FullENUMCost(c,n,INPUT_NONSQUARED)));
            //logfec = ComputeLogTargetFEC(n,beta-1);
            
            opts.str("");
            if (n>50) {
                quad_float pred_maxcost = 2.0 * 0.25 * beta * ENUMCost(c,n-beta+1,beta,prob,alpha,'G',0,INPUT_NONSQUARED,1);
                opts << " maxcost=" << pred_maxcost;
            }
            
            opts << " bs=" << beta << " prob=" << prob << " alpha=" << alpha << " maxtour=500";
            opts << " minlogfec=" << logfec << " optpfunc preprocess extendblocksize";
            
            opts << " heuristicbreak";
             
            if (beta >= 70) {
                opts << " verbose=3"; 
            } else
            {
                opts << " verbose=1"; 
            }
            
            cout << "options=" << opts.str() << endl;
            fixedparamBKZ(L,opts.str());
        }
        
    }    


    
    
    
}   // end of namespace progressive bkz






#endif