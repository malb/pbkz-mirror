#ifndef _inc_pbkzproperty
#define _inc_pbkzproperty


//memory allocation, working property, and contains some returned values
struct BKZproperty {
    int beta[2];  //blocksize in j-th recursive
    int enumlim[2];       //#limit of enumerationt in j-th recursive
    bool preprocess;
    int verboselevel;
    double pruning_prob;

    //# holding vectors
    int holdvecs;
    
    //initial pruning radius
    double init_radius;
    char init_mode;
    
    //preprocess
    long int preprocess_at_least;
    int preprocess_strategy;
    bool preprocessonly;    //Abort the main routine after the first preprocessing
    int preprocess_max_time;
    
    //optimizing pruning function
    bool optimizepf;
    int optimizepf_at_least;
    
    //multithread process in ENUM subroutine
    char multithread;
    int numthreads;
    long int MTlimit;   // #ucost < MIlimit, then ENUM is single-threaded
    
    //tour setting
    int tourlim;
    int breakindex;
    
    //Aborting condition
    int abortbysec;
    double abortbyr;
    double abortbyfec;
    int abortbyindex;
    double abortbyapprox;
    
    //performance
    double enum_speed;  //Speed in the last enum subroutine

    //log
    std::string tlname; //time-log
    int tllevel;
    
    //skip
    bool ghskip;
    double ghskipcoeff;

    
    //For finding short vector
    bool process_entire_basis;
    double ec_boundalpha;
    
    bool extend_blocksize;
    int extend_blocksizemult;
    double extend_blocksize_initmax;
    bool modifyprobability;
    
    int startindex;
    int endindex;
    int equivrcomputes,equivrcomputee;  //start and end
    int blocktimelimit;
    
    //record found vectors
    bool recvec;
    std::vector<std::vector<int> > recvecdata;
    bool supportVSflag;
    std::vector<std::vector<int> >* supportvs;    

    //Returned values
    quad_float return_tour_max_cost;
    quad_float current_approx_ratio;
    int processed_tours;
    bool record_min_fec;
    mat_ZZ min_fec_basis;
    
    //Heuristic break
    bool heuristicbreak;
    
    //Constructor
    BKZproperty() {
        //default values
        preprocess = false;
        preprocessonly = false;
        preprocess_max_time = -1; 
        verboselevel = 1;
        pruning_prob=0.01;
        
        holdvecs=16;
        init_radius = 1.05;
        init_mode='G';
        
        multithread=false;
        numthreads=1;
        
        tourlim=1;
        breakindex=-1;
        
        extend_blocksize = false;
        extend_blocksizemult = -1;
        extend_blocksize_initmax = 0;

        ghskip = false;
        ghskipcoeff = 1.00;

        blocktimelimit = -1;
        
        enum_speed=-1;

        process_entire_basis = false;
        ec_boundalpha = 1.05;
        
        optimizepf=false;
        optimizepf_at_least = -1;
        modifyprobability=false;
        
        enumlim[0] = 1000000;    // limit of #processed nodes (/10^8) 
        enumlim[1] = 1; 

        beta[0] = 20;
        beta[1] = -1;

        breakindex = -1;
        startindex = 1;
        endindex = -1;
        equivrcomputes = 1;
        heuristicbreak=false;
        
        tlname="";
        tllevel=1;
        
        abortbysec=-1;
        abortbyr = -1;
        abortbyfec = -1;
        abortbyindex = -1;
        abortbyapprox = -1;
         
        processed_tours = 0;
        record_min_fec = false;

        recvec = false;
        supportVSflag = false;
        supportvs=0;
        processed_tours = 0;

        
        
    }
};


namespace progressive_bkz {

    void extractoptions(BKZproperty& BP,std::string opts) {
        //extracting option string
      string::size_type pos;
        std::string com,rcom;

        while (opts.length()!=0) {
            com = "";
            rcom = "";

            pos = opts.find(" ");
            if (pos!=string::npos) {
                com = opts.substr(0,pos);
                opts = opts.substr(pos+1);
            } else {
                com = opts;
                opts = "";
            }
            //cout << "com=" << com << endl;
            pos = com.find("=");
            if (pos!=string::npos) {
                rcom = com.substr(pos+1);
                com = com.substr(0,pos);
                //cout << com << " " << rcom << endl;
            } else {
                //opts = "bs=30 prob=0.01 alpha=1.05 maxtour=10000 verbose=1 timelog=bkzlog100_0.txt extendlog=gsalog100_0.txt preprocess extendblocksize";
            }

            if (com=="bs") BP.beta[0] = atoi(rcom.c_str());     //block-size in main tour
            if (com=="beta") BP.beta[0] = atoi(rcom.c_str());     //block-size in main tour
            if (com=="prob") BP.pruning_prob = atof(rcom.c_str());   
            if (com=="alpha") {
                BP.init_radius = atof(rcom.c_str());
                BP.init_mode = 'G';  //R=alpha*GH(L)
            }
            if (com=="maxtour") BP.tourlim = atoi(rcom.c_str());   
            if (com=="blocktimelimit") BP.blocktimelimit = atoi(rcom.c_str());   
            if (com=="verbose") BP.verboselevel = atoi(rcom.c_str());   
            if (com=="vl") BP.verboselevel = atoi(rcom.c_str());   

            if (com=="logfile") BP.tlname = rcom;   
            if (com=="loglevel") BP.tllevel = atoi(rcom.c_str());   


            if (com=="maxtime") BP.abortbysec = atoi(rcom.c_str());   
            if (com=="maxr") BP.abortbyr = atof(rcom.c_str());   
            if (com=="maxcost") {
                BP.extend_blocksize_initmax = atof(rcom.c_str());   
            }
            if (com=="minlogfec") BP.abortbyfec = atof(rcom.c_str());   
            if (com=="startindex") BP.startindex = atoi(rcom.c_str());   
            if (com=="endindex") BP.endindex = atoi(rcom.c_str());   
            if (com=="abortindex") BP.abortbyindex = atoi(rcom.c_str());   
            if (com=="istart") BP.startindex = atoi(rcom.c_str());   
            if (com=="iend") BP.endindex = atoi(rcom.c_str());   
            if (com=="eq-r-compute") BP.equivrcomputes = atoi(rcom.c_str());   
            if (com=="heuristicbreak") BP.heuristicbreak = true;   
            if (com=="targetapprox") BP.abortbyapprox = atoi(rcom.c_str());

            if ((com=="preprocess") || (com=="preprocess2")){
                BP.preprocess = true;
                BP.preprocess_at_least = 5;   
                BP.preprocess_strategy = 1;
                if (rcom=="l2") BP.preprocess_strategy = 2;
            }
            if (com=="preprocessonly") {
                BP.preprocessonly = true;
                BP.preprocess = true;
                BP.preprocess_max_time = atoi(rcom.c_str());
            }

            if (com=="extendblocksize") {
                BP.extend_blocksize=true;
                BP.extend_blocksizemult = 2;   //extend blocksize while expected cost <= max cost of the current tour
                //if (rcom!="") BP.extend_blocksize_initmax =  atof(rcom.c_str());   
            }
            if (com=="breakindex") {
                BP.breakindex = atoi(rcom.c_str());
            }
            if (com=="multithread") {
                BP.multithread = true;
                BP.numthreads = atoi(rcom.c_str());   
                BP.MTlimit = 50000000;      //1sec
            }

            if (com=="optpfunc") {
                BP.optimizepf = true;
                BP.optimizepf_at_least = 1000;
            }
            if (com=="modifyprobability") BP.modifyprobability=true;
        }
    }
}
#endif
