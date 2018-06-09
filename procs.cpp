#ifdef _include_idealchallenge_generator
void GenerateChallenge(mat_ZZ& L, stringmap& arg) {
    
    std::string t;
    int dim,index;
    ZZ seed;

    //Understanding options
    if (arg.find("t")==arg.end()) {
        t = arg["t"];
    }
    if ((t=="") && (arg.find("type")!=arg.end()))  {
        t = arg["type"];
    }
    dim = 0;
    if (arg.find("dim")!=arg.end()) {
        dim = atoi(arg["dim"].c_str());
    } 
    if ((dim<=0) && (arg.find("d")!=arg.end())) {
        dim = atoi(arg["d"].c_str());
    }
    index = 0;
    if (arg.find("index")!=arg.end()) {
        index = atoi(arg["index"].c_str());
    }
    if ((index<=0) && (arg.find("ind")!=arg.end() ))  dim = atoi(arg["ind"].c_str());
    seed = 0;
    if (arg.find("seed")!=arg.end()) {
        seed = to_ZZ(arg["seed"].c_str());
    }
    
    //generate
    if (t=="") {
        cout << "Type is not fixed." << endl;
        exit(0);
    } else {
        if (seed<0) {
            cout << "seed error: s=" << seed << endl;
        }
        //generate
        if ((t=="svpc") || (t=="svp") || (t=="s")) {
            if (dim<=0) {
                cout << "dimension error: d=" << dim << endl;
                exit(0);
            }
            gen_svpchallenge(L,dim,seed);
            return;
        }
        if ((t=="isvpc") || (t=="ideal") || (t=="i")) {
            if (index<=0) {
                cout << "index error: i=" << index << endl;
                exit(0);
            }
            vec_ZZ phivec;
            gen_idealsvpchallenge(L,index,seed,phivec);
            cout << "index=" << index << " dim=" << L.NumRows() << endl;
            cout << "mod. poly=" << phivec << endl; 
            return;
        }
    }
    
}
#endif

#include <lattice/gen_uni_mat.cpp>

#include "lattice/latticetools.cpp"

void GenerateMatrix(mat_ZZ& L, stringmap& arg) {
    
    std::string t;
    int dim,index;
    int bit;
    int seed;

    //Understanding options
    if (arg.find("t")!=arg.end()) {
        t = arg["t"];
    }
    if ((t=="") && (arg.find("type")!=arg.end()))  {
        t = arg["type"];
    }
    dim = 0;
    if (arg.find("dim")!=arg.end()) {
        dim = atoi(arg["dim"].c_str());
    } 
    if ((dim<=0) && (arg.find("d")!=arg.end())) {
        dim = atoi(arg["d"].c_str());
    }
    index = 0;
    if (arg.find("index")!=arg.end()) {
        index = atoi(arg["index"].c_str());
    }
    if ((index<=0) && (arg.find("ind")!=arg.end() ))  dim = atoi(arg["ind"].c_str());

    seed = 0;
    if (arg.find("seed")!=arg.end()) {
        seed = atoi(arg["seed"].c_str());
    }
    
    bit = 10;
    if (arg.find("bit")!=arg.end()) {
        bit = atoi(arg["bit"].c_str());
        if (bit<0) {
            cout << "bit error: " << bit << endl;
            exit(0);
        }
    }

    //generate
    if (t=="") {
        cout << "Type is not fixed." << endl;
        exit(0);
    }
    if (dim<=0) {
        cout << "dimension error: d=" << dim << endl;
        exit(0);
    }
    if (seed<0) {
        cout << "seed error: s=" << seed << endl;
    }
    
    //generate
    if ((t=="unimodular") || (t=="uni")) {
        gen_random_unimodular2(L,dim,seed,bit,VL0);
        return;
    }
    
}


void ApplyLLL(mat_ZZ& L, stringmap& arg) {
    
    double LLLSIback = LLLStatusInterval;
    if (arg["si"]!="") LLLStatusInterval = atof(arg["si"].c_str()); 
    
    LLL_RR(L,0.5,0,0,1);
    
    if (L.NumCols()<=450) {
        LLL_QP(L,0.999,0,0,1);
    } else {
        LLL_RR(L,0.999,0,0,1);
    }
    LLLStatusInterval = LLLSIback;
}

void RandomizeBasis(mat_ZZ& L, stringmap& arg) {

    int n = L.NumCols();
    
    if (n<=0) {
        cout << "RandomizeBasis: lattice is empty!" << endl;
        return;
    }

    mat_ZZ U;
    arg["dim"] = to_string(n);
    arg["type"] = "unimodular";
    GenerateMatrix(U,arg);
    L = U * L;
    
    
}
       
void gen_bkz_strategy_sub(bkzstrategy& BS,int sb,int eb,int n,RR& returncost) {
    
    //detect reduced level of L

    int model;
    if (n <= 400) model = modelQF;
    if (n > 400) model = modelRR;
    
    //Generating Strategy
    double A1,W;
    if (model==modelQF) {
        A1 = 1e-10;
        W = 1.5e-8;
    }
    if (model==modelRR) {
        A1 = 1e-6;    
        W = 1.5e-8;
    }

    RR cost,mincost,pcost;
    int i,j;
    int beta;
    int betaalg;
    int beta_s;

    RR* bkztime;
    bkztime = new RR[400];
    std::ostringstream* bkzstr;
    bkzstr = new std::ostringstream[400];
    for (i=0;i<400;i++) bkztime[i] = to_RR(0);
    
    
    RR mintotal;
    mintotal = -1;

    verbose = 0;
    
    for (beta=sb+1;beta<=eb;beta++) {
        //computing minimum

        //starting from LLL
        if (beta<30) {
            mincost = progressive_bkz::TimeBKZmin(n,sb,beta,A1,W,model,betaalg);
            //bkzstr[beta] << "LLL->(" << beta << "," << betaalg << ")"; 
            //bkzstr[beta] << bkzstr[sb].str() << "->" << beta << "(using BKZ-" << betaalg << ")";
            bkzstr[beta] << bkzstr[sb].str() << " " << beta << " " << betaalg;
        } else {
            mincost = to_RR(-1);
            bkzstr[beta].str("");
        }        

        beta_s = max(20,max(beta-20,sb));
        for (int lb=beta_s;lb<beta;lb++) {
            pcost = progressive_bkz::TimeBKZmin(n,lb,beta,A1,W,model,betaalg);
            if (pcost <= 0) {
                break;
            } else {
                cost = pcost + bkztime[lb];

                if ((mincost > cost) || (mincost==-1)) {
                    mincost = cost;
                    bkzstr[beta].str("");
                    //bkzstr[beta] << bkzstr[lb].str() << "{\\scriptsize $\\PROCESSBKZX{ " << beta << "}{" << betaalg << "}$} & ";
                    //bkzstr[beta] << bkzstr[lb].str() << "->" << beta  << "(using BKZ-" << betaalg << ")";
                    bkzstr[beta] << bkzstr[lb].str() << " " << beta << " " << betaalg;
                }
            }
        }
        bkztime[beta] = mincost;
        //cout << "BKZtime[" << n << "," << beta << "]=" << mincost << endl;
        //cout << "strategy[" << beta << "]: " << bkzstr[beta].str() << endl;
    }
    returncost = mincost;
    
    //cout << "strategy[" << eb << "]: " << bkzstr[eb].str() << endl;
    
    
    //Output BKZ strategy
    std::vector<std::string> bstr;
    std::string bb = bkzstr[eb].str();
    boost::algorithm::split(bstr, bb, boost::is_any_of(" "));
    quad_float* c;
    c = new quad_float[n+1];
    
    stringstream bkzss;  //bkz-strategy-string
    bkzss << "# strategy: LLL";
    for (size_t i=1;i<bstr.size();i+=2) {
        bkzss << " -> (" << bstr[i+1] << "," << bstr[i] << ")";
    }
    
    for (size_t i=1;i<bstr.size();i+=2) {
        cout << i << ": " << bstr[i] <<  ":" << bstr[i+1] << endl;
        bkzstr[i].str("");  //reuse
        
        int tbeta = atoi(bstr[i].c_str());  //target beta
        int abeta = atoi(bstr[i+1].c_str()); //beta alg.
        
        //one line strategy 
        double r,abar,alpha,prob;
        double logfec;
        r= compute_target_r(abeta);
        abar = 1.0 / lattice_tools::gaussian_heuristic_constant[abeta] * exp(-0.25*(abeta-1.0)*log(r));
        alpha = abar * ((1.0+abeta)/abeta);
        if (alpha<1.00) alpha = 1.0;
        prob = 2.0 * exp(-1.0*abeta*log(alpha));

        progressive_bkz::BKZGSsim(c,n,tbeta);  //outputs GS length with det=1    
        logfec = to_double(log(lattice_tools::FullENUMCost(c,n,INPUT_NONSQUARED)));
        
        if (n>50) {
            quad_float pred_maxcost = 2.0 * 0.25 * tbeta * ENUMCost(c,n-tbeta+1,tbeta,prob,alpha,'G',0,INPUT_NONSQUARED,1);
            bkzstr[i] << " maxcost=" << pred_maxcost;
        }
        bkzstr[i] << " bs=" << abeta << " prob=" << prob << " alpha=" << alpha << " maxtour=500";
        bkzstr[i] << " minlogfec=" << logfec << " optpfunc preprocess extendblocksize";
        bkzstr[i] << " modifyprobability";
	if (abeta<=50) bkzstr[i] << " heuristicbreak"; 
       if (abeta >= 70) {
            bkzstr[i] << " verbose=3"; 
        } else
        {
            bkzstr[i] << " verbose=1"; 
        }
    }
    BS.resize(1 + (bstr.size()-1)/2);
    j = 1;
    BS[0] = bkzss.str();
    for (size_t i=1;i<bstr.size();i+=2) {
        cout << "str[" << i << "]=" << bkzstr[i].str() << endl;
        BS[j] = bkzstr[i].str();
        j++;
    }
    

}
template <typename T> double nextparamntl(std::vector<std::vector<T> >& ct ) {
    //Return candidate of parameter to minimize cost

    sort(ct.begin(),ct.end());
    T ret;

    int i,j,n;
    n=ct.size();
    if (n==2) {
        ret = (2.0*ct[0][0]+ct[1][0])/3.0;
    }
    if (n==3) {
        ret = (ct[0][0]+2.0*ct[2][0])/3.0;
    }
    if (n>=4) {

        //detect current mincost
        j = 0;
        ret = ct[0][1];
        for (i=1;i<n;i++) {
            if (ret > ct[i][1]) {
                ret = ct[i][1];
                j=i;
            }
        }

        if (n%2==0) {
            j = min(j,n-2);
            ret = (ct[j][0] + ct[j+1][0]) * 0.5;
        } else {
            j = max(j,1);
            ret = (ct[j][0] + ct[j-1][0]) * 0.5;
        }

        //cout << "ret=" << ret << endl;

    }
    return to_double(ret);

}


void gen_shortvec_strategy_sub(bkzstrategy& BS,long int& min_M,int n,double alpha,RR& cost) {
    
    int sb = 10;
    int eb;
    bkzstrategy tempbkzstr;
    RR bkzcost;
    RR mintotal;
    mintotal = -1;
    
    quad_float* c;
    c = new quad_float[n+1];
    
    
        min_M = -1;

    for (eb=sb+1;eb<n;eb++) {

        gen_bkz_strategy_sub(tempbkzstr,sb,eb,n,bkzcost);
        
        if (bkzcost>0) {

            //TimeENUM
            double prob = 1.0 * exp(-1.0 * n * log(alpha));
            progressive_bkz::BKZGSsim(c,n,eb);  //outputs GS length with det=1    

            
            
            
            RR enumcost;
            RR total;

            //Detecting optimal # of randomized basis
            int M=1;
            int localM=min_M;
            RR localmin = to_RR(-1);

            std::vector<std::vector<RR> > ct;      //costtable
            std::vector<RR> ctt;
            ctt.resize(2);
            int lp=0;
		int bM;

	if (eb>n/2) {
            while (1) {
                if (lp==1) {
                    M = 100000;
                    cout << "maxM=" << M << endl;
                } 
                if (lp==2) M=2;
                bM = M;
                if (lp>=3) {
                    M = nextparamntl(ct);
                }

                char ef=0;
                // looks strange because i is not used in the loop
                // is this "if(ctt[i] == M)" instead of "if(ctt[0]==M)" ? 
                // for (size_t i=0;i<ctt.size();i++) {
                    if (ctt[0]==M) {
                        ef=1;
                        total=ctt[1];
                        break;
                    }
                    // }

                if (ef==0) {
                    enumcost = to_RR(ENUMCost(c,1,n,prob/1.0/M,1.05,'G',0,INPUT_NONSQUARED,1));
                    //cout << "ec=" << enumcost << endl;
                    total = M * bkzcost + to_RR(M * enumcost / 6e+7); 
                    cout << "beta=" << eb << " M=" << M << " cost=" << total << endl;
                    if ((localmin==-1) || (localmin>total)) {
                        localmin = total;
                        localM = M;
                    }
                }
                if ((lp>=3) && (abs(bM-M)<=1))  break;
                
		ctt[0] = M;
                ctt[1] = total;
                ct.push_back(ctt);
                lp++;
            }

            cout << "enumtime[" << n << "," << eb << "]=" << enumcost / 6e+7 << endl;
            cout << "totaltime[" << n << "," << eb << "]=" << total << endl;
            if ((mintotal==-1) || (mintotal > total)) {
                mintotal = total;
                min_M = localM;
                BS = tempbkzstr;
            }
            if (3 * mintotal < total) break;
		}
        }
    }
    cost = mintotal;
}


void GenerateBKZStrategy(bkzstrategy& BS,mat_ZZ &L, stringmap& arg) {
    
    std::string type;
    if (arg["type"]!="") type =arg["type"]; 
    if (arg["t"]!="") type =arg["t"];
    
    if (type=="") {
        cout << "type is not fixed." << endl;
        exit(0);
    }

    std::string outfilename =  arg["of"];
    std::string infilename =  arg["if"];

    int nt=1;   //# threads in ENUM subroutine
    if (arg["nt"]!="") nt = atoi(arg["nt"].c_str()); 
    if (arg["nthreads"]!="") nt = atoi(arg["nthreads"].c_str()); 
    if (arg["numthreads"]!="") nt = atoi(arg["numthreads"].c_str()); 
    if (nt<1) {
        cout << "Invalid setting of #Threads: " << nt << endl;
        return;
    }
    if (nt > omp_get_num_procs()) {
        cout << "***Warning*** Setting #Threads=" << nt << " > #Cores: " << omp_get_num_procs() << endl;
    }

    if (type=="bkz") {

        if (L.NumRows()==0) {
            cout << "GenerateBKZStrategy: lattice is empty!" << endl;
            return;
        }

        int n = L.NumCols();
        cout << "GenerateBKZStrategy: dim=" << n;
        if (n<=400) cout << "(simulation with QP-model)" << endl;
        if (n>400) cout << "(simulation with RR-model)" << endl;
        
        
        int sb = progressive_bkz::detectBKZLevel(L);
        cout << "CurrentBKZ-level=" << sb << ")" << endl;
        int eb;
        if (arg["eb"]=="") {
            cout << "end BKZ-level is not fixed";
            exit(0);
        } else {
            eb = atoi(arg["eb"].c_str());
        }
        RR cost;

        double curlogfec = to_double(log(lattice_tools::FullENUMCost(L)));
        double endlogfec = progressive_bkz::ComputeLogTargetFEC(n,eb);

        cout << "generating progressive bkz strategy:" << endl;
        cout << "CurrentLogFEC=" << curlogfec << " (PBKZ-level=" << sb << ")" << endl;
        cout << "TargetLogFEC=" << endlogfec << "(PBKZ-level=" << eb << ")" << endl;

        gen_bkz_strategy_sub(BS,sb,eb,n,cost);
        ofstream logstream;
        if (outfilename!="") {
            cout << "Output the strategy to " << outfilename << endl;
            logstream.open(outfilename.c_str(),ios_base::trunc);
            for (size_t i=0;i<BS.size();i++) {
                if (nt>1) BS[i] += " multithread=" + to_string(nt);
                logstream << BS[i] << endl;
            }
            logstream.close();
        } 
        return;
    }
    
    if (type=="shortvec") {
        //Fixing target radius
        double tradius=0;
        double gh = lattice_tools::LatticeGH(L);
        double alpha=0;
        
        if (arg["len"]!="")  tradius = atof(arg["det"].c_str());
        if (tradius>0) {
            cout << "Searching strategy to find a vector shorter than " << tradius << endl;
        }

        if (arg["alpha"]!="")  alpha = atof(arg["alpha"].c_str());
        if (arg["a"]!="")  alpha = atof(arg["a"].c_str());
        if (alpha>0) {
            tradius = alpha * gh; 
            cout << "Searching strategy to find a vector shorter than " << alpha << "*GH(L)" << endl;
        }

        double det;
        if (arg["det"]!="")  det = atof(arg["det"].c_str());
        if (det>0) {
            tradius = det * lattice_tools::LatticeVolumeRoot(L); 
            cout << "Searching strategy to find a vector shorter than " << det << "*det(L)^(1/n)" << endl;
        }

        cout << "Target length=" << tradius << endl;
        
        if (arg["withenum"]=="true") {

            //This routine generates a shell file
            std::string logfilename;
            int loglevel=1;
            if (arg["lf"]!="") logfilename = arg["lf"]; 
            if (arg["ll"]!="") loglevel = atoi(arg["ll"].c_str()); 
            if (loglevel<1) loglevel=1;
            
            int n = L.NumCols();
            if (n<=0) n = atof(arg["dim"].c_str());
            if (n<=0) {
                cout << "lattice dimension error:" << n << endl;
                return;
            }
            RR cost;
            long int minM;

            gen_shortvec_strategy_sub(BS,minM,n,tradius/gh,cost);

            cout << "strategy: " << endl;
            for (size_t i=0;i<BS.size();i++) {
                if (nt>1) BS[i] += " multithread=" + to_string(nt);
                cout << BS[i] << endl;
            }
            cout << "min_M=" << minM << endl;

            if (outfilename!="") {
                ofstream logstream;
                cout << "Output the strategy to:" << outfilename << endl;
                logstream.open(outfilename.c_str(),ios_base::trunc);
                for (size_t i=0;i<BS.size();i++) {
                    logstream << BS[i];
                    if (logfilename!="") {
                        logstream << " logfile=" << logfilename << " loglevel=" << loglevel;
                    }
			if (i==BS.size()-1) logstream << " heuristicbreak";
                    logstream << endl;
                }
                logstream.close();
            } 

            if (outfilename!="") {
                ofstream logstream;
                cout << "Output the shellscript to:" << outfilename << ".sh" << endl;
                logstream.open((outfilename+".sh").c_str(),ios_base::trunc);
                logstream << "#bash script to find a vector shorter than " << alpha << "*GH(L)" << endl;
                logstream << "Expected # Num. of randomized bases=" << minM << endl;
                logstream << "# Expected total CPUtime (BKZ+ENUM) in sec=" << cost << endl;
                logstream << BS[0] << endl;    //TBW

                for (int i=0;i<3*minM;i++) {
                    logstream << arg["arg0"] << " ";
                    logstream << "-if " << infilename << " ";
                    if (i!=0) {
                        logstream << "-randomizebasis -seed " << i << " ";
                    }
                    logstream << "-bkz -sf " << outfilename  << " -of " << outfilename << ".bkz.log." << i << " ";
                    logstream << "-enum -vl 3 -alpha " << alpha << " -prob " << 1.0 * exp(-1.0 * n * log(alpha)) / minM << " -of " << outfilename << ".enum." << i << " ";
                    logstream << " -lf " << outfilename << ".enum.log." << i << " ";
                    if (nt>1) logstream << " -nt " << nt;

                    logstream << endl;
                    if (arg["stopatfound"]=="true") {
                        logstream << "if [ -e " << outfilename << ".enum." << i << " ]; then" <<  endl;
                        logstream << "    exit 0" << endl;
                        logstream << "fi" << endl;
                        logstream << endl;
                    }

                }
                logstream.close();
            } 
        } else {
            //A strategy without enum (i.e., only BKZ)
            //Detecting necessary BKZ-level
            int n = L.NumCols();
            quad_float* c;
            c = new quad_float[n+1];
            int sb = progressive_bkz::detectBKZLevel(L);
            int eb = sb;
            double vr = lattice_tools::LatticeVolumeRoot(L); 
            double gh = lattice_tools::LatticeGH(L); 
            cout << "Current BKZ-level=" << sb << " (|b1|=" << LengthOf(L[0]) << endl;
            while (1) {
                cout << "BKZ-" << eb << ": ";
                progressive_bkz::BKZGSsim(c,n,eb);  //outputs GS length with det=1    
                cout << "expected |b1|=" << c[1] * vr << "=" << c[1] << "*det(L)^(1/n)="<< c[1]*vr/gh << "*GH(L)" <<endl;
                if (c[1]*vr < tradius) break;
                eb++;
                if (eb==n) break;
            }
            delete [] c;

            cout << "generating BKZ-" << sb << " to BKZ-" << eb << endl;
            RR cost;
            gen_bkz_strategy_sub(BS,sb,eb,n,cost);
            ofstream logstream;
            if (outfilename!="") {
                cout << "Output the strategy to " << outfilename << endl;
                logstream.open(outfilename.c_str(),ios_base::trunc);
                for (size_t i=0;i<BS.size();i++) {
                    if (nt>1) BS[i] += " multithread=" + to_string(nt);
                    logstream << BS[i] << endl;
                }
                logstream.close();
            } 
            return;
            
                
            
            
        }        
        
        return;
        
        
    }
    

}

void ApplyBKZ(mat_ZZ& L, bkzstrategy& BS, stringmap& arg) {
    
    if (L.NumCols()<=1) {
        cout << "Lattice is empty!" << endl;
        return;
    }
    //ApplyLLL(L,arg);
    int n = L.NumCols();
    
    //Read strategy file
    std::string strfilename;  
    std::vector<std::string> bkzstrategy;
    strfilename = arg["sf"];
    
    if ((strfilename=="") || (FileExists(strfilename)==false)) {
        //BS.resize(0);
    } else {
        //Read 
        ifstream logstream;
        logstream.open(strfilename.c_str());

        std::string t;
        while (getline(logstream,t)) {
            //cout << "read: " << t << endl;
            BS.push_back(t);
        }
        logstream.close();
    }


    int nt=1;   //# threads in ENUM subroutine
    if (arg["nt"]!="") nt = atoi(arg["nt"].c_str()); 
    if (arg["nthreads"]!="") nt = atoi(arg["nthreads"].c_str()); 
    if (arg["numthreads"]!="") nt = atoi(arg["numthreads"].c_str()); 
    if (nt<1) {
        cout << "Invalid setting of #Threads: " << nt << endl;
        return;
    }
    if (nt > omp_get_num_procs()) {
        cout << "***Warning*** Setting #Threads=" << nt << " > #Cores: " << omp_get_num_procs() << endl;
    }

    std::string logfilename;
    int loglevel=1;
    if (arg["lf"]!="") logfilename = arg["lf"]; 
    if (arg["ll"]!="") loglevel = atoi(arg["ll"].c_str()); 
    if (loglevel<1) loglevel=1;
    
    if (arg["fix"]!="") {    
        //fixed strategy given in the arguments;
        std::stringstream opts;
        int beta;
        double alpha,prob;
        int maxtour;
        
        beta = 30;  alpha = 1.0;    prob = 1.0; //default values
        maxtour = 500;

        if (arg["maxtour"]!="") maxtour = atoi(arg["maxtour"].c_str());
        if (arg["beta"]!="") beta = atoi(arg["beta"].c_str());
        if (arg["alpha"]!="") alpha = atof(arg["alpha"].c_str());

        opts << "bs=" << beta << " prob=" << prob << " alpha=" << alpha << " maxtour=" << maxtour;
        if (logfilename!="") opts << " logfile=" << logfilename << " loglevel=" << loglevel;

        if (arg["preprocess"]!="") opts << " preprocess=" << arg["preprocess"];
        if (arg["vl"]!="") opts << " vl=" << arg["vl"];
        
        cout << "Fixed parameter BKZ with options=[" << opts.str() << "]" << endl;

        verbose=1;
        LLLStatusInterval=10;
        
        //BKZ_FP(L,0.999,beta,0,0,1);
        cout << lattice_tools::FullENUMCost(L) << endl;
        progressive_bkz::fixedparamBKZ(L,opts.str()); 
        
    } else {
    
        if (BS.size()==0) {    
            cout << "BKZ strategy file not found: generating one-by-one strategy" << endl;

            int sb;
            if (arg["sb"]!="") {
                sb = atoi(arg["sb"].c_str());
            } else {
                sb = progressive_bkz::detectBKZLevel(L);
            }
            int eb =  atoi(arg["eb"].c_str());
            
            int beta;

            double r,abar,alpha,prob;
            std::stringstream opts;
            quad_float* c;
            c = new quad_float[n+1];
            double logfec;

            for (beta=sb;beta<=min(eb,n);beta++) {
                r= compute_target_r(beta);
                abar = 1.0 / lattice_tools::gaussian_heuristic_constant[beta] * exp(-0.25*(beta-1.0)*log(r));
                alpha = abar * ((1.0+beta)/beta);
                if (alpha<1.00) alpha = 1.0;
                prob = 2.0 * exp(-1.0*beta*log(alpha));

                int is=1;
                int ie=n;
                int rn=n; //really processed dimension
                if (arg["ignoreflat"]=="true") {
                    //For special lattices such as q-ary
                    cout << "Detect non-flat range: ";
                    lattice_tools::LatticeGSLengths(L,c);

                    for (int i=1;i<=n;i++) cout << c[i] << " ";
                    cout << endl;
                    is = 1;
                    while ((is<=n) && (0.999<c[is]/c[1]) && (c[is]/c[1]<1.001)) is++;
                    ie = n;
                    while ((ie>=is) && (0.999<c[ie]/c[n]) && (c[ie]/c[n]<1.001)) {
                        //cout << ie << " " << c[ie] << " " << c[n] << endl;
                        ie--;
                    }
                    //cout << ie << " " << c[ie] << " " << c[n] << endl;
                    
                    is = max(is-3,1);
                    ie = min(ie+3,n);
                    cout << is << "-" << ie << endl;
                    rn = ie-is+1;
                }


                progressive_bkz::BKZGSsim(c,rn,beta-1);  //outputs GS length with det=1    
                logfec = to_double(log(lattice_tools::FullENUMCost(c,rn,INPUT_NONSQUARED)));
                //logfec = ComputeLogTargetFEC(n,beta-1);

                opts.str("");
                if (rn>50) {
                    quad_float pred_maxcost = 2.0 * 0.25 * beta * ENUMCost(c,rn-beta+1,beta,prob,alpha,'G',0,INPUT_NONSQUARED,1);
                    opts << " maxcost=" << pred_maxcost;
                }

                opts << " bs=" << beta << " prob=" << prob << " alpha=" << alpha << " maxtour=500";
                opts << " minlogfec=" << logfec << " optpfunc preprocess extendblocksize";
                if (nt>1) opts << " multithread=" << nt;
                opts << " modifyprobability ";
		if (beta<=50) opts << " heuristicbreak";
                opts << " istart=" << is << " iend=" << ie;
                if (logfilename!="") {
                    opts << " logfile=" << logfilename << " loglevel=" << loglevel;
                }
                if (beta >= 70) {
                    opts << " verbose=3"; 
                } else
                {
                    opts << " verbose=1"; 
                }
                
                cout << opts.str() << endl;
                progressive_bkz::fixedparamBKZ(L,opts.str());
            }

        } else {
            //When the strategy is given
            for (size_t i=0;i<BS.size();i++) {
                if (BS[i][0]!='#') {    //skip comment lines
                    if (nt>1) BS[i] = BS[i] + " multithread=" + to_string(nt); 
                    if (logfilename!="") {
                        BS[i] += " logfile=" + logfilename + " loglevel=" + to_string(loglevel);
                    }
                    cout << BS[i] << endl;
                    progressive_bkz::fixedparamBKZ(L,BS[i]);
                }
            }
        }
        
    }
}

void PBKZ_QP(mat_ZZ& L,int sb,int eb) {
    bkzstrategy BS; //fake;
    stringmap arg;
    arg["sb"]=to_string(sb);
    arg["eb"]=to_string(eb);
    ApplyBKZ(L, BS, arg);
}


void GeneratePruningFunction(pruningfunction& ret,mat_ZZ &L, stringmap& arg) {
    //Generating normalized pruning function pf[0..n-1] in monotonically increasing order
    //The values are squared
    //Target is input lattice or Gram-Schmidt length (given by -gsf option)

    int n = L.NumCols();
    char gsflag = 0;
    quad_float* c;

    if (n<=0) {
        //try the -gsf option
        if (arg["gsf"]=="") {
            return;
        }
        if (FileExists(arg["gsf"])==false) return;
        gsflag = 1;
        //read the file
        ifstream gstream;
        gstream.open( arg["gsf"].c_str(),ios_base::in);
        c = new quad_float[1000];   
        int k = 1;
        while (gstream.eof()==false) {
            std::string line;
            gstream >> line;
            if (line.length()>0) {
                if ((line[0]!='#') && (atof(line.c_str())>0)) {
                    c[k] = atof(line.c_str());
                    c[k] *= c[k];   //squared
                    k++;
                }
            }
        }
        n = k-1;
    } else {
        c = new quad_float[n+1];
        lattice_tools::LatticeGSLengths(L,c);
    }
    
    //default parameters
    double prob = 1.0;
    double alpha = 1.0;
    double radius;
    int vl = 1;
    int optimize=0;
    std::string outfilename;
    std::string pffilename;
    std::string logfilename;
    

    //Reading parameters
    if (arg["p"]!="") prob = atof(arg["p"].c_str());
    if (arg["prob"]!="") prob = atof(arg["prob"].c_str());
    if ((prob<=0) || (prob>1)) {
        cout << "probability error: p=" << prob << endl;
        exit(0);
    }
    
    double gh;
    if (gsflag==0) gh=lattice_tools::LatticeGH(L);
    if (gsflag==1) gh=lattice_tools::LatticeGH(c,n,INPUT_SQUARED);

    if (arg["a"]!="") alpha = atof(arg["a"].c_str());
    if (arg["alpha"]!="") alpha = atof(arg["alpha"].c_str());
    if (alpha<=0) {
        cout << "alpha error: a=" << alpha << endl;
        exit(0);
    }
    radius = alpha * gh;
    if (arg["radius"]!="") {
        radius = atof(arg["radius"].c_str());
        alpha = radius / gh;
    }
    
    //output level
    if (arg["vl"]!="") vl = atoi(arg["vl"].c_str());
    if (arg["ol"]!="") vl = atoi(arg["ol"].c_str());

    //filename
    if (arg["of"]!="") outfilename = arg["of"];
    if (arg["pf"]!="") pffilename = arg["pf"];
    if (arg["lf"]!="") logfilename = arg["lf"];

    //Optimizing pruning function
    if (arg["optimize"]!="") optimize = atoi(arg["optimize"].c_str());

    if (vl>=1) cout << "Generate pruning function dim=" << n << " radius=" << radius << "=" << alpha << "*GH, prob=" << prob << endl; 
    if (vl>=1) {
        cout << "GS-lengths=[";
        for (int i=1;i<=n;i++) cout << sqrt(c[i]) << " ";
        cout << "]" << endl;
    }    
    double* pf;
    pf = new double[n+1];
    
    if ((pffilename!="") && (FileExists(pffilename)==true)) {
        loadarray(pf-1,n,pffilename);
        if (vl>=2) {
            cout << "input Pruning function=[";
            for (int i=0;i<n;i++) cout << sqrt(pf[i])*radius << " ";
            cout << "]" << endl;
        }
    } else {
        set_approximate_extreme_pruning(pf,c+1,n,prob,0);
        pruning_func::invert(pf,0,n-1);
        if (vl>=2) {
            cout << "approx. Pruning function=[";
            for (int i=0;i<n;i++) cout << sqrt(pf[i])*radius << " ";
            cout << "]" << endl;
        }
    }
    
    RR cost;
    if (vl>=1) {
        //display properties
        cout << "target=" << prob << endl;
        cout << "LowerBound=" << pruning_func::Rigid_lower_prob(pf-1,n) << endl;
        cout << "UpperBound=" << pruning_func::Rigid_upper_prob(pf-1,n) << endl;
        if (vl>=3) cout << "Approximate=" << pruning_func::Approx_prob(pf-1,n) << endl;
        progressive_bkz::sharememalloc();
        quad_float* V = progressive_bkz::bkz_share_quad_float2[0];
        pruning_func::init_pruning_func();
        pruning_func::computeVfromR(V,pf-1,n);
        cost = to_RR(pruning_func::Rigid_upper_cost_pcvf(V,c,n,radius));
        cout << "Cost=" << cost << endl;
    }
    
    if (optimize>0) {
        cost = pruning_func::optimize_pruning_function(pf-1,radius,c,n,prob,vl,optimize,orderMI);
        cout << "Pruning function=[";
        for (int i=0;i<n;i++) cout << sqrt(pf[i])*radius << " ";
        cout << "]" << endl;
        cout << "Optimized Cost=" << cost << endl;
    }    

    if (logfilename!="") {
        ofstream lstream;
        lstream.open(logfilename.c_str(),ios_base::trunc);
        
        lstream << "#Log of GeneratePruningFunction()" << endl;
        lstream << "lattice_dimension=" << n << endl;
        lstream << "target_radius=" << radius << "=" << alpha << "*GH" << endl;
        lstream << "target_probability=" << prob << endl;
        lstream << "GS-lengths=[";
        for (int i=1;i<=n;i++) lstream << sqrt(c[i]) << " ";
        lstream << "]" << endl;

        lstream << "pruning func to bound |pi_i(v)| =[";
        for (int i=1;i<=n;i++) lstream << sqrt(pf[i])*radius << " ";
        lstream << "]" << endl;
        
        lstream << "Cost=" << cost << endl;
        
    }

    //returned value
    ret.resize(n);
    for (int i=0;i<n;i++) ret[i] = pf[i];

    if (outfilename!="") {
        savearray((double*)ret.data()-1,n,outfilename);
    }
}

void SimulateGSLengths(stringmap& arg) {

    //Simulating Gram-Schmidt lengths
    quad_float* c;
    int dim=0;
    double gsa_r=0;
    int pbkz=0;
    int hkz=0;

    int vl=0;
    char sqflag=0;  //1 if the output is squared values
    quad_float normdet; //normalize determinant as this value
    quad_float normgh; //normalize gaussian heuristic as this value
    std::string outfilename;
    
    if (arg["dim"]!="") dim = atoi(arg["dim"].c_str());
    if (arg["hkz"]!="") hkz = 1;
    if (arg["gsa"]!="") gsa_r = atof(arg["gsa"].c_str());
    if (arg["pbkz"]!="") pbkz = atoi(arg["pbkz"].c_str());
    if (arg["det"]!="") normdet = conv<quad_float>(atof(arg["det"].c_str()));
    if (arg["gh"]!="") normgh = conv<quad_float>(atof(arg["gh"].c_str()));
    if (arg["vl"]!="") vl = atoi(arg["vl"].c_str());
    if (arg["of"]!="") outfilename = arg["of"];
    
    if (dim<=1) {
        cout << "Dimension error: " << dim << endl;
        exit(0);
    }
    c = new quad_float[dim+1];
    if (gsa_r>0) {
        //Generating Schnorr's GSA lengths
        c[1] = 1.0;
        for (int i=2;i<=dim;i++) c[i] = c[i-1] * gsa_r;
        pbkz=0;
    }
    if (pbkz>1) {
        //Generating simulated PBKZ basis
        //output is non-squared
        progressive_bkz::BKZGSsim(c,dim,pbkz);
    }
    
    if (hkz!=0) {
        //Generating simulated HKZ basis 
        //that satisfies |b*i|=GH(B_i) for all i
        //the last components are from the table in Chen-Nguyen's BKZ 2.0 paper
        gen_hkzlattice(c,dim,1.00);
    }
    
    //normalization
    if (normgh>0) {
        lattice_tools::initialize();
        normdet = normgh / to_double(lattice_tools::gaussian_heuristic_constant[dim]);
    }
    if (normdet>0) {
        double ldet=0;
        for (int i=1;i<=dim;i++) ldet += to_double(log(c[i]));
        ldet = exp(ldet / dim); //vol^(1/n)
        for (int i=1;i<=dim;i++) c[i] = c[i] / ldet * normdet;
    }
    
    if (sqflag==1) {
        for (int i=1;i<=dim;i++) c[i] *= c[i];
    }
    
    //display
    double rootdet = lattice_tools::LatticeVolumeRoot(c,dim,INPUT_NONSQUARED);
    double gh = lattice_tools::LatticeGH(c,dim,INPUT_NONSQUARED);
    if (vl>=1) {
        cout << "GS-lengths";
        if (sqflag==1) cout << "(squared)";
        cout << "=[";
        for (int i=1;i<=dim;i++) cout << c[i] << " ";
        cout << "]" << endl;
        cout << "det^(1/dim)=" << rootdet << endl;
        cout << "GH=" << gh << endl;
    }

    
    if (outfilename!="") {
        ofstream gstream;
        gstream.open(outfilename.c_str(),ios_base::trunc);
        
        gstream << "#Simulated Gram-Schmidt lengths, dim=" << dim;
        if (sqflag==1) gstream << "(squared)";
        gstream << endl;
        gstream << "#Type=";
        if (gsa_r>0) gstream << "GSA(r=" << gsa_r << ")";
        if (pbkz>0) gstream << "progressive_BKZ(beta=" << pbkz << ")";

        gstream << endl;
        gstream << "#det^(1/dim)=" << rootdet << " gh=" << gh << endl;
        for (int i=1;i<=dim;i++) gstream << c[i] << endl;
        gstream.close();
    }
    

}

void LatticeENUM(mat_ZZ &L,pruningfunction& argpf, stringmap& arg) {

    int n = L.NumCols();
    if (n<=0) return;
    
    mat_ZZ foundmat;
    //default parameters
    double prob = 1.0;
    double alpha = 1.0;
    int elim = 100000000;
    int vl = 1;
    int parallel=1;
    std::string outfilename;
    std::string pffilename;
    quad_float* c;
    quad_float** mu;
    c = new quad_float[n+1];
    mu = new quad_float*[n+1];
    for (int i=0;i<n+1;i++) mu[i] = new quad_float[n+1];
    GramSchmidt(L,mu,c,n,1,1);
    lattice_enum::initialize();
    
    double gh = lattice_tools::LatticeGH(L);
    
    //Reading parameters
    if (arg["p"]!="") prob = atof(arg["p"].c_str());
    if (arg["prob"]!="") prob = atof(arg["prob"].c_str());
    if ((prob<=0) || (prob>1)) {
        cout << "probability error: p=" << prob << endl;
        exit(0);
    }
    
    if (arg["a"]!="") alpha = atof(arg["a"].c_str());
    if (arg["alpha"]!="") alpha = atof(arg["alpha"].c_str());
    if (alpha<=0) {
        cout << "alpha error: a=" << alpha << endl;
        exit(0);
    }
    double radius = alpha * lattice_tools::LatticeGH(L);
    if (arg["radius"]!="") {
        radius = atof(arg["radius"].c_str());
        alpha = radius / gh;
    }

    //output level
    if (arg["vl"]!="") vl = atoi(arg["vl"].c_str());
    if (arg["ol"]!="") vl = atoi(arg["ol"].c_str());

    if (arg["parallel"]!="") parallel = atoi(arg["parallel"].c_str());
    if (arg["th"]!="") parallel = atoi(arg["th"].c_str());
    if (arg["nt"]!="") parallel = atoi(arg["nt"].c_str());

    if (arg["of"]!="") outfilename = arg["of"];
    if (arg["pf"]!="") pffilename = arg["pf"];

    std::string logfilename =  arg["lf"];

    //pruning function
    double* pf;
    pf = new double[n+1];
    if ((int)argpf.size()==n) {
        for (int i=0;i<n;i++) pf[i] = argpf[i];
    } else 
    if (pffilename!="") {
        //values are non-squared value
        loadarray(pf-1,n,pffilename);
     }  else {
        set_approximate_extreme_pruning(pf,c+1,n,prob,0);
        pruning_func::invert(pf,0,n-1);
     }
    if (vl>=2) {
        cout << "Pruning function=[";
        for (int i=0;i<n;i++) cout << pf[i] << " ";
        cout << "]" << endl;
    }
    
    if (vl>=1) cout << "ENUM dim=" << n << " radius=" << radius << "=" << alpha << "*GH, prob=" << prob << endl; 
    if (vl>=2) cout << "#threads=" << parallel << endl;
    
    //preprocess the basis	
    if (arg["preprocess"]!="") {
   	int sec = atoi(arg["preprocess"].c_str());
        std::stringstream opts;
        
        opts << "bs=" << n << " alpha=" << alpha << " prob=" << prob << " vl=" << vl;
        opts << " preprocessonly=" << sec;
        cout << "Options=[" << opts.str() << endl;

        progressive_bkz::fixedparamBKZ(L,opts.str());
        GramSchmidt(L,mu,c,n,1,1);  //mu[1..n][1..n] and c[1..n]
        
        set_approximate_extreme_pruning(pf,c+1,n,prob,0);
        pruning_func::invert(pf,0,n-1);

    }

	

    if (arg["optimize"]!="") {

   	 int optimizesec = atoi(arg["optimize"].c_str());
         
        pruning_func::init_pruning_func();
        RR cost = pruning_func::optimize_pruning_function(pf-1,radius,c,n,prob,vl,optimizesec,orderMI);
        if (vl>=3) {
            cout << "Pruning function=[";
            for (int i=0;i<n;i++) cout << pf[i] << " ";
            cout << "]" << endl;
        }
        cout << "Optimized Cost=" << cost << endl;
    }


    int buff=maxfoundvec;       //maxfoundvec is defined in the header file
    double* foundlen = lattice_enum::enum_share_double[9];  //9 is alloced for foundlen
    int** foundvec = lattice_enum::enum_share_foundvec;
    int findex=0;   //start index of found vector
    int slide=0;
    int enummode = enum_mode_all_vectors;
    int finishmode = finish_mode_exact;
    pruning_func::invert(pf,0,n-1);

    char alg=0; //0=find short vectors, 1=find close vectors
    
    if (arg["enummode"]=="shortest") enummode |= enum_mode_find_shortest;
    if (arg["enummode"]=="all") enummode |= enum_mode_all_vectors;
    if (arg["enummode"]=="abort") enummode |= enum_mode_find_abort;
    if (arg["enummode"]=="shortproj") enummode |= enum_mode_find_short_projection;
    if (arg["enummode"]=="excepttrivial") enummode |= enum_mode_except_trivial;
    if (arg["finishmode"]=="nowaste") finishmode |= finish_mode_nowaste;
    if (arg["finishmode"]=="exact") finishmode |= finish_mode_exact;

    if (arg["maxbuff"]!="") buff = atoi(arg["maxbuff"].c_str());

    vec_ZZ tv;
    if (arg["cvp"]!="") {
        alg = 1;
        std::string vecname = arg["cvp"];
        if (FileExists(vecname)==true) {
            LoadElement(tv,vecname);
        }
        if (vl>=1) cout << "target_vector=" << tv << endl;
    }
    
    double es;
    es = gettimeofday_sec();   //timer
    double startcputime = clock();

    if (alg==0) {
        //SVP
        if (parallel>=2) {
            lattice_enum::ENUMCoremt(foundvec,foundlen,buff,findex,slide,1,n,elim,radius*radius,pf,mu,c,vl,parallel,enummode,finishmode);
        } else {
            lattice_enum::ENUMCore(foundvec,foundlen,buff,findex,slide,1,n,elim,radius*radius,pf,mu,c,vl,enummode);
        }
    } else 
    if (alg==1) {
        //CVP
        //target_mu
        quad_float* target_mu;      //represent so that target = \sum_{i=0}^{n-1} target_mu[i].b*[i]
        target_mu = new quad_float[n+1];
        lattice_tools::alpha_decompose(target_mu,tv,L,1);   //decompose tv=\sum mu[i]. b*[i]]
        if (vl>=3) {
            cout << "target_mu=[";
            for (int i=0;i<n;i++) cout << target_mu[i] << " ";
            cout <<"]" <<  endl;
        }
        for (int i=n;i>=1;i--) target_mu[i] = target_mu[i-1]; 
        lattice_enum::ENUMCVPCoremt(foundvec,foundlen,buff,findex,slide,1,n,elim,radius*radius,pf,target_mu,mu,c,vl,parallel,enummode,finishmode);
    }        
    if (vl>=1) cout << "# found vector=" << findex << endl;

    convert_intcoefftomat_zz(foundmat,L,foundvec,findex,0);
    if (vl>=3) {  
        for (int i=0;i<findex;i++) {
            //cout << sqrt(foundlen[i]) << " (doublecheck:" << LengthOf(foundmat[i]) <<  endl;
            cout << foundmat[i] << endl;
        }
    }
    if (outfilename!="") {
        if (foundmat.NumCols()>0) 
            if (foundmat.NumRows()>0) {
                SaveLattice(foundmat,outfilename);
                if (vl>=1) {
                    cout << "found " << foundmat.NumRows() << " vectors are saved to"  << outfilename.c_str() << endl;
                }
            } 
    }
    if (logfilename!="") {
        //Output log file
        ofstream logstream;
        logstream.open(logfilename.c_str(),ios_base::app);
        logstream << timestr() << "\t"; //Current time
        logstream << n << "\t"; //Dimension
        logstream << radius << "\t"; //Searching radius
        logstream << alpha << "\t"; //Radius/GH(L)
        logstream << prob << "\t"; //Pruning probability
        logstream << parallel << "\t";	//# threads
        logstream << gettimeofday_sec() - es << "\t";   //Elaplsed time
        logstream << (clock()-startcputime) / CLOCKS_PER_SEC << "\t";   //CPUTime time
        logstream << lattice_enum::current_processed << "\t";   //# processed nodes
        logstream << lattice_enum::current_enum_speed << "\t";   //# processed nodes / Elaplsed time
        logstream.close();
    }
}

void ReorderBasis(mat_ZZ &L, stringmap& arg) {

    int i;
    int n = L.NumRows();
    
    if (arg["type"]=="inverse") {
        for (i=0;i<n/2;i++) swap(L[i],L[n-i-1]);
    }
}

void DisplayLatticeProperty(mat_ZZ &L, stringmap& arg) {
    //TBD
    std::string outfilename;
    int vl = 0;
    int i;
    int n;
    lattice_tools::initialize();
    
    if (arg["vl"]!="") vl = atoi(arg["vl"].c_str());
    
    if (vl>=1) cout << L << endl;
    n = L.NumCols();
    cout << "LatticeDimension=" << n << endl;

    cout << "GS-lengths=[";
    quad_float* c;
    c = new quad_float[n+1];
    lattice_tools::LatticeGSLengths(L,c);
    for (i=1;i<=n;i++) cout << sqrt(c[i]) << " ";
    cout << "]" << endl;
    
    cout << "FEC=" << lattice_tools::FullENUMCost(c,n) << endl;
    cout << "GH=" << lattice_tools::LatticeGH(c,n,INPUT_SQUARED) << endl;
    cout << "det^(1/n)=" << lattice_tools::LatticeVolumeRoot(c,n,INPUT_SQUARED) << endl;
     
    
}


void processcommand(std::vector<std::string>& args) {
    
    mat_ZZ L;
    bkzstrategy BS;
    pruningfunction pf;
    
    std::string command;
    std::string infilename,outfilename;
    stringmap addarg; // to be used in subroutines

    size_t i;
    
    //Switch
    i = 1;
    stringmap currentcommand; // to be used in subroutines
    do {
        
        if ((i==args.size()) || (args[i].substr(0,4)=="com=")) {
            //execute stacked command 
    #ifdef __develop
            processotherdevel(L,currentcommand,infilename,outfilename);
    #endif
            
            
    #ifdef _include_svpchallenge_generator
            if (currentcommand["command"]=="genbasis") {
                GenerateChallenge(L,currentcommand);
                if (outfilename!="") {
                    SaveLattice(L,outfilename);
                } else {
                    cout << L << endl;
                }
                infilename = "";
                outfilename = "";
            }
    #endif

            if (currentcommand["command"]=="genmatrix") {
                GenerateMatrix(L,currentcommand);
                if (outfilename!="") {
                    SaveLattice(L,outfilename);
                } else {
                    cout << L << endl;
                }
                infilename = "";
                outfilename = "";
            }

            if (currentcommand["command"]=="lll") {
                if (infilename!="") LoadLattice(L,infilename);
                ApplyLLL(L,currentcommand);
                if (outfilename!="") {
                    SaveLattice(L,outfilename);
                } else {
                    cout << L << endl;
                }
                infilename = "";
                outfilename = "";
            }

            if (currentcommand["command"]=="reorder") {
                if (infilename!="") LoadLattice(L,infilename);
                ReorderBasis(L,currentcommand);
                if (outfilename!="") {
                    SaveLattice(L,outfilename);
                } else {
                    cout << L << endl;
                }
                infilename = "";
                outfilename = "";
            }

            if (currentcommand["command"]=="bkz") {
                if (infilename!="") LoadLattice(L,infilename);
                ApplyBKZ(L,BS,currentcommand);
                if (outfilename!="") {
                    SaveLattice(L,outfilename);
                } else {
                    cout << L << endl;
                }
                infilename = "";
                outfilename = "";
            }

            if (currentcommand["command"]=="randomizebasis") {
                if (infilename!="") LoadLattice(L,infilename);
                RandomizeBasis(L,currentcommand);
                if (outfilename!="") {
                    SaveLattice(L,outfilename);
                } else {
                    cout << L << endl;
                }
                infilename = "";
                outfilename = "";
            }

            if (currentcommand["command"]=="genstrategy") {
                currentcommand["arg0"] = args[0];
                if (infilename!="") {
                    LoadLattice(L,infilename);
                    currentcommand["if"] = infilename;
                }
                GenerateBKZStrategy(BS,L,currentcommand);
                infilename = "";
                outfilename = "";
            }
            
            if (currentcommand["command"]=="enum") {
                if (infilename!="") LoadLattice(L,infilename);
                LatticeENUM(L,pf,currentcommand);
                infilename = "";
                outfilename = "";
            }
            if (currentcommand["command"]=="pfunc") {
                if (infilename!="") LoadLattice(L,infilename);
                GeneratePruningFunction(pf,L,currentcommand);
                infilename = "";
                outfilename = "";
            }
            if (currentcommand["command"]=="simgslength") {
                SimulateGSLengths(currentcommand);
                infilename = "";
                outfilename = "";
            }

            if (currentcommand["command"]=="disp") {
                if (infilename!="") LoadLattice(L,infilename);
                DisplayLatticeProperty(L,currentcommand);
                infilename = "";
                outfilename = "";
            }


            currentcommand.clear();
            if (i<args.size()) {
                currentcommand["command"] = args[i].substr(4,args[i].length());
                //cout << "nextcommand = " << currentcommand["command"] << endl; 
            }
        } else {
            if (i==args.size()-1) {
                //last one argument
                currentcommand[args[i]] = "true";
                if (args[i][0]=='-') {
                    currentcommand[args[i].substr(1,args[i].length()-1)] = "true";
                }
                //cout << "reg: " << args[i] << " true";
            } else {

                //file
                if ((args[i]=="-f") || (args[i]=="-if")) {
                    infilename = args[i+1];
                }
                if (args[i]=="-of") {
                    outfilename = args[i+1];
                }

                //other commands
                if ((args[i+1][0]=='-') && (args[i+1].find(" ")==std::string::npos )) {
                    //next one is not argument
                    currentcommand[args[i].substr(1,args[i].length())] = "true" ;
                    //cout << "reg: " << args[i] << " true" << endl;
                } else {
                    currentcommand[args[i].substr(1,args[i].length())] = args[i+1];
                    //cout << "reg: " << args[i] << " " << args[i+1] << endl;
                    i++;
                }
            }
        }
        i++;
    } while (i<=args.size());
    
    
    
    
}
