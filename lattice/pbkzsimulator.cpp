#ifndef _inc_pbkz_simulator
#define _inc_pbkz_simulator

#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#define modelQF 0x01
//QF=quad_float is for small dimensions < 400

#define modelRR 0x02
//RR is for large dimensions > 400

namespace progressive_bkz {
    double total_simulating_time = 0;
}

namespace progressive_bkz {

    
    void simulateLLL(quad_float* c,int n) {

        //output is c[1...n] not squared
        double delta = 1.022;
        double q = exp(-4.0*n/(n-1.0) * log(delta));
        int i;
        c[1] = 1.0;
        for (i=2;i<=n;i++) c[i] = c[i-1] * q;

    }

    void BKZsim_first(quad_float* c,int n,int beta) {
       //first round 
       //outputs GS length with det=1 
        // c[i]=|b*i| (i=1,...,n)

       int i,j;
       double alpha2;
       
       if (bkz_verbose>=1) cout << "Gram-Schmidt sim (first step): n=" << n << " beta=" << beta << endl; 
       double r = compute_target_r(beta);
       if (bkz_verbose>=1) cout << "r=" << r << endl;

       double abar,alpha;
       abar = 1.0 / lattice_tools::gaussian_heuristic_constant[beta] * exp(-0.25*(beta-1.0)*log(r));        //expected finding radius
       alpha = abar * ((1.0+beta)/beta);   //searching radius
       if (bkz_verbose>=1) cout << "alpha=" << alpha << " abar=" << abar << endl;

       c[1] = 1.0;
       
        for (i=2;i<=n;i++) {
            //Generating |b*i| s.t. |b*i|=alpha*GH(L')

            int bs = min(i,beta);

            alpha2 = alpha * (bs)/(1.0+bs);
            if (bs<=50) alpha2 = max(alpha2,lattice_tools::smallghconst[bs]);

            c[i] = 0;
            for (j=i-bs+1;j<i;j++) {
                c[i] += 1.0*log(c[j]) / (bs-1.0);
            }
            c[i] -= to_quad_float(1.0*log(VolumeBall(bs,1.0))/(bs-1.0));
            c[i] += bs/(bs-1.0)*log(alpha2);
            c[i] = exp(c[i]);
        }

        for (i=1;i<=n/2;i++) swap(c[i],c[n-i+1]);
        for (i=2;i<=n;i++) c[i] /= c[1];
        c[1] = 1.0;
   }    //end of BKZ20GSsim_first()

   void BKZsim_second(quad_float* c,int n,int beta) {
       //outputs GS length with det=1 
       // c[i]=|b*i| (i=1,...,n)
       //Assume c[] is the temporal simulation
       
       int i,j;
       double alpha2;
       
       if (bkz_verbose>=1) cout << "Gram-Schmidt sim (second step): n=" << n << " beta=" << beta << endl; 
       double r = compute_target_r(beta);
       //cout << "r=" << r << endl;

       double abar,alpha,prob;
       abar = 1.0 / lattice_tools::gaussian_heuristic_constant[beta] * exp(-0.25*(beta-1.0)*log(r));        //expected finding radius
       alpha = abar * ((1.0+beta)/beta);   //searching radius
       if (bkz_verbose>=1) cout << "alpha=" << alpha << " abar=" << abar << endl;

       prob = min(1.0,2.0*exp(-beta*log(alpha)));
       if (bkz_verbose>=1) cout << "prob=" << prob << endl; 
       //cout << "max_enum_cost=" << ENUMCost(c,1,beta,prob,alpha,'G',0,INPUT_NONSQUARED) << endl;
       RR maxcost = to_RR(ENUMCost(c,n-beta+1,beta,prob,alpha,'G',0,INPUT_NONSQUARED,1));
       if (bkz_verbose>=1) cout << "sim_max_enum_cost=" << maxcost << endl;
       if (bkz_verbose>=1) cout << "Modify last indexes" << endl;
       
       //computing modified probability
       double* malpha;
       malpha = new double[beta+1];
       double tprob = 1.0;
       double talpha;
       for (i=2;i<beta;i++) {
           if ( lattice_tools::FullENUMCost(c+n-i,i,INPUT_NONSQUARED) < 2 * maxcost) {
               tprob = 1.0;
               talpha = 1.0;
           } else {
               tprob = min(1.0,tprob*1.1);
               do {
                   tprob /= 1.1;
                   talpha = max(1.0,exp(1.0 / i * log(2.0/tprob)));
               } while ( (to_RR(ENUMCost(c,n-i+1,i,tprob*1.1,talpha,'G',0,INPUT_NONSQUARED,1)) > 2*maxcost) && (tprob>prob) );
           }
           tprob = max(prob,tprob);
           if (bkz_verbose>=1) cout << "modprob[" << i << "]=" << tprob << " mod_alpha=" << talpha << endl;
           malpha[i] = talpha;
       }

       //Re-compute GS-basis with modified last probabilities
       c[1] = 1.0;
        for (i=2;i<=n;i++) {
            //Generating |b*i| s.t. |b*i|=alpha*GH(L')

            int bs = min(i,beta);
            if (bs<beta) {
                alpha2 = 1.0 * malpha[bs] * bs/(1.0+bs);
            } else {
                alpha2 = 1.0 * alpha * (bs)/(1.0+bs);
            }
            if (bs<=50) alpha2 = max(alpha2,lattice_tools::smallghconst[bs]);

            c[i] = 0;
            for (j=i-bs+1;j<i;j++) {
                c[i] += 1.0*log(c[j]) / (bs-1.0);
            }
            c[i] -= to_quad_float(1.0*log(VolumeBall(bs,1.0))/(bs-1.0));
            c[i] += bs/(bs-1.0)*log(alpha2);
            c[i] = exp(c[i]);
        }
       
        for (i=1;i<=n/2;i++) swap(c[i],c[n-i+1]);
        for (i=2;i<=n;i++) c[i] /= c[1];
        c[1] = 1.0;
       
        //modifying first exb indexes
        int exb = 1;    //index 1..i are modified
        quad_float* bitemp;
        RR cost;
        bitemp = new quad_float[n+1];
        while (beta + exb < n) {

            for (i=1;i<=n;i++) {
                bitemp[i] = c[i];
            }
            //generating "test" lengths of dim=exb
            //cout << "exb=" << exb << endl;
            for (i=exb;i>=1;i--) {
                int bb = beta + max((exb-i+1)/2,exb-2*(i-1));
                //bb = b;
                quad_float q;
                for (j=i+1;j<i+bb;j++) {
                    q += log(c[j]) / (bb-1.0);
                }
                q -= to_quad_float(1.0*log(VolumeBall(bb,1.0))/(bb-1.0));
                q += 1.0*bb/(bb-1.0) * log(alpha*bb/(1.0+bb));
                
                c[i] = to_double(exp(q));    //|b*_i|
            }
            cost = to_RR(ENUMCost(bitemp,1,exb+beta,prob,abar,'G',0,INPUT_NONSQUARED,1));
            if (bkz_verbose>=1) cout << "cost(" << exb << ")=" << cost << endl;
            if (maxcost < cost) {
                exb--;
                for (i=1;i<=n;i++) {
                    c[i] = bitemp[i];
                }
                break;
            }
            exb++;
        }

        //Adopting det=1
        double det=0;
        for (i=1;i<=n;i++) det += log(to_double(c[i]));
        det = exp(1.0/n * det);
        for (i=1;i<=n;i++) c[i] /= det;

   }    //end of BKZ20GSsim_second()
   
   void BKZGSsim(quad_float* c,int n,int beta) {

       sharememalloc();
       
        double startcputime = clock();
        double tptback = progressive_bkz::total_simulating_time; 
       
        std::ostringstream fname;
        fname << simcachedir << "/data" << n << "_" << beta << ".dat";
        if (FileExists(fname)==true) {
            //cout << fname.str() << endl;
            //load the data
            ifstream textlog; 
            textlog.open(fname.str().c_str(),ios::in);
            for (int i=1;i<=n;i++) {
                if (textlog.eof()==true) break;
                textlog >> c[i];
            }
            textlog.close();
            return;
        }

        BKZsim_first(c,n,beta);
        BKZsim_second(c,n,beta);
        ofstream textlog; 
        textlog.open(fname.str().c_str(),ios::trunc);
        for (int i=1;i<=n;i++) {
            textlog << c[i] << endl;
        }
        textlog.close();
        
        progressive_bkz::total_simulating_time = tptback + (clock()-startcputime) / CLOCKS_PER_SEC; 
   } // end of BKZ20GSsim()

    
    double ComputeLogTargetFEC(int n,int bs) {
        
        //Todo: will use cache files
        if (n<bs) return -1;
        
        sharememalloc();
        if (logfec[n][bs]==-1) {
            BKZGSsim(temp_c,n,bs);  //outputs GS length with det=1    
            logfec[n][bs] = to_double(log(lattice_tools::FullENUMCost(temp_c,n,INPUT_NONSQUARED)));
            //cout << "sim: logfec[" << n << "][" << bs << "]=" << logfec[n][bs] << endl;
            
        }
        return logfec[n][bs];
    }


    RR simulatetour(quad_float* c,int n,int beta,RR maxcost=to_RR(0)) {
        //returns cost
        //computes c[1...n] non squared
        //returns cost

        int i;
        lattice_tools::initialize();
        quad_float gh;
        RR largefactor;
        largefactor = 1;
        RR totalcost=to_RR(0);
        
        
        double r,abar,alpha,prob;
        r= compute_target_r(beta);
        abar = 1.0 / lattice_tools::gaussian_heuristic_constant[beta] * exp(-0.25*(beta-1.0)*log(r));
        alpha = abar * ((1.0+beta)/beta);
        if (alpha<1.00) alpha = 1.0;
        prob = 2.0 * exp(-1.0*beta*log(alpha));

        int exb=0;
        RR cost;
        quad_float* c2;
        c2 = new quad_float[n+1];
        for (i=0;i<=n;i++) c2[i] = c[i];        //backup

        //Detecting exb
        if (maxcost==0) {
            exb = 0;
        } else {
            if (beta<n) {
                do {
                    exb++;
                    cost = to_RR(ENUMCost(c,1,beta+exb,prob,abar,'G',0,INPUT_NONSQUARED,0));
                    //cout << "exb=" << exb << " cost=" << cost << "/" << maxcost << endl;
                } while ((cost < 2 * maxcost) && (beta+exb<n));
            } else {
                exb = 0;
            }
        }
        
        double modabar = abar;
        double modalpha = alpha;
        double modprob = prob;

        //cout << "alpha=" << alpha << " prob=" << prob << endl;
        
        int p02flag = 0;

        for (i=1;i<=n-min(50,beta);i++) {
            
            int bs = beta;  
            bs += max(0,max((exb-i+1)/2,exb-2*(i-1)));
            int b = min(bs,n-i+1);
            //cout << i << " " << b << endl;
            
            //gh
            if (b<=50) {
                gh = max(abar,lattice_tools::smallghconst[b]) * lattice_tools::LatticeGH(c+i-1,b,INPUT_NONSQUARED);
                //gh = LatticeGH(c+i-1,b,INPUT_NONSQUARED);
            } else {
                gh = abar * lattice_tools::LatticeGH(c+i-1,b,INPUT_NONSQUARED);
            }
            //cout << "gh1=" << gh << endl;
            //modify last indexes
            //cout << "maxcost=" << maxcost << endl;
            if (i>n-beta) {
                if (p02flag==0) {
                    modalpha *= 1.1; 
                    //if (modalpha > 1.0) modalpha = 1.05; 
                    do {
                        modalpha *= 0.99; 
                        modprob = 2.0 * exp(-1.0*b*log(modalpha));
                        cost = to_RR(ENUMCost(c2,n-b+1,b,modprob,modalpha,'G',0,INPUT_NONSQUARED,0));
                        //if (beta>80) cout << "i=" << i <<   " alpha=" << modalpha << " prob=" << modprob << " cost=" << cost << endl;
                    } while ((cost < maxcost) && (modprob<0.2));
                } else {
                    modprob = 0.2;
                    modalpha = exp(1.0/b*log(2.0/modprob));
                    cost = to_RR(ENUMCost(c2,n-b+1,b,modprob,modalpha,'G',0,INPUT_NONSQUARED,0));
                }
                if (modprob>0.2) p02flag  =1;
                
                modabar = modalpha / ((1.0+b)/b);
                //gh = modabar * LatticeGH(c+i-1,b,INPUT_NONSQUARED);
                if (b<50) {
                        gh = max(modabar,lattice_tools::smallghconst[b]) * lattice_tools::LatticeGH(c+i-1,b,INPUT_NONSQUARED);
                } else {
                        gh = max(modabar,1.0) * lattice_tools::LatticeGH(c+i-1,b,INPUT_NONSQUARED);
                }
                //cout << "lgh=" << lattice_tools::LatticeGH(c+i-1,b,INPUT_NONSQUARED) << endl;
                //cout << "gh2=" << gh << endl;   //Nan
                if (cost < 5*maxcost) totalcost += cost;
            } else {
                //cout << "b=" << b << endl;
                //cout << "prob=" << prob << endl;
                //cout << "alpha=" << alpha << endl;
                if (i==1) {
                    totalcost = 2 * maxcost;    //by definition
                } else {
                    //totalcost = 2 * maxcost;    //by definition
                   //totalcost += min(maxcost,to_RR(ENUMCost(c,i,b,prob,alpha,'G',0,INPUT_NONSQUARED,1)));
                   totalcost += to_RR(ENUMCost(c,i,b,prob,alpha,'G',0,INPUT_NONSQUARED,0));
                }
            }
            gh *= to_quad_float(exp(1.0 * log(largefactor) / b));
            if ((to_RR(c[i]) * largefactor) > to_RR(gh)) {
                largefactor *= to_RR(c[i] / gh);
                c[i] = gh;
            }
            //if (beta>140) cout << i << " " << totalcost << endl;
        }       
        
        //modify n-49...n
        quad_float det;
        for (i=n-min(50,beta)+1;i<=n;i++) det += log(c[i]) - log(lattice_tools::simhkzconst[i-n+50]);
        det += to_quad_float(log(largefactor));
        det = exp(1.0 * det / min(50,beta));
        for (i=n-min(50,beta)+1;i<=n;i++) c[i] = det * lattice_tools::simhkzconst[i-n+50];
/*
        for (i=1;i<=n;i++) {
            cout << "c[" << i << "]=" << c[i] << endl;
        }
 */
        //cout << "last=" << totalcost << endl;
        //exit(0);
        return totalcost;
    }
    
    RR TimeBKZ(int n,int sb,int eb,int beta,double A,double W,char model) {
        //n=dim, sb=start_beta, eb=end_beta, beta=applied beta
        //expected cost to generage BKZ-eb basis starting from BKZ-sb basis
        //using BKZ-beta algorithm

        std::ostringstream fname;
        if (model==modelQF) {
            mkdir( (simcachedir+"/timesim").c_str() ,0777); 
            fname << simcachedir << "/timesim/" << n << "_" << sb << "_" << eb << "_" << beta << ".txt";
        }
        if (model==modelRR) {
            mkdir( (simcachedir + "/timesim_modelRR").c_str(),0777);  //Todo: will be modified
            fname << simcachedir << "/timesim_modelRR/" << n << "_" << sb << "_" << eb << "_" << beta << "rr.txt";
        }

        if (FileExists(fname)==true) {
            RR* scost;
            scost = new RR[2];
            loadarray(scost,1,fname.str());
            if (scost[1]>0) {
                    //cout << "cost (BKZ-" << sb << "->" << eb<< ", usingBKZ-" << beta << ") in sec=" << scost[1] << endl;
                    //cout << "LastFEC=" << progressive_bkz::ComputeLogTargetFEC(n,eb) << endl;
                    return scost[1];
            }
        }

        quad_float* c;
        c = new quad_float[n+1];
        RR cost = to_RR(0);

        double tours;

        double targetfec;
        double currentfec;
        RR localcost,currentcost;;
        RR tcost;
        RR lllcost;
        char esflag = 0;
    #ifdef useeasysimulator
        double* esalpha,*talpha;
        talpha = new double[n+1];
        esalpha = new double[n+1];
    #endif
        if (sb==0) {
            simulateLLL(c,n);
        } else {
    #ifdef useeasysimulator
            if ((beta>=150) && (n>2*beta)) {
                esflag  = 1;
                pbkz_easysim_alpha(c,esalpha,n,beta); //esalpha[] is from beta
                pbkz_easysim_alpha(c,talpha,n,sb); //c[] is from sb
            } else {
                BKZGSsim(c,n,sb);  //outputs GS length with det=1    
            }
    #else
            BKZGSsim(c,n,sb);  //outputs GS length with det=1    
    #endif
        }


        //BKZ parameters
        double r,abar,alpha,prob;
        r= compute_target_r(beta);
        abar = 1.0 / lattice_tools::gaussian_heuristic_constant[beta] * exp(-0.25*(beta-1.0)*log(r));
        alpha = abar * ((1.0+beta)/beta);
        if (alpha<1.00) alpha = 1.0;
        prob = 2.0 * exp(-1.0*beta*log(alpha));


        RR maxcost = to_RR(0);
        if (esflag==0) {
            //Detecting exb in BKZ-beta
            for (int i=1;i<=n-beta;i++) {
                    maxcost = max(maxcost,to_RR(ENUMCost(c,i,beta,prob,alpha,'G',0,INPUT_NONSQUARED,1)));
            }
            //cout << "maxcost = " << maxcost << endl;
        } else {
            maxcost = exp(0.00229916935128517*beta*beta -0.217033716436143*beta + 16.2223328348192);
        }
        int tlim = 1000;

        //Simulate tours

        targetfec = to_double(progressive_bkz::ComputeLogTargetFEC(n,eb));
        //cout << "tfec=" << targetfec << endl;
        tours = 0;
        localcost = 0;
        RR simcost;
        double logfec;
        do {
            tours++;
            //backups
            currentfec = to_double(log(lattice_tools::FullENUMCost(c,n,INPUT_NONSQUARED)));
            //cout << "currentfec=" << currentfec << endl;
            currentcost = localcost;
            simcost = simulatetour(c,n,beta,maxcost);
            
            //cout << "simcost=" << simcost << endl;
            localcost += simcost;
            logfec = to_double(log(lattice_tools::FullENUMCost(c,n,INPUT_NONSQUARED)));
            //cout << "beta=" << beta  << "tour=" << tours << " tfec=" << logfec << "/" << targetfec << endl;

            if (abs(currentfec-logfec)<0.01) {
                //converged?
                tours=tlim+10;
                break;
            }
        } while ((logfec > targetfec) && (tours<tlim));

        tours--;        
        //imaginary last tours
        //cout << "imaginary #tours=" << tours << endl;
        if (tours<tlim-1) {
            tours += (currentfec-targetfec) / (currentfec - logfec);
            currentcost += simcost * (currentfec-targetfec) / (currentfec - logfec);
            //cout << "imaginary #tours=" << tours << endl;
        }
        if (tours<0) {
            cost = 1e+300;
        }
        if (tours>tlim) {
            cost = 1e+300;
        }

        tcost += currentcost * beta;
        if (model==modelQF) {
                lllcost += (double)tours * beta*beta*n*n*n;
        } else {
            double T;
            T = (n-beta)/(250.0-beta);
            //lllcost += (double)tours * T*(n-beta)*(n-beta);
            lllcost += (double)tours * T*250.0 * n*n;
        }

        if (tours<tlim) {
            cost = A * lllcost + W*tcost;
        } else {
            cost = 1e+300;
        }
        if (bkz_verbose>=1) cout << "cost (BKZ-" << sb << "->" << eb<< ", usingBKZ-" << beta << ") in sec=" << cost << endl;

        RR* scost;
        scost = new RR[2];
        scost[1] = cost;
        savearray(scost,1,fname.str());
        if (cost<0) cost=0;

        return cost;
    }

    RR TimeBKZmin(int n,int sb,int eb,double A,double W,char model,int& betaalg) {

        std::ostringstream fname;
        fname << "timesimB/" << n << "_" << sb << "_" << eb << ".txt";
        if (FileExists(fname)==true) {
            RR* scost;
            scost = new RR[3];

            savearray(scost,2,fname.str());
            conv(betaalg,scost[2]);
            return scost[1];
        }

        //n=dim, sb=start_beta, eb=end_beta
        //Detecting minimum cost to get BKZ-eb basis starting from BKZ-sb basis
        RR cost;
        RR mincost = to_RR(-1);
        int beta = eb;
        int minbeta = eb;

        if (eb<70) beta = eb+1; 

        while (1) {
            cost = TimeBKZ(n,sb,eb,beta,A,W,model);
            if (mincost==-1) {
                mincost = cost;
            }
            if (mincost >= cost) {
                mincost = cost;
                minbeta = beta;
            }
            beta++;
            if ((sb>0) &&  ((cost > mincost * 1.5) || (beta>eb*2.5) || (beta>=n))) break;
            if ((sb==0) &&  ((cost > mincost * 1.5) || (beta>eb*2.5) || (beta>=n))) break;
            if (cost<=0) {
                mincost = 0;    //error?
                break;
            }
        }

        RR* scost;
        scost = new RR[3];
        scost[1] = cost;
        scost[2] = minbeta;
        savearray(scost,2,fname.str());

        cout << "cost (BKZ-" << sb << "->" << eb<< ") in sec=" << mincost << " (using BKZ-" << minbeta << ")" << endl;
        betaalg = minbeta;
        return mincost;

    }
    
    void exp_svpchallenge_optbeta_multiple(int n,char model) {

        using namespace progressive_bkz;
        //compute the optimal strategy of beta
        double A1,W;
        if (model==modelQF) {
            A1 = 1e-10;
            W = 1.5e-8;
        }
        if (model==modelRR) {
            A1 = 1e-6;    
            W = 1.5e-8;
        }

        RR cost,mincost;
        int sb = 0;
        int i;
        int beta;
        int betaalg;
        int beta_s;

        //cost = TimeBKZmin(n,sb,eb,A1,W,model);

        //Computing Table of TimeBKZ(n,beta)
        RR* bkztime;
        bkztime = new RR[400];
        std::ostringstream* bkzstr;
        bkzstr = new std::ostringstream[400];
        std::ostringstream minstr;

        for (i=0;i<400;i++) bkztime[i] = to_RR(0);

        RR mintotal;
        mintotal = -1;

        for (beta=20;beta<n;beta++) {
            //computing minimum

            //starting from LLL
            if (beta<30) {
                mincost = TimeBKZmin(n,0,beta,A1,W,model,betaalg);
                bkzstr[beta] << bkzstr[sb].str() << "{\\scriptsize $\\PROCESSBKZX{ " << beta << "}{" << betaalg << "}$} & ";
            } else {
                mincost = to_RR(-1);
                bkzstr[beta].str("");
            }        

            beta_s = max(20,beta-20);
            if (beta>80) beta_s = beta-15;
            if (beta>100) beta_s = beta-5;

            for (int sb=beta_s;sb<beta;sb++) {
                cost = bkztime[sb] + TimeBKZmin(n,sb,beta,A1,W,model,betaalg);

                if ((mincost > cost) || (mincost==-1)) {
                    mincost = cost;
                    bkzstr[beta].str("");
                    //bkzstr[beta] << bkzstr[sb].str() << "->(" << beta << "," << betaalg << ")";
                    bkzstr[beta] << bkzstr[sb].str() << "{\\scriptsize $\\PROCESSBKZX{ " << beta << "}{" << betaalg << "}$} & ";
                }
            }
            bkztime[beta] = mincost;
            //cout << "BKZtime[" << n << "," << beta << "]=" << mincost << endl;
            //cout << "strategy[" << beta << "]=" << bkzstr[beta].str() << endl;

        }
    }
    
    int  detectBKZLevel(mat_ZZ& L) {

        int n = L.NumCols();
        double curlogfec = to_double(log(lattice_tools::FullENUMCost(L)));

        int sb;
        for (sb=2;sb<=n;sb++) {
            if (progressive_bkz::ComputeLogTargetFEC(n,sb) < curlogfec) break;
        }    
        return sb;
    }

    
}


#endif
