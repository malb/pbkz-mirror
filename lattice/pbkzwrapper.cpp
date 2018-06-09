#ifndef _inc_pbkzwrapper
#define _inc_pbkzwrapper

//Progressive BKZ main

#include "pbkzsimboost.cpp"
#include "bkzlibwrapper.hpp"
#include "enumwrapper.cpp"

template <typename T> inline void EraseZeroVectors(LatticeBasis<T>& B) {
    int j = 0;
    for (int i=0;i<B.L.NumRows();i++) {
        char nonzeroflag=0;
        int len = B.L[i].length();
        for (int k=0;k<len;k++) {
            if (B.L[i][k] != 0) {
                nonzeroflag = 1;
                break;
            }
        }
        if (nonzeroflag==1) {
            if (i!=j) swap(B.L[j],B.L[i]);
            j++;
        } else {
        }
    }
    B.L.SetDims(j,B.L.NumCols());
    
}


mat_ZZ BUCandidates;

template <typename T> int BasisUpdate(LatticeBasis<T>& B,mat_ZZ*U,int istart,int iend,FoundENUMVectors& EV,int vl,int option,char memoryflag=0) {
    //Insert vectors after B[istart]
    EV.duplicate();
    int nonzero = EV.countnonzero();
    if (nonzero==0) return 0;

    double ss = gettimeofday_sec();
    istart = EV.getminproj();
    
    int n = B.L.NumCols();
    int m = B.L.NumRows();
    
    if ((nonzero > BUBuffer.NumRows()) || (n > BUBuffer.NumCols())) {
       BUBuffer.SetDims(nonzero,n);
    }
    if ((nonzero > BUCandidates.NumRows()) || (n > BUCandidates.NumCols())) {
       BUCandidates.SetDims(nonzero,n);
    }

    if (B.L.NumRows() < iend+nonzero ) {
        B.L.SetDims(iend + nonzero,n);
    }


    //escape to buffer memory
    for (int i=0;i<nonzero;i++) {
        BUBuffer[i] = B.L[iend+i];
    }
    T* cback= (T*)shared_memory::allocate1<T>(750,iend+nonzero+1);
    for (int i=istart;i<=iend+nonzero;i++) cback[i] = B.gs.c[i];

    //make zero buffer
    for (int i=iend;i>=istart;i--) {
        swap(B.L[i+nonzero-1],B.L[i-1]);
    }

    //Here, B[istart...istart+nonzero-1] is free
    B.GScomputed = istart-1;
    
    //insert vectors (greedy)
    int j = 0;
    int numinsert=0;
    for (int i=0;i<EV.data.size();i++) {
        if (EV.data[i].norm > 0) {
            coeff_to_zz(BUCandidates[j++],EV.data[i].coeffs,B.L,EV.data[i].projection-1 + nonzero );
        }
    }
    int numcandidates = j;
    
    j = istart;
    T localdet = 1;
    T newdet=1;

    while (1) {
        double minproj = -1;
        int index;
        for (int i=0;i<numcandidates;i++) {
            double proj = ProjectedLength(B,1,j-1,BUCandidates[i]);
            if (proj>0.5) {
                if ((minproj == -1) || (proj < minproj)){
                    minproj = proj;
                    index = i;
                }
            }
        }

        if (minproj>0) {
            //the first inserting vector must improve the basis
            localdet *= cback[j];
            newdet *= minproj;
            if (newdet > localdet * 0.999) {
                break;
            }
            B.L[j-1] = BUCandidates[index];
            B.updateGSBasis(j,j);
            numinsert++;
            j++;
        } else {
            break;
        }
    }

    for (int i=j;i<istart+nonzero;i++) {
        clear(B.L[i-1]);
    }
    double delta = 0.99;
    int final_iend = local_LLL(B,0,delta,istart,iend+nonzero);
    debug_display( cout << "fiend=" << final_iend << " " << iend << endl;)
    
    if (final_iend != iend) {
        B.updateGSBasis();
        B.gs.displaymu();
        B.gs.displayc();
        exit(0);
    }


    B.GScomputed = final_iend;
    for (int i=0;i<nonzero;i++) {
        B.L[iend+i] = BUBuffer[i] ;
    }
    for (int i=0;i<nonzero;i++) {
        B.gs.c[iend+i+1] = cback[iend+i+1];
        B.gs.cd[iend+i+1] = B.gs.c[iend+i+1]*B.gs.c[iend+i+1];
    }
    //Todo: this re-computation of GS can be replaced by shifting 
    if (memoryflag==0) B.L.SetDims(m ,n);
    B.updateGSBasis(istart,iend);
    return numinsert;
}


#include "pbkz_smalllattice.cpp"
#include "pfuncwrapper.cpp"

template <typename PFLOAT> int EnumthenUpdate(LatticeBasis<PFLOAT>& B,mat_ZZ*U,int istart,int iend,double radius,bkzfloat prob,int pfoption,int vl=0,std::string stringoptions="") {
        //Note: after update, B.gs.c and B.gs.cd have correct values for 1..dim
        //      but mu after index iend+1 are in old values
    
        double ss = gettimeofday_sec(); //timer

        std::map<std::string,std::string> options;
        if (stringoptions!="") ExtractOptions(options,stringoptions);

        B.updateGSBasis(1,iend);        //Assume that GS-basis is computed 
        for (int j=istart;j<=iend;j++) SizeReduce(B,U,j,nogsupdate);

        int parallel = 1;
        if (options["parallel"]!="") {
            parallel = atoi(options["parallel"].c_str());
        }

        if ((vl>=1) && (pfoption!=0)) {
            cout << "set pruning function [" << istart << ":" << iend << "]    \r";
            cout.flush();
        }
        radius = min(radius,(double)B.gs.c[istart]);
        debug_display(cout << "radius=" << radius  << endl;);
        debug_display(B.gs.displayc(istart,iend););

        long int elim = 10000000;
        if (options["elim"]!="") {
            elim = atoi(options["elim"].c_str());
        }
        
        //elim = 1;   //temporal for debug

        PruningFunction PF;
        bkzfloat realprob = pruning_func::SetPruningFunction(B,PF, istart,iend,radius , prob, pfoption,vl);
        vec_ZZ zerovec;
        
        //PF.display();
        
        FoundENUMVectors EV;
        EV.clearmemory(16);
        
        int finishmode = finish_mode_nowaste;

        
        int enumoption = enum_mode_find_short_projection;
        int slide = iend-istart;
        
        int dim = iend-istart+1;
        
        if ( dim <= 40) {
            parallel = 1;
        }
        
        if ((vl>=1) && (pfoption!=0)) {
            cout << "process enum [" << istart << ":" << iend << "]    \r";
            cout.flush();
        }
        
        if (parallel==1) {
            if (dim <= 140) {
                lattice_enum::ENUMCoreboost_double(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
            } else {
                lattice_enum::ENUMCoreboost<long double,long double>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
            }
        } else {
            if (dim <= 140) {
                lattice_enum::ENUMCoreboostmt_double(EV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
            } else {
                lattice_enum::ENUMCoreboostmt<long double,long double>(EV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
            }
        }
        debug_display(cout << "enum_time=" << gettimeofday_sec() -ss << endl;)

        if ((vl>=1) && (pfoption!=0)) {
            cout << "basis update [" << istart << ":" << iend << "]    \r";
            cout.flush();
        }

        int ret=0;
        if (iend-istart+1 <= 60) {
            ret =  BasisUpdateGSsmall<double,PFLOAT>(B,U,istart,iend,EV,VL0,0,1);
        } else 
        if (iend-istart+1 <= 90) {
           ret =  BasisUpdateGSsmall<long double,PFLOAT>(B,U,istart,iend,EV,VL0,0,1);
        } else {
           ret =  BasisUpdateGSsmall<float15,PFLOAT>(B,U,istart,iend,EV,VL0,0,1);
        }
        debug_display(cout << "update_time=" << gettimeofday_sec() -ss << endl;)
        return ret;

}
#include "bkzpreprocessboost.cpp"
#include "pbkzsimtimeboost.cpp"

template <typename PFLOAT> void PrunedBKZ(LatticeBasis<PFLOAT>& B,mat_ZZ*U,int beta,double alpha,double prob,double tour,int vl=0,std::string stringoptions="") {


        std::map<std::string,std::string> options;
        ExtractOptions(options,stringoptions);

        int istart=1;
        int iend=B.dim;
        if (options["istart"]!="") {
            istart = atoi( options["istart"].c_str());
            if (istart < 1) istart = 1;
            if (istart > iend) istart = iend;
        }
        if (options["iend"]!="") {
            iend = atoi( options["iend"].c_str());
            if (iend > B.dim) iend = B.dim;
        }

        if (options["doLLL"]!="") {
            BigLLL(B.L,0,0.999,VL1);
        }

        B.updateGSBasis();
        
        
        int n = B.dim;
        double ss = gettimeofday_sec();

        int t = 0;
        while (1) {
            int totalinsert = 0;
            for (int i = istart;i<iend-1;i++) {
                if (vl>=2) cout << "time=" <<   gettimeofday_sec() - ss << " tour=" << t+1 << "/" << tour << " index=" << i << " |b1|=" << LengthOf(B.L[0]) << endl;
                int bstart = i;
                int bend = min(bstart+beta-1,iend);
                double radius = alpha * lattice_tools::LatticeGH(B.gs.c,bstart,bend);
                
                //Returned value is number of inserted vectors
                totalinsert += EnumthenUpdate<PFLOAT>(B,U,bstart,bend,radius,prob,0,vl,"parallel=1");
                for (int i=bstart;i<=bend;i++) {
                    SizeReduce(B,U,i);
                }
                int swap;
                local_LLL(B,0,0.999,bstart,bend,stnormal,VL0,swap);
            }
            if (vl>=1) cout << "time=" <<   gettimeofday_sec() - ss << " tour=" << t+1;
            if (tour>0) cout << "/" << tour;
            cout << " |b1|=" << LengthOf(B.L[0]) << endl;
            t++;
            if ((tour>0) && (t>=tour)) break;
            if (totalinsert==0) break;
        }
}


template <typename PFLOAT> int DetectBeta(LatticeBasis<PFLOAT>& B,int istart,int iend) {

    //Find PBKZ level: FEC(beta) > FEC(B) > FEC(beta+1)
    
    int dim = iend-istart+1;
    bkzfloat radius = lattice_tools::LatticeGH(B.gs.c,istart,iend);
    bkzfloat fec = FullENUMCost<bkzfloat,PFLOAT>(B.gs.c,istart,iend,(PFLOAT)radius);
    bkzfloat logfec = log(fec);
    for (int beta=10;beta<dim;beta++) {
        bkzfloat simfec = bkzconstants::simlogfec(dim,beta); 
        if (logfec > simfec) return beta-1; //level of log-fec between is sim[beta] and sim[beta-1] 
    }
    return dim;
}

//This subroutine is called by progressive BKZ routine
template <typename PFLOAT> int processoneblock(LatticeBasis<PFLOAT>& B,mat_ZZ*U,int beta,int bstart,int bend,double falpha,double ualpha,bkzfloat uprob,int parallel,
        int vl, double& totalpptime,double& totalpftime,double& totalenumtime,double& totalenumcputime ) {
    //preprocess
    int pfoption=0;
    int bs = bend-bstart+1; //blocksize 
    int updateflag=0;

    bool enumflag=true;
    
    if ((beta >= 50) && (bs >= 30)) {
        int Ca = 4;
        bkzfloat costlimit_lb = bkzconstants::target_cost_lb(beta) * Ca;
        bkzfloat costlb;
        double lss = gettimeofday_sec(); //timer
        if (bend-bstart<=80) {
            costlb = ProgressiveBKZPreprocess<PFLOAT,double>(B,U,bstart,bend,costlimit_lb,falpha,ualpha,uprob,vl-1);
        } else {
            costlb = ProgressiveBKZPreprocess<PFLOAT,long double>(B,U,bstart,bend,costlimit_lb,falpha,ualpha,uprob,vl-1);
        }
        totalpptime += (gettimeofday_sec()-lss);
        if (costlb < 1e+6) parallel=1;  //to prevent slowdown from overhead
        
        if (costlb > costlimit_lb) enumflag = false;
        
    }
    if (bend-bstart+1<=50) parallel=1;

    //main process
    if (vl>=2) cout << "ghcheck:" << B.gs.c[bstart] << " " << lattice_tools::LatticeGH(B.gs.c,bstart,bend)  << " " << B.gs.c[bstart]/lattice_tools::LatticeGH(B.gs.c,bstart,bend) << endl;
    if (B.gs.c[bstart] < falpha * lattice_tools::LatticeGH(B.gs.c,bstart,bend)) enumflag = false; 
    
    
    if (enumflag==true) {
        double lss = gettimeofday_sec(); //timer
        auto gh = lattice_tools::LatticeGH(B.gs.c,bstart,bend);
        double radius = ualpha * (double)gh;
        auto prev = B.gs.c[bstart];
        updateflag += EnumthenUpdate<PFLOAT>(B,U,bstart,bend,radius,uprob,pfoption,vl,"parallel=" + to_stdstring(parallel));
        if (vl>=1) {
            cout << "enum(" << bstart << "," << bend << ")                       \r";
        }

        totalenumtime += lattice_enum::current_etime;
        totalenumcputime += lattice_enum::current_cputime;

        //Post process
        for (int i=bstart;i<=bend;i++) {
            SizeReduce(B,U,i);
        }
    }
    return updateflag;
}

template <typename PFLOAT> void UpdateSubWindow(mat_ZZ& WU,LatticeBasis<PFLOAT>& B,mat_ZZ* U,int wstart) {
    
    int dim = WU.NumRows();

    mat_ZZ Btemp;
    Btemp.SetDims(dim,B.L.NumCols());

    for (int i=0;i<dim;i++) {
        for (int j=0;j<dim;j++) {
            Btemp[i] += WU[i][j] * B.L[j+wstart-1];
        }
    }
    for (int i=0;i<dim;i++) {
        B.L[i+wstart-1] = Btemp[i];
    }
    B.GScomputed = min(B.GScomputed,wstart-1);
    B.updateGSBasis(1,wstart+dim-1);;

    if (U!=0) {
        Btemp.SetDims(dim,(*U).NumCols());
        for (int i=0;i<dim;i++) {
            clear(Btemp[i]);
        }
        for (int i=0;i<dim;i++) {
            for (int j=0;j<dim;j++) {
                Btemp[i] += WU[i][j] * (*U)[j+wstart-1];
            }
        }
        for (int i=0;i<dim;i++) {
            (*U)[i+wstart-1] = Btemp[i];
        }
    }
}

template <typename PFLOAT> void StrategiedBKZ(LatticeBasis<PFLOAT>& B,mat_ZZ*U,BKZStrategy& BS,int vl=0,std::string stringoptions="") {

    B.updateGSBasis();
    bkzconstants::loadcachetable();
    
    //ExtractOptions
    std::map<std::string,std::string> options;
    ExtractOptions(options,stringoptions);
    bool flag_betashiftupdate = true;
    if (options["betashiftupdate"]=="false") {
        flag_betashiftupdate = false;
    }

    bool flag_ignoreflat = false;
    if (options["ignoreflat"]=="true") {
        flag_ignoreflat = true;
    }
    
    int endflag=0;
    int beta = 10;
    int iistart = 1;
    int iiend = B.dim;
    if (options["istart"]!="") {
        iistart = atoi( options["istart"].c_str());
        if (iistart < 1) iistart = 1;
        if (iistart > iiend) iistart = iiend;
    }
    if (options["iend"]!="") {
        iiend = atoi( options["iend"].c_str());
        if (iiend > B.dim) iiend = B.dim;
    }

    std::ofstream of;
    if (options["logfile"]!="") {
        of.open(options["logfile"],ios::app);
    }


    int parallel = 1;
    if (options["parallel"]!="") {
        parallel = atoi( options["parallel"].c_str());
        if (parallel < 1) parallel = 1;
    }

    PFLOAT target_length = 0;
    if (options["targetlen"]!="") {
        target_length = boost::lexical_cast<PFLOAT>(options["targetlen"]);
    }
    
    char ignoreflat = false;
    if (options["ignoreflat"]!="") {
        ignoreflat = true;
    }

    double start = gettimeofday_sec();
    double ssim = 0;
    int tour = 1;
    B.updateGSBasis();
    PFLOAT gh = lattice_tools::LatticeGH(B.gs.c,iistart,iiend);
    double totalenumtime = 0;
    double totalenumcputime = 0;
    double totalpptime = 0; //Time for BKZ-preprocess
    double totalpftime = 0;
    totalsimulatingtime=0;
    
    if (gh<=0) {
        cout << "gh error: " << gh << endl;
        B.gs.displayc(iistart,iiend);
        exit(0);
    }
    
    bkzfloat prevfec =  log(FullENUMCost<bkzfloat,PFLOAT>(B.gs.c,iistart,iiend,(PFLOAT)gh));
    if (prevfec<0) {
        cout << "fec error: " << prevfec << endl;
        B.gs.displayc(iistart,iiend);
        exit(0);
    }
    
    int achievedbeta = 0;
    int bsline = 0;
    int betashift = 0;
    
    while (endflag==0) {
        int istart = iistart;
        int iend = iiend;
        if (flag_ignoreflat==true) {
            //narrow the processing range
            while ((istart < iend-1) && ( abs(B.gs.c[istart] -  B.gs.c[istart+1])/B.gs.c[istart] < 0.001)  ) istart++;
            while ((istart < iend-1) && (  abs(B.gs.c[iend] -  B.gs.c[iend-1])/B.gs.c[iend] < 0.001) )  iend--;
            istart = max(iistart,istart-5);
            iend = min(iiend,iend+5);
            if ((istart!=iistart) || (iend!=iiend)) {
                if (vl>=1) cout << "reset processing range[" << iistart << "," << iiend << "] -> [" << istart << "," << iend << "]                       " << endl;
            }
            gh = lattice_tools::LatticeGH(B.gs.c,istart,iend);
        }

        //Detect blocksize
        int betadetect = DetectBeta(B,istart,iend);
        if (betadetect>=iend-istart+1) return;    //already reduced

        achievedbeta = max(achievedbeta,betadetect);

        while ((bsline < BS.process.size()) && (achievedbeta >= BS.process[bsline].ebeta)) bsline++;
        if (bsline == BS.process.size()) {
            endflag = 1;
            break;
        }
        beta = min(BS.process[bsline].usebeta , iend-istart+1);
        bkzconstants::simlogfec(iend-istart+1,beta); //call to make (alpha,prob) table
        
        
        double alpha,r;
        bkzfloat prob;
        PBKZParam(alpha,r,prob,beta);

        int updateflag=0;

        bkzfloat fec = FullENUMCost<bkzfloat,PFLOAT>(B.gs.c,istart,iend,(PFLOAT)gh);
        bkzfloat logfec = log(fec);
        bkzfloat simfec = bkzconstants::simlogfec(iend-istart+1,achievedbeta + 1);    //next target fec 
        bkzconstants::savecachetable();

        if ((flag_betashiftupdate==true) && (flag_ignoreflat==false))  {
            if (logfec > prevfec) {
                betashift++;
                beta = min(BS.process[bsline].usebeta + betashift , iend-istart+1);
            }
        }
        prevfec = logfec;

        if (vl>=1) {
            cout << "time=" <<   gettimeofday_sec() - start - (totalsimulatingtime - ssim) <<  " tour=" << tour << " applybeta=" << beta; 
            cout << " basisbeta=" <<  betadetect << " |b1|=" << B.gs.c[istart];
            cout << "=" << B.gs.c[istart] / gh << "GH";
            cout << " logfec=" << logfec << "/" << simfec;
            cout << endl;
        }
        if (of.is_open()==true) {
            of << gettimeofday_sec() - start - (totalsimulatingtime - ssim)  << "\t";
            of << tour << "\t";
            of << beta << "\t";
            of << betadetect << "\t";
            of << B.gs.c[istart] << "\t";
            of << B.gs.c[istart] / gh<< "\t";
            of << logfec << "\t";
            of << simfec << "\t";
            of << totalenumtime << "\t";
            of << totalenumcputime << "\t";
            of << totalpptime << "\t";
            of << totalsimulatingtime - ssim << endl;
        }

        if (endflag==1) break;
        tour++;

        std::vector<double> falphatable,ualphatable;
        std::vector<bkzfloat> probtable;
        SimulateAlphatable(falphatable,ualphatable,iend-istart+1,beta);
        SimulateProbtable(probtable,iend-istart+1,beta);

        //start tour
        int wstart,wend;    //window_start and window_end
        wstart = -1;
        wend = -1;
        int wmaxsize = 2*beta; //max window size (heuristic)
        int wsize;

        LatticeBasis<long double> WB;   //window basis
        mat_ZZ WU;  //window unitary

        bool usewindow=true;
        
        for (int i = istart;i<iend-2;i++) {

            int bstart = i;
            int bend = min(bstart+beta-1,iend);
            int bs = bend-bstart+1; //blocksize 
            double falpha =  falphatable[i-istart+1];    //alpha expected to find
            double ualpha = ualphatable[i-istart+1];; //alpha used to search
            bkzfloat uprob = probtable[i-istart+1];
            if ((vl>=1) && (beta>=45)) {
                cout << "process [" << bstart << "," << bend <<  "]  \r";
                cout.flush();
            }            

            if (usewindow==false) {
                //Not use window
                B.updateGSBasis(istart,bend);
                updateflag += processoneblock(B,U,beta,bstart,bend,falpha,ualpha,uprob,parallel,vl,totalpptime,totalpftime,totalenumtime,totalenumcputime);
                int swapcount;
                local_LLL(B,U,0.999,1,bend,stnormal,min(vl,VL0),swapcount);
            } else {
                if (bend > wend) {
                    if (wend>0) {
                        if (vl>=2) cout << "update window[" << wstart << ":" << wend << "]                              \r" << endl;
                        UpdateSubWindow(WU,B,U,wstart); //update previous window
                        //Blockwised LLL started at i=wstart
                        int llsize = 160;
                        int lstart = wstart;
                        int lend;
                        int localswap,maxshift;
                        while (1) {
                            lend = min(lstart+llsize -1,wend);
                            SmallLocalLLL<long double>(B,U,(long double)0.999,lstart,lend,VL0,localswap,maxshift);
                            for (int i=maxshift;i<=lend;i++) {
                                SizeReduce(B,U,i);
                            }
                            if (localswap > 0) {
                                if ((lstart==istart) && (lend==wend)) break;
                                if ((lstart < maxshift) && (lend==wend)) break;
                                lstart = maxshift - llsize/2;
                                if (lstart < istart) lstart = istart;
                            } else {
                                if (lend == wend) break;
                                lstart += llsize;
                            }
                            lend = min(lstart+llsize -1,wend);
                        }
                    }    

                    //define new window
                    wstart = max(istart,bstart-10);
                    wend = min(iend,bstart + wmaxsize-1);
                    wsize = wend-wstart+1;
                    if (vl>=2) cout << "new window[" << wstart << ":" << wend << "]" << endl;

                    if ((tour!=1) || (wstart!=istart)) {
                        //if else, size reduced by the previous loop
                        for (int i=wstart;i<=wend;i++) {
                            cout << "size reduce:" << i << "/" << wend << "                    \r";
                            cout.flush();
                            SizeReduce(B,U,i);
                        }
                    }

                    WB.L.SetDims(wsize,wsize);
                    WU.SetDims(wsize,wsize);
                    ExtractSubBlock(WB,B,wstart,wend);  //force min|b*i|=10000
                    SetIdentity(WU);
                    WB.GScomputed=-1;
                    WB.updateGSBasis();
                    if (WB.dim>=140) {
                        BigLLL(WB,&WU,0.999,1,WB.dim,min(vl,1));
                    } else {
                        int swapcount;
                        local_LLL(WB,&WU,0.999,1,WB.dim,stnormal,min(vl,VL1),swapcount);
                    }
                }
                updateflag += processoneblock(WB,&WU,beta,bstart-wstart+1,bend-wstart+1,falpha,ualpha,uprob,parallel,vl,totalpptime,totalpftime,totalenumtime,totalenumcputime);
                int swapcount;
                local_LLL(WB,&WU,0.999,1,bend-wstart+1,stnormal,min(vl,VL0),swapcount);
            }
            //display
            if (vl>=2) cout << "time=" <<   gettimeofday_sec() - start - (totalsimulatingtime - ssim)<<  " index=" << i << " beta=" << beta << " alpha=" << ualpha << " prob=" << uprob << endl;
            bkzconstants::savecachetable();
        }       //end of tour
        if (wend>0)  {
            if (vl>=2) cout << "update window[" << wstart << ":" << wend << "]" << endl;
            UpdateSubWindow(WU,B,U,wstart); //update previous window
            for (int i=wstart;i<=wend;i++) {
                if (vl>=1) cout << "size reduce:" << i << "/" << wend << "                    \r";
                cout.flush();
                SizeReduce(B,U,i);
            }
        }
        if (updateflag==0) beta++;
        if (vl>=2) cout << "time=" <<  gettimeofday_sec() - start - (totalsimulatingtime - ssim) <<  " tour=" << tour << " finished" << endl;
        if (options["temporal"]!="") {
            SaveLattice(B,options["temporal"]);
        }

        if (target_length>0) {
            if (B.gs.c[iistart] < target_length) {
                break;
            }
        }
    }
    for (int i=iistart;i<=iiend;i++) {
        cout << "LastSizeReduce i=" << i << "              \r";
        cout.flush();
        SizeReduce(B,U,i);
    }
    if (vl>=1) cout << "progressive bkz finish. time=" <<   gettimeofday_sec() - start - (totalsimulatingtime - ssim) << "                          " << endl;
    bkzconstants::savecachetable();
}

template <typename PFLOAT> void ProgressiveBKZ(LatticeBasis<PFLOAT>& B,mat_ZZ*U,int targetlevel,int vl=0,std::string stringoptions="") {

    //ExtractOptions
    std::map<std::string,std::string> options;
    ExtractOptions(options,stringoptions);

    int betashift = 6;
    if (options["betashift"]!="") {
        betashift = atoi( options["betashift"].c_str());
        if (betashift < 1) betashift = 1;
    }

    int startbeta=15;
    if (options["startbeta"]!="") {
        startbeta = atoi( options["startbeta"].c_str());
        if (startbeta < 1) startbeta = 1;
    }

    BKZStrategy BS;
    for (int i=startbeta-betashift+1;i<=targetlevel;i++) {
        BS.addstrategy(i-1,i,i+betashift-1,0,0,0);  //Last 0's are dummy
    }
    StrategiedBKZ(B,U,BS,vl,stringoptions);    
}

#endif

