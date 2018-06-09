#ifndef _inc_latticeenum_close
#define _inc_latticeenum_close

namespace lattice_enum {


    void ENUMCVPCoremt(int** foundvec,double* foundlen,int fnum,int& findex,int& slide,int jj,int bs,long int elim,double clim,double* normal_pruning_func,quad_float* target_mu,quad_float** mu,quad_float* c,int vl,int numthreads,char enummode,char finishmode) {

        //Find a close vector to the target point

        //stop after elim*10^8 nodes are processed
        //ENUM0 = enumeration for L(b0,...,bk)
        //ENUM1a = enumeration for L(b[k+1],...,b[n-1])
        //ENUM1b = sub-enumeration for vectors from ENUM1a

        initialize();
        int i,j;

        //memories for enumeration {\phi,q,so}{ctilda,vvec,yvec,uvec,Deltavec,deltavec} are defined in ENUMsharememory

        //memory for ENUM0
        double *ctilda = enum_share_double[0];
        double *vvec = enum_share_double[1];
        double *yvec = enum_share_double[2];
        double *uvec = enum_share_double[3];
        double *utildavec = enum_share_double[4];
        long *Deltavec = enum_share_long[0];
        long *deltavec = enum_share_long[1];

        //memory for ENUM1a
        double *qctilda = enum_share_double[10];
        double *qvvec = enum_share_double[11];
        double *qyvec = enum_share_double[12];
        double *quvec = enum_share_double[13];
        double *qutildavec = enum_share_double[14];
        long *qDeltavec = enum_share_long[4];
        long *qdeltavec = enum_share_long[5];
        int sbk,tbk;    //backup of s and t in ENUM1a

        //memory for ENUM1b
        for (i=0;i<slotmax;i++) sos[i] = -1;   //default: the slot is not used

        double es;
        es = gettimeofday_sec();   //timer


        int kk = jj+bs-1;
        slide = bs;
        if (enummode==enum_mode_all_vectors) slide=0;


        double* pruning_func = enum_share_double[7];
        double* pftemp = pruning_func - jj;
        for (i=0;i<bs;i++) {
            pruning_func[i] = normal_pruning_func[i] * clim;
            //cout << "pf[" << i << "]=" << pruning_func[i] << endl;
        }
    /*
        if ((enummode==enum_mode_find_shortest) || (enummode==enum_mode_find_abort)) {
            for (i=0;i<0.25*(kk-jj);i++) pruning_func[i] = max(pruning_func[i],to_double(c[jj+i])); 
        }
     */ 
        //initialize memory for vector holding
        for (j=0;j<fnum;j++) {
            foundlen[j] = -1;
            foundvec[j][0] = jj + bs + 1;
            for (i=jj;i<=kk;i++) foundvec[j][i] = 0; 
        }

        //transmu
        double** transmu = enum_share_transmu;
        double* cd = enum_share_double[5];
        double* target_mud = enum_share_double[6];
        for (i=0;i<bs;i++) {
            for (j=i;j<bs;j++) {
                transmu[i][j] = to_double(mu[j+jj][i+jj]);
            }
        }
        //cd[i]=|b*_i|
        for (i=0;i<bs;i++) cd[i+jj] = sqrt(to_double(c[i+jj]));
        for (i=0;i<bs;i++) target_mud[i+jj] = to_double(target_mu[i+jj]);

        //initial values for ENUM0
        utildavec[jj] = uvec[jj] = 1;
        yvec[jj] = vvec[jj] = 0;
        Deltavec[jj] = 0;
        deltavec[jj] = 1;   
        for (i = jj+1; i <= kk+1; i++) {
            ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
            Deltavec[i] = 0;
            vvec[i] = 0;
            deltavec[i] = 1;
        }

        //initial values for ENUM1a
        int jjs;    //depth of ENUM1a

        //automatic search for jjs
        double cost;

        //memo: pruning func[0..bs-1]

        for (i=0;i<bs/2;i++) swap(normal_pruning_func[i] , normal_pruning_func[bs-1-i]);

        cost =  to_double(pruning_func::Rigid_upper_cost(normal_pruning_func-1,c+jj-1,bs,sqrt(clim)));
        if (vl>=3) cout << "upper # points in cylinder-intersection=" << cost << endl;
        cost *= 2.0;
        if (vl>=3) cout << "upper # processed nodes=" << cost << endl;
        double* costvec = enum_share_double[8];  //temporary used
        for (i=2;i<=min(50,bs-2);i+=2) {
            costvec[i] = to_double(Rigid_upper_cost_enuma1(normal_pruning_func-1,c+jj-1+(bs-i),i,sqrt(clim)));
        }
        for (i=3;i<=19;i+=2) costvec[i] = sqrt(costvec[i-1]*costvec[i+1]);  

        i=2;
        while (i<=50) {
            if (costvec[i]>1000) {
                jjs = i;
                break;
            }
            if (cost/costvec[i]<1e+6) { 
                jjs = i;
                break;
            }
            i++;
        }
        if (vl>=3) {
            cout << "expected # subtree (depth=" << i << "): " << costvec[i] << endl;
            //cout << "expected average # nodes in each subtree (depth=" << i << "): " << cost/costvec[i] << endl;
            //cout << "setting depth of tree in ENUM1a = " << jjs << endl;
        }
        for (i=0;i<bs/2;i++) swap(normal_pruning_func[i] , normal_pruning_func[bs-1-i]);
        jjs = kk-(jjs-1);


        //computing bias=nearest_plane(L,t)
        //this algorithm computes vectors close to t-bias
        double* bias;
        bias = new double[kk+1];
        //cout << "bias=";
        for (i=kk;i>=jj;i--) {
             bias[i] = round(to_double(target_mud[i]));
             //cout << bias[i] << " ";
             target_mud[i] -= bias[i];
             for (j=jj;j<i;j++) {
                 target_mud[j] -= bias[i] * to_double(mu[i][j]);
             }
             //cout << i << " : " << target_mud[i] << " ";
         }
        //cout << endl;

        qutildavec[jjs] =1;
        quvec[jjs] = 1;
        qyvec[jjs] = qvvec[jjs] = 0;
        qDeltavec[jjs] = 0;
        sbk = tbk = jjs;
        qdeltavec[jjs] = 1;   
        for (i = jj+1; i <= kk+1; i++) {
            qctilda[i] = quvec[i] = qutildavec[i] = qyvec[i] = 0;
            qDeltavec[i] = 0;
            qvvec[i] = 0;
            qdeltavec[i] = 1;
        }

        //initial_value=a distance when a[i]=0 (projective norm of target point)
        qctilda[kk] = target_mud[kk]*target_mud[kk]*cd[kk]* cd[kk];
        for (i=kk-1;i>=jj;i--) {
              qctilda[i] = qctilda[i+1] + target_mud[i]*target_mud[i] * cd[i]* cd[i];
              qutildavec[i] = round(target_mud[i]);
        }

        //initial_value=a distance when a[i]=0 (projective norm of target point)
        ctilda[kk] = target_mud[kk]*target_mud[kk]*cd[kk]* cd[kk];
        for (i=kk-1;i>=jj;i--) {
              ctilda[i] = ctilda[i+1] + target_mud[i]*target_mud[i] * cd[i]* cd[i];
              utildavec[i] = round(target_mud[i]);
        }

        unsigned int lmct;  //local enum count
        long long int gmct;   //global enum count
        gmct = 0;


        int ea1thread=-1;  //thread number working for 

        //flags used in subroutines
        int endflag=0;  // set endflag=kk to abort all enums
        char ea1flag=0;   //ea1flag=0 <=> No thread works for ENUM1a
                          //       =1 <=>  a thread works
                          //       =2 <=> ENUM1a is finished
        char ea1breakflag=0;    //sufficient number of vectors may be slotted in ENUM1a



        //start of parallel-ENUM
        int vcount=0;   //# processed nodes in ENUM-1a

        #pragma omp parallel  num_threads(numthreads) private(i,j,lmct) 
        {

            //local variables (not shared)
            int current_used_slot;  //# used slot 
            int vindex; //for searching empty slot in ENUM1a
            char cflag;
            double t1d;
            double* ttmu;
            double* tut;
            double dtemp;
            int s,t;


            double mlen;
            int mythread = omp_get_thread_num();
            lmct = 0;  


            if (mythread==0) {
                //Thread0 processes ENUM0

                if (vl>=3) {
                    cout << "Start ENUM-0 (initial radius=" << sqrt(clim) << ")" << endl;
                }

                s=t=jj;

                //the below is the command to skip the ENUM-0
                //t = jjs;

                while (t <= jjs-1) {
                    dtemp = (yvec[t]+utildavec[t]- target_mud[t]) * cd[t];
                    dtemp *= dtemp;
                    ctilda[t] = ctilda[t+1] + dtemp;
                    cflag = 0;
/*
                    for (i=jj;i<=s;i++) {
                        cout << utildavec[i] << " ";
                        if (i==t) cout << " * ";
                    }
                    cout << ctilda[t+1] << " " << ctilda[t] << " " << pftemp[t] << " " << s << endl;
  */  
                    if (ctilda[t] <= pftemp[t]) {
                        if (t<=jj+slide) {
                            //a vector s.t. pi_{i}|v|<|b*_i| (i=jj+slide+4) is found
    /*
                            cout << "t1:" << ctilda[t] << " "<< pftemp[t] << endl;
                            cout<< "jj" << jj<<  " slide=" << slide<< endl;
                            for (i=jj;i<=kk;i++) {
                                cout << ctilda[i] << " ";
                            }
                            cout << endl;
      */                      

                            if ((enummode & enum_mode_find_shortest) || (enummode & enum_mode_find_abort)) {
                                    //Store the found vector to local memory 
                                    mlen = storetomemory(foundlen, foundvec, fnum, utildavec, ctilda[t], t, kk,yvec,mu,jj,deltavec);
                                    if (mlen >= 0) {
                                        if (t-jj < slide) {
                                            slide = t-jj;
                                        }
                                    }
                                    if ((mlen != 0 ) && (mlen < clim) && (slide==0)) {
                                        //update pruning function if the pruning radius is updated
                                        #pragma omp critical
                                        {
                                            clim = mlen;
                                            for (i=0;i<kk-jj+1;i++) pruning_func[i] = normal_pruning_func[i] * clim;
                                            //for (i=0;i<0.25*(kk-jj);i++) pruning_func[i] = max(pruning_func[i],to_double(c[jj+i]));
                                            if (enummode==enum_mode_find_abort) {
                                                endflag = kk;
                                            }
                                        }
                                    }
                               // }
                            }
                            if (t==jj) {
                                //cout << t << " " << ctilda[t] << " " << pftemp[t] << endl;
                                if (enummode & enum_mode_all_vectors) {
                                //holding
                                    storetomemory_allmode(foundlen,foundvec,fnum,findex,utildavec,ctilda[t],jj,kk);
                                    //store v + k*b_1 with small norm
                                    double rr = yvec[t]+utildavec[t] - target_mud[t];
                                    double sign = +1;
                                    double utback = utildavec[t];
                                    if ((rr>0) && (s>jj)) sign=-1;
                                    i=2;
                                    while (1) {
                                        dtemp = yvec[t]+utback - target_mud[t]+ sign*(int)(i/2);
                                        dtemp *= cd[t];
                                        dtemp *= dtemp;
                                        ctilda[t] = ctilda[t+1]+dtemp;
                                        if (ctilda[t]< clim) {
                                            utildavec[t] = utback +  sign*(int)(i/2);
                                            storetomemory_allmode(foundlen,foundvec,fnum,findex,utildavec,ctilda[t],jj,kk);
                                        } else {
                                            break;
                                        }
                                        if (s>jj) {
                                            i++;
                                            sign=-sign;
                                        } else {
                                            i += 2;
                                        }
                                        if (findex==fnum) break;
                                    }
                                    if (findex==fnum) {
                                        endflag = kk;
                                        break;    //buffer memory is full
                                    }                                    
                                }
                            }
                        }

                        if (t > jj) {
                            cflag = 1;
                            t--;
                            t1d = 0;
                            ttmu = transmu[t-jj];
                            tut = utildavec + t + 1;
                            ttmu += t+1-jj;
                            i = s - t + 1;

    #ifdef AVX_PRD
                            __m256d ymm;
                            double atmp[4] __attribute__((aligned(32)));
                            i = s - t;
                            ymm = _mm256_setzero_pd();
                            while(i & 3){
                                t1d += (*tut) * (*ttmu);
                                ++tut, ++ttmu, --i;
                            }
                            while(i > 0){
                                ymm = _mm256_add_pd(ymm, _mm256_mul_pd(_mm256_loadu_pd(tut), _mm256_loadu_pd(ttmu)));
                            // ymm = _mm256_fmadd_pd(_mm256_loadu_pd(tut), _mm256_loadu_pd(ttmu), ymm);
                                tut += 4, ttmu += 4, i -= 4;
                            }
                            _mm256_store_pd(atmp, ymm);
                            for(int j = 0; j < 4; ++j) t1d += atmp[j];
    #else
                            while (--i) {
                                t1d += (*tut)*(*ttmu);
                                ttmu++;
                                tut++;
                            }
    #endif

                            yvec[t] = t1d;
    #ifdef SSE_ROUND
                            //ASM part is not fixed
                            __m128d tmp;
                            tmp = _mm_load_sd(&t1d);
                            tmp = _mm_round_sd(tmp, tmp, 0);
                            _mm_store_sd(&t1d, tmp);
                            t1d = -t1d;
    #else
                            //t1d = -round(t1d);
                            t1d = -round(t1d- target_mud[t]);
    #endif

                            utildavec[t] = vvec[t] = t1d;
                            Deltavec[t] = 0;
                            if (utildavec[t] > -yvec[t]) { 
                                deltavec[t] = -1;
                            } else {
                                deltavec[t] = 1;
                            }
                        } else {
                        }
                    }



                    if (cflag==0) {
                        t++;
                        if (s>t) {
                            Deltavec[t] = -Deltavec[t];
                        } else {
                            s = t;
                        }
                        Deltavec[t] += (1+((Deltavec[t]*deltavec[t]) >> 31))*deltavec[t];
                        utildavec[t] = vvec[t] + Deltavec[t];
                    }

                    //counter
                    if ((++lmct)^100000000) {
                    } else {
                        #pragma omp critical
                        {
                            gmct += lmct / 1000000;
                            lmct %= 1000000;
                            if (vl>=3) {
                                cout << "ENUM-0 count=" << gmct*0.01 << "*10^8 ";
                                cout << " clim=" << sqrt(clim) ;
                                cout << " time="  << gettimeofday_sec() - es;
                                cout << " speed=" << 1.0*gmct / (gettimeofday_sec() - es) << "Mnodes/sec ea1=" << (int)ea1flag << "  \r";
                                cout.flush();
                            }
                            if (gmct /100 > elim) {
                                endflag = kk;
                            }
                        }
                        t += endflag;
                    }


                } // end of ENUM0-core
                #pragma omp critical
                {
                    gmct += lmct / 1000000;
                    lmct %= 1000000;
                }
                if (vl>=3) {
                    cout << "ENUM-0 END count=" << gmct*0.01 << "*10^8 ";
                    cout << " clim=" << sqrt(clim) ;
                    cout << " time="  << gettimeofday_sec() - es;
                    cout << " speed=" << 1.0*gmct / (gettimeofday_sec() - es) << "Mnodes/sec ea1=" << (int)ea1flag << "  \r";
                    cout.flush();
                }
            }   // end of ENUM0


            //start of ENUM1a
            do {     
                if (ea1flag==0) {
                    #pragma omp critical
                    {
                        if (ea1thread==-1) {
                            ea1thread = mythread; 
                        }
                    }
                }            
                if ((mythread==ea1thread) && (ea1flag==0)) {
                    current_used_slot=0;
                    for (i=0;i<slotmax;i++) {
                        if (sos[i]!=-1) current_used_slot++;
                    }
                    if (current_used_slot > 0.9 * slotmax) {
                        ea1flag=0;  //no thread works for ENUM-1a
                    } else {
                        //Start of ENUM1a
                        //Searching vector candidates while current_used_slot < 0.9*slotmax
                        t = tbk;
                        s = sbk;
                        ea1flag = 1;

                        vindex=0;
                        //Start of ENUM1a-core
                        while (t <= kk) {
                            dtemp = (qyvec[t]+qutildavec[t]- target_mud[t]) * cd[t];
                            dtemp *= dtemp;
                            qctilda[t] = qctilda[t+1] + dtemp;
                            cflag = 0;
     /*                       
                            int abssum = 0;
                            for (i=t;i<=s;i++) abssum += qutildavec[i]*qutildavec[i];
                            //cout << s << " " << t << " " << abssum << endl;
                            if ((qctilda[t] <=  pftemp[t]) && (abssum<=4)) {
        */                    
                            if (qctilda[t] <=  pftemp[t]) {
                                if (t >= jjs) {
                                    cflag = 1;
                                    t--;
                                    t1d = 0;
                                    ttmu = transmu[t-jj];
                                    tut = qutildavec + t + 1;
                                    ttmu += t+1-jj;
                                    i = s - t + 1;
                                    while (--i) {
                                        t1d += (*tut)*(*ttmu);
                                        ttmu++;
                                        tut++;
                                    }
                                    qyvec[t] = t1d;
                                    //t1d = -round(t1d);   
                                    t1d = -round(t1d- target_mud[t]);
                                    qutildavec[t] = qvvec[t] = t1d;
                                    qDeltavec[t] = 0;
                                    if (qutildavec[t] > -qyvec[t]) { 
                                        qdeltavec[t] = -1;
                                    } else {
                                        qdeltavec[t] = 1;
                                    }

                                    if (t+1==jjs) {
                                        //a vector whose projection is under the pruning function is found
                                        //search the empty slot
                                        
                                        //all-zero-check
                                        int isum = 0;
                                        for (i=jjs;i<=kk;i++) isum += abs(qutildavec[i]);
                                        //if isum==0 then it is the same as the sublattice
                                        
                                        if ((qctilda[t+1]!=0) && (isum>0)) {
                                            //the subtree of qctilda(projective length)=0 is already searched in ENUM0
                                            while (sos[vindex]!=-1){
                                                vindex++;
                                                if (vindex==slotmax) {
                                                    vindex = 0;
                                                    ea1breakflag=1;
                                                }
                                            }
                                            //store the state into slot
                                            for (i=jj;i<=kk+1;i++) {
                                                soctilda[vindex][i-jj+1] = qctilda[i];
                                                sovvec[vindex][i-jj+1] = qvvec[i];
                                                soyvec[vindex][i-jj+1] = qyvec[i];
                                                souvec[vindex][i-jj+1] = quvec[i];
                                                soutildavec[vindex][i-jj+1] = qutildavec[i];
                                                sodeltavec[vindex][i-jj+1] = qdeltavec[i];
                                                soDeltavec[vindex][i-jj+1] = qDeltavec[i];
                                            }

                                            sos[vindex] = s;    //set the slag for ENUM-1b process
                                            vindex++;
                                            if (vindex==slotmax) ea1breakflag=1;
                                        }
                                        cflag=0;
                                    }
                                }
                            } 

                            if (cflag==0) {
                                t++;
                                if (s>t) {
                                    qDeltavec[t] = -qDeltavec[t];
                                } else {
                                    if (s>=t) { 
                                    } else {
                                    }
                                    s = t;
                                }
                                qDeltavec[t] += (1+((qDeltavec[t]*qdeltavec[t]) >> 31))*qdeltavec[t];
                                qutildavec[t] = qvvec[t] + qDeltavec[t];
                            }
                            lmct++; //counter

                            //Does sufficient number of vectors stored?
                            if (ea1breakflag==1) {
                                current_used_slot=0;
                                for (i=0;i<slotmax;i++) {
                                    if (sos[i]!=-1) current_used_slot++;
                                }
                                if (current_used_slot>=slotmax*0.9) {
                                    break;
                                }
                                vindex = 0;
                                ea1breakflag=0;
                            }
                            t += endflag;
                        } //End of ENUM1a-core

                        sbk = s;
                        tbk = t;
                        if (t==kk+1) {
                            //ENUM1a is finished
                            ea1flag=2;
                        } else {
                            ea1flag=0;  //no thread works for ENUM-1a
                        }
                        #pragma omp critical
                        {
                            ea1thread=-1;
                            #pragma omp atomic    
                            gmct += lmct / 1000000;
                            lmct %= 1000000;
                        }            
                    }  //End of ENUM-1a
                }            

                //Start of ENUM-1b
                vindex = mythread;

                //sos=-1 <=> slot is empty
                //sos>0 <=> slot is waiting for process
                //-10000<sos<-100 <=> slot is reserved for process
                //sos=-29999 <=> slot is processed

                do {
                    if (ea1flag==2) {
                        //if ENUM-1a is finished, search an empty slot incrementaly
                        #pragma omp critical
                        {
                            vindex = 0;
                            while (vindex<slotmax) {
                                if (sos[vindex]>0) {
                                    sos[vindex]-=10000;
                                    break;
                                }
                                vindex++;
                            }
                        }
                    } else {

                        if (sos[vindex]>0) {
                            sos[vindex] -=10000;
                        }
                    }

                    if ((vindex<slotmax) && (-10000<sos[vindex]) && (sos[vindex]<-100) && (endflag==0)) {  

                        //Start of ENUM-1b for the slot
                        //cout << "start slot " << vindex << endl;
                        t = jjs-1;
                        s = sos[vindex]+10000;
                        sos[vindex]=-29999;


                        #pragma omp critical
                        {
                        vcount++;
                        }

                        //memory copy
                        double* putildavec = soutildavec[vindex]-jj+1; 
                        double* pyvec = soyvec[vindex]-jj+1; 
                        double* pvvec = sovvec[vindex]-jj+1; 
                        double* pctilda = soctilda[vindex]-jj+1; 
                        long int* pdeltavec = sodeltavec[vindex]-jj+1;
                        long int* pDeltavec = soDeltavec[vindex]-jj+1;

                        //Start of ENUM-1b-core
                        while (t < jjs) {
                            cflag = 0;
                            dtemp = (pyvec[t]+putildavec[t]-target_mud[t]) * cd[t];
                            dtemp *= dtemp;
                            pctilda[t] = pctilda[t+1] + dtemp;
/*
                            for (i=jj;i<=s;i++) {
                                cout << putildavec[i] << " ";
                                if (i==t) cout << " * ";
                            }
                            cout << pctilda[t+1] << " " << pctilda[t] << " " << pftemp[t] << " " << s << endl;
*/
                            if (pctilda[t] <= pftemp[t]) {
                                if (t<=jj+slide) {
                                    //cout << "t2:" << pctilda[t] << " "<< pftemp[t] << endl;
                                    if ((enummode==enum_mode_find_shortest) || (enummode==enum_mode_find_abort)) {
                                        //if (pctilda[t]<cd[t]*cd[t]*0.9999) {
                                            //Store the found vector to local memory 
                                            mlen = storetomemory(foundlen, foundvec, fnum, putildavec, pctilda[t], t, kk,yvec,mu,jj,pdeltavec);
                                            if (mlen >= 0) {
                                                if (t-jj < slide) {
                                                    slide = t-jj;
                                                }
                                            }
                                            if ((mlen != 0 ) && (mlen < clim) && (slide==0)) {
                                                //update pruning function if the pruning radius is updated
                                                #pragma omp critical
                                                {
                                                    clim = mlen;
                                                    for (i=0;i<kk-jj+1;i++) pruning_func[i] = normal_pruning_func[i] * clim;
                                                    //for (i=0;i<0.25*(kk-jj);i++) pruning_func[i] = max(pruning_func[i],to_double(c[jj+i]));
                                                    if (enummode==enum_mode_find_abort) {
                                                        endflag = kk;
                                                    }
                                                }
                                            }
                                        //}
                                    }
                                    if (t==jj) {
                                        if (enummode==enum_mode_all_vectors) {
                                        //holding
                                            storetomemory_allmode(foundlen,foundvec,fnum,findex,putildavec,pctilda[t],jj,kk);
                                            //store v + k*b_1 with small norm
                                            double rr = pyvec[t]+putildavec[t] - target_mud[t];
                                            double sign = +1;
                                            double utback = putildavec[t];
                                            if ((rr>0) && (s>jj)) sign=-1;
                                            i=2;
                                            while (1) {
                                                dtemp = pyvec[t]+utback - target_mud[t]+ sign*(int)(i/2);
                                                dtemp *= cd[t];
                                                dtemp *= dtemp;
                                                pctilda[t] = pctilda[t+1]+dtemp;
                                                if (pctilda[t]< clim) {
                                                    putildavec[t] = utback - target_mud[t]+  sign*(int)(i/2);
                                                    //storetomemory_allmode(foundlen,foundvec,fnum,findex,putildavec,pctilda[t],jj,kk);
                                                } else {
                                                    break;
                                                }
                                                if (s>jj) {
                                                    i++;
                                                    sign=-sign;
                                                } else {
                                                    i += 2;
                                                }
                                                if (findex==fnum) break;
                                            }
                                            if (findex==fnum) {
                                                endflag = kk;
                                                break;    //buffer memory is full
                                            }
                                        }
                                    }
                                }

                                if (t != jj) {
                                    cflag = 1;
                                    t--;
                                    t1d = 0;
                                    ttmu = transmu[t-jj];
                                    tut = putildavec + t + 1;
                                    ttmu += t+1-jj;
                                    i = s - t + 1;
    #ifdef AVX_PRD
                                    __m256d ymm;
                                    double atmp[4] __attribute__((aligned(32)));
                                    i = s - t;
                                    ymm = _mm256_setzero_pd();
                                    while(i & 3){
                                        t1d += (*tut) * (*ttmu);
                                        ++tut, ++ttmu, --i;
                                    }
                                    while(i > 0){
                                        ymm = _mm256_add_pd(ymm, _mm256_mul_pd(_mm256_loadu_pd(tut), _mm256_loadu_pd(ttmu)));
                                    // ymm = _mm256_fmadd_pd(_mm256_loadu_pd(tut), _mm256_loadu_pd(ttmu), ymm);
                                        tut += 4, ttmu += 4, i -= 4;
                                    }
                                    _mm256_store_pd(atmp, ymm);
                                    for(int j = 0; j < 4; ++j) t1d += atmp[j];
    #else
                                    while (--i) {
                                        t1d += (*tut)*(*ttmu);
                                        ttmu++;
                                        tut++;
                                    }
    #endif
                                    pyvec[t] = t1d;
    #ifdef SSE_ROUND
                                    __m128d tmp;
                                    tmp = _mm_load_sd(&t1d);
                                    tmp = _mm_round_sd(tmp, tmp, 0);
                                    _mm_store_sd(&t1d, tmp);
                                    t1d = -t1d;
    #else
                                    t1d = -round(t1d- target_mud[t]);
    #endif
                                    putildavec[t] = pvvec[t] = t1d;
                                    pDeltavec[t] = 0;

                                    if (putildavec[t] > -pyvec[t]) { 
                                        pdeltavec[t] = -1;
                                    } else {
                                        pdeltavec[t] = 1;
                                    }
                                }
                            }

                            if (cflag==0) {
                                t++;
                                pDeltavec[t] = -pDeltavec[t];
                                pDeltavec[t] += (1+((pDeltavec[t]*pdeltavec[t]) >> 31))*pdeltavec[t];
                                putildavec[t] = pvvec[t] + pDeltavec[t];
                            }

                            //counter 
                            if ((++lmct)^100000000) {
                            } else {
                                #pragma omp critical
                                {
                                    gmct += lmct / 1000000;
                                    lmct %= 1000000;
                                    if (vl>=3) {
                                        cout << "ENUM-1b count=" << gmct*0.01 << "*10^8 subtree=" << vcount;
                                        cout << " clim=" << sqrt(clim) ;
                                        cout << " time="  << gettimeofday_sec() - es;
                                        cout << " speed=" << 1.0*gmct / (gettimeofday_sec() - es);
                                        cout << "\r";
                                        cout.flush();
                                    }
                                    if (gmct /100 > elim) {
                                        endflag = kk;
                                    }
                                }
                                t +=  endflag;
                            }

                        } // End of ENUM-1b core
                        sos[vindex]=-1;
                    }

                    if (ea1flag<=1) {
                        //ENUM-a1 is not finished
                        vindex += numthreads;
                    } else {
                        //ENUM-a1 is finished
                        vindex++;
                    }
                } while (vindex<slotmax);


                //if ENUM1a was finished
                if (ea1flag==2) {
                    if (finishmode==finish_mode_nowaste) {
                        //counting current "waiting" slot
                        current_used_slot=0;
                        for (i=0;i<slotmax;i++) {
                            if (sos[i]>0) current_used_slot++;
                        }
                        if (current_used_slot==0) endflag = kk;
                    }
                    if (finishmode==finish_mode_exact) {
                        //counting current "empty, waiting and under processing" slot
                        current_used_slot=0;
                        for (i=0;i<slotmax;i++) {
                            if (sos[i]!=-1) current_used_slot++;
                        }
                        if (current_used_slot==0) endflag = kk;
                    }
                } 
            } while (endflag==0);
            endflag = kk;

            //summerize # processed nodes
            #pragma omp critical
            {
                gmct += lmct / 1000000;
                lmct %= 1000000;
            }

        } //End of paralleled part

        //Add the bias
        for (j=0;j<findex;j++) {
            //cout << "f[" << j << "]: ";
            for (i=jj;i<=kk;i++) {
                    foundvec[j][i] += bias[i];
                    ///cout << foundvec[j][i] << " ";
             }
             //cout << endl;
        }

        if (vl>=3) {
            cout << "ENUM count=" << gmct*0.01 << "*10^8 subtree=" << vcount;
            cout << " clim=" << sqrt(clim) ;
            cout << " time="  << gettimeofday_sec() - es;
            cout << " speed=" << 1.0*gmct / (gettimeofday_sec() - es) << "M/s ea1=" << (int)ea1flag << "   " << endl;
        }
        current_enum_speed = 1.0*gmct / (gettimeofday_sec() - es);
        current_processed = gmct;
    }



}

#endif
