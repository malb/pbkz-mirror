#ifndef _inc_vector_enumeration
#define _inc_vector_enumeration

//Definitions to control the behaviours of enumeration

#define enum_mode_find_shortest 0x01
//update pruning radius during enumeration

#define enum_mode_all_vectors 0x02  
//enumerate all vectors s.t. |v|<bound and under the pruning function

#define enum_mode_find_abort 0x04
//Abort the enumeration when a vector is found

#define enum_mode_find_short_projection 0x08 

#define enum_mode_except_trivial 0x10

#define finish_mode_nowaste 0x01
//In multithread-mode, when a non-working thread exists, then halt the subroutine

#define finish_mode_exact 0x02

//#slots (#holding vectors found in ENUM1a)        
#define slotmax 1000

#ifdef SSE_ROUND
#include <emmintrin.h>
#include <smmintrin.h>
#endif

#ifdef AVX_PRD
#include <immintrin.h>
#endif

extern void set_approximate_extreme_pruning(double* pf,quad_float* c,int bs,double prob,double delta=0);

void convert_intcoefftomat_zz(mat_ZZ& fm,mat_ZZ& L,int** foundvec,int findex,int offset=-1) {
    
    int i,j,k;
    int n = L.NumCols();
    fm.SetDims(findex,n);
    
    for (i=0;i<findex;i++) {
        for (k=0;k<n;k++) {
            fm[i][k] = 0;
            //cout << foundvec[i][k+1] << " ";
        }
        for (j=1;j<=n;j++) {
            for (k=0;k<n;k++) {
                fm[i][k] += foundvec[i][j+offset] * L[j-1][k];
            }
        }
    }
}



namespace lattice_enum {
    
    int initval=0;
    
    double** enum_share_double;
    long** enum_share_long;
    int* enumpara_share_int;
    double** enum_share_transmu;
    int** enum_share_foundvec;
    
    //memory for ENUM0
    double *ctilda, *vvec, *yvec, *uvec, *utildavec;
    long *Deltavec, *deltavec;
    
    //memory for ENUM1a
    double *qctilda, *qvvec, *qyvec, *quvec, *qutildavec;
    long *qDeltavec, *qdeltavec;
    
    //memory for ENUM1b
    double **soctilda, **sovvec, **soyvec, **souvec, **soutildavec;
    long **soDeltavec, **sodeltavec;
    int* sos;
    
    //some returned values
    double current_enum_speed;
    double* enum_speed;
    double current_processed;
    
    void initialize() {

        if (initval!=0) return;
        
        int i;

        ctilda = new double[latticemaxdim];
        vvec = new double[latticemaxdim];
        yvec = new double[latticemaxdim];
        uvec = new double[latticemaxdim];
        utildavec = new double[latticemaxdim];
        Deltavec = new long[latticemaxdim];
        deltavec = new long[latticemaxdim];
        
        qctilda = new double[latticemaxdim];
        qvvec = new double[latticemaxdim];
        qyvec = new double[latticemaxdim];
        quvec = new double[latticemaxdim];
        qutildavec = new double[latticemaxdim];
        qDeltavec = new long[latticemaxdim];
        qdeltavec = new long[latticemaxdim];
        
        soctilda = new double*[slotmax];
        sovvec = new double*[slotmax];
        soyvec = new double*[slotmax];
        souvec = new double*[slotmax];
        soutildavec = new double*[slotmax];
        soDeltavec = new long*[slotmax];
        sodeltavec = new long*[slotmax];
        for (i=0;i<slotmax;i++) {
            soctilda[i] = new double[latticemaxdim];
            sovvec[i] = new double[latticemaxdim];
            soyvec[i] = new double[latticemaxdim];
            souvec[i] = new double[latticemaxdim];
            soutildavec[i] = new double[latticemaxdim];
            soDeltavec[i] = new long[latticemaxdim];
            sodeltavec[i] = new long[latticemaxdim];
        }

        sos = new int[slotmax];


        enum_share_double = new double*[16];
        enum_share_long = new long*[8];
        for (i=0;i<16;i++) enum_share_double[i] = new double[latticemaxdim];
        for (i=0;i<8;i++) enum_share_long[i] = new long[latticemaxdim];

        enum_share_transmu = new double*[latticemaxdim];
        for (i=0;i<latticemaxdim;i++) enum_share_transmu[i] = new double[latticemaxdim];

        enum_share_foundvec = new int*[maxfoundvec];
        for (i=0;i<maxfoundvec;i++) enum_share_foundvec[i] = new int[latticemaxdim];

        enum_speed = new double[128];
        for (i=0;i<128;i++) enum_speed[i] = -1;

        initval = 1;
    }

    void storetomemory_allmode(double* foundlen,int** foundvec,int fnum,int& findex,double* utildavec,double ctilda,int jj,int kk) {

            if (ctilda==0) return;
            if (findex == fnum) return;

            int i;
            for (i=jj;i<=kk;i++) {
                foundvec[findex][i] = utildavec[i];
            }
            foundlen[findex] = ctilda;

            findex++;
    }

    int msback = 0;

    double storetomemory(double* foundlen, int** foundvec, int fnum, double* utildavec, double ctilda, int jj, int kk,double* yvec, quad_float** mu,int jjorg,long int* deltavec,char enummode = 0) {

            if( ctilda==0) return -1;

            int i,j;
            double cmu;

            //variable name has been changed, need to execute to check 2016.3
            double* utt = enum_share_double[15];
            long* fvcoeff = enum_share_long[6];

            for (i=jjorg;i<kk+1;i++) utt[i] = 0;
            for (i=jjorg;i<kk+1;i++) fvcoeff[i] = 0;

            for (i=kk;i>=jj;i--) {
                for (j=0;j<i;j++) {
                    utt[j] += utildavec[i] * to_double(mu[i][j]);
                }
                utt[i] += utildavec[i];
                fvcoeff[i] = utildavec[i];
            }
            for (i=jj-1;i>=jjorg;i--) {
                cmu = round(utt[i]);
                for (j=0;j<i;j++) {
                    utt[j] -= cmu * to_double(mu[i][j]);
                }
                utt[i] -= cmu;
                fvcoeff[i] -= cmu;      
            }

            //is it a trivial vector?
            if ((enummode & enum_mode_except_trivial) != 0 ) {
                j=0;
                for (i=jj;i<=kk;i++) {
                    if (fvcoeff[i]!=0) j++;
                }
                if (j==1) return -1;
            }

            //search empty slot
            int ms = msback;
            ms %= fnum;
            char found=0;
            for (i=0;i<fnum;i++) {
                if (foundlen[ms] == -1) {
                    found = 1;
                    break;
                }
                ms = (ms+1)%fnum;
            }
            msback = ms;
            double maxlen = 0;
            int maxjj;

            if (found==0) {
                maxlen = foundlen[0];
                maxjj = foundvec[0][0];
                ms = 0;
                for (i = 1; i < fnum; i++) {
                    if ( (foundvec[i][0] > maxjj) || ((foundvec[i][0] == maxjj) && (maxlen < foundlen[i]))) {
                        maxjj = foundvec[i][0];
                        maxlen = foundlen[i];
                        ms = i;
                    }
                }
            }

            for (i = jjorg; i <= kk; i++) {
                foundvec[ms][i] = fvcoeff[i];
            }
            foundvec[ms][0] = jj;
            foundlen[ms] = ctilda;
            for (i = 1; i < jjorg; i++) foundvec[ms][i] = 0;

            maxlen = 0;
            for (i=0;i<fnum;i++) {
                if (foundvec[i][0]==jjorg) {
                    if ((maxlen > foundlen[i]) ||(maxlen==0)) {
                        maxlen = foundlen[i];
                    }
                }
            }
            return maxlen;
    }

    RR Rigid_upper_cost_enuma1(double* R,quad_float* c,int n,double radius_c) {
        //pruning_func[1,...,n]
        //c[1,...,n]
        //return upper # vectors found in ENUM-1a
        RR sum;
        RR t,t2;
        int i,k;
        double* R2;
        R2 = new double[n/2+2];
        t2 = 1;
        k = n;
        if (R[k]!=0) {
            for (i=1;i<=k/2;i++) R2[i] = R[2*i]/R[k];
            t = EvenSimplex::EvenSimplexVolumeqf(k/2,R2,opt_volume_prob);
            t *= VolumeBall(k,sqrt(R[k]));
            for (k=2;k<=n;k+=2) {
                t2 *= radius_c * radius_c / sqrt(to_double(c[n-k+1])*to_double(c[n-k+2]));
            }
            sum = t * t2;
        }
        delete [] R2;
        return 0.5 * sum; 
    }



    void ENUMCore(int** foundvec,double* foundlen,int fnum,int& findex,int& slide,int jj,int bs,long int elim,double clim,double* normal_pruning_func,quad_float** mu,quad_float* c,int vl,char enummode) {
        //stop after elim*10^8 nodes are processed


         initialize();
         char cflag;
         double t1d;
         double* ttmu;
         double* tut;
         double dtemp;

         double *ctilda = enum_share_double[0];
         double *vvec = enum_share_double[1];
         double *yvec = enum_share_double[2];
         double *uvec = enum_share_double[3];
         double *utildavec = enum_share_double[4];
         long *Deltavec = enum_share_long[0];
         long *deltavec = enum_share_long[1];


         double es;
         es = gettimeofday_sec();   //timer


         int kk = jj+bs-1;
         int s,t;
         int i,j;
         slide = bs;
         if (enummode &  enum_mode_all_vectors) slide=0;
         msback = 0;

         int enum_cnt=0;
         long int totalenum = 0;
         double mlen;       //minimum length of found vector

         double* pruning_func = enum_share_double[7];
         double* pftemp = pruning_func - jj;
         for (i=0;i<kk-jj+1;i++) {
             pruning_func[i] = normal_pruning_func[i] * clim;
             //cout << i << " " << pruning_func[i] << endl;
         }

         if (enummode & enum_mode_find_shortest)  {
             for (i=0;i<0.25*(kk-jj);i++) {
                 pruning_func[i] = max(pruning_func[i],to_double(c[jj+i]));
             }
         }
         //initialize memory for vector holding
         for (j=0;j<fnum;j++) {
             foundlen[j] = -1;
             foundvec[j][0] = jj + bs + 1;
             for (i=jj;i<=kk;i++) foundvec[j][i] = 0; 
         }

         //transmu
         double** transmu = enum_share_transmu;
         double* cd = enum_share_double[5];
         for (i=0;i<bs;i++) {
             for (j=i;j<bs;j++) {
                 transmu[i][j] = to_double(mu[j+jj][i+jj]);
             }
         }
         //cd[i]=|b*_i|
         for (i=0;i<bs;i++) cd[i+jj] = sqrt(to_double(c[i+jj]));


         //initialize variables for lattice enumeration
         utildavec[jj] = uvec[jj] = 1;
         yvec[jj] = vvec[jj] = 0;
         Deltavec[jj] = 0;
         s = t = jj;
         deltavec[jj] = 1;   
         for (i = jj+1; i <= kk+1; i++) {
            ctilda[i] = uvec[i] = utildavec[i] = yvec[i] = 0;
            Deltavec[i] = 0;
            vvec[i] = 0;
            deltavec[i] = 1;
         }
         //utildavec[0]=0;
        if (vl>=3) {
            cout << "Start ENUM-0 (initial radius=" << sqrt(clim) << ")" << endl;
        }
         //ENUM main
         while (t <= kk) {
            dtemp = (yvec[t]+utildavec[t]) * cd[t];
            dtemp *= dtemp;
            ctilda[t] = ctilda[t+1] + dtemp;
            cflag = 0;

            //cout << t << " " << ctilda[t] << " " <<  pftemp[t] << " " << enum_cnt << endl;
            if (ctilda[t] <= pftemp[t]) {
                if (t<=jj+slide) {
                    //a vector s.t. pi_{i}|v|<|b*_i| (i=jj+slide+4) is found
		  if ( ((enummode & enum_mode_find_shortest)+ enum_mode_find_abort + enum_mode_find_short_projection) != 0) {
                        if ( (enummode & enum_mode_find_short_projection) != 0) {
                            //cout << "store: " << t << endl;
                            #pragma omp critical
                            {
                                mlen = storetomemory(foundlen, foundvec, fnum, utildavec, ctilda[t], t, kk,yvec,mu,jj,deltavec,enummode);
                            }
                            if (mlen >= 0) {
                                //cout << "find: slide=" << t << endl;
                                if (t-jj+6 < slide) {
                                    slide = max(0,t-jj+6);
                                    //cout << "newslide=" << slide << endl;
                                }
                            }
                        } else                         
                        if (ctilda[t]<cd[t]*cd[t]*0.9999) {
                            //Store the found vector to local memory 
                            mlen = storetomemory(foundlen, foundvec, fnum, utildavec, ctilda[t], t, kk,yvec,mu,jj,deltavec,enummode);
                            if (mlen >= 0) {
                                if (t-jj < slide) {
                                    slide = t-jj;
                                }
                            }
                            if (mlen>0) {
                                if ((mlen < clim) && (slide==0)) {
                                    //update pruning function if the pruning radius is updated
                                    clim = mlen;
                                    //cout << "mlen=" << mlen << endl;
                                    for (i=0;i<kk-jj+1;i++) {
                                        pruning_func[i] = normal_pruning_func[i] * clim;
                                        //cout << i << " " << pruning_func[i] << endl;
                                    }
                                    for (i=0;i<0.25*(kk-jj);i++) {
                                        pruning_func[i] = max(pruning_func[i],to_double(c[jj+i]));
                                        //cout << i << " " << pruning_func[i] << endl;
                                    }
                                    if (enummode & enum_mode_find_abort) {
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    if (t==jj) {
                        if (enummode & enum_mode_find_abort) {
                            mlen = storetomemory(foundlen, foundvec, fnum, utildavec, ctilda[t], t, kk,yvec,mu,jj,deltavec,enummode);
                            break;
                        }
                        
                        if (enummode & enum_mode_all_vectors) {
                            storetomemory_allmode(foundlen,foundvec,fnum,findex,utildavec,ctilda[t],jj,kk);
                            //store v + k*b_1 with small norm
                            double rr = yvec[t]+utildavec[t];
                            double sign = +1;
                            double utback = utildavec[t];
                            if ((rr>0) && (s>jj)) sign=-1;
                            i=2;
                            while (1) {
                                dtemp = yvec[t]+utback + sign*(int)(i/2);
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
                            }
                            if (findex==fnum) break;    //buffer memory is full
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
                    __m128d tmp;
                    tmp = _mm_load_sd(&t1d);
                    tmp = _mm_round_sd(tmp, tmp, 0);
                    _mm_store_sd(&t1d, tmp);
                    t1d = -t1d;
    #else
                    t1d = -round(t1d);
    #endif
                  //the above is equivalent to the following in LLL_QP.c
                  //t1 = -t1;
                  //  if (t1 >= 0)
                  //t1 = ceil(t1-0.5);
                  //  else
                  //t1 = floor(t1+0.5);

                  utildavec[t] = vvec[t] = t1d;
                  Deltavec[t] = 0;

                  //To avoid conditional jump

                // __m128d tmp1, tmp2;
                // tmp1 = _mm_loadh_pd(tmp1, utildavec + t);
                // tmp2 = _mm_loadh_pd(tmp2, yvec + t);
                // tmp1 = _mm_add_pd(tmp1, tmp2);
                // deltavec[t] = 1 - (~_mm_movemask_pd(tmp1) & 2);

                  if (utildavec[t] > -yvec[t]) 
                     deltavec[t] = -1;
                  else
                     deltavec[t] = 1;
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
               //The above is equivalent to the following in LLL_QP.c (aboiding if-sentence)
               //if (Deltavec[t]*deltavec[t] >= 0) Deltavec[t] += deltavec[t];

               utildavec[t] = vvec[t] + Deltavec[t];
            }

            //counter
            if ((++enum_cnt)^100000000) {

            } else {
                enum_cnt = 0;
                totalenum++;
                if (totalenum>elim) {
                    break;
                }
                if ((enummode & enum_mode_all_vectors) && (vl>=1)) {
                    cout << "ENUM subdim=" << s-jj+1 << "/" << kk-jj+1 << " count=" << totalenum << "*10^8 ";
                    cout << "#found vector=" << findex << " time="  << gettimeofday_sec() - es << "\r";
                    cout.flush();
                }
                if (vl>=3) {
                    cout << "ENUM count=" << totalenum << "*10^8 ";
                    cout << " clim=" << sqrt(clim);
                    cout << " time="  << gettimeofday_sec() - es;
                    cout << " speed=" << 100.0*totalenum / (gettimeofday_sec() - es) << "Mnodes/sec" << "\r";
                    cout.flush();
                }
            }
         }  //end of enum
        if ((enummode & enum_mode_all_vectors) && (vl>=1)) {
            cout << endl;
        }
         
        //Count # found vector
        findex = 0; 
        for (i=0;i<fnum;i++) {
           if (foundlen[i] != -1) {
               findex++;
           }
        }

        es = gettimeofday_sec() - es;
        if (vl>=3) {
            cout << "#nodes=" << 0.1*totalenum + 0.000000001* enum_cnt << "Gnodes ";
            cout << "Time=" << es << " enum_speed=" << (100.0*totalenum+enum_cnt*0.000001) / es << "Mnodes/sec ";
            cout << "final_radius=" << sqrt(clim) << endl;
        }
        current_enum_speed = (100.0*totalenum+enum_cnt*0.000001) / es;
        current_processed = 100.0*totalenum+enum_cnt*0.000001;
    }





    void ENUMCoremt(int** foundvec,double* foundlen,int fnum,int& findex,int& slide,int jj,int bs,long int elim,double clim,double* normal_pruning_func,quad_float** mu,quad_float* c,int vl,int numthreads,char enummode,char finishmode) {

        //pruningfunction[1...bs] in monotonically decreasing

        //stop after elim*10^8 nodes are processed
        //ENUM0 = enumeration for L(b1,...,bk)
        //ENUM1a = enumeration for L(b[k+1],...,bn)
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
        if (enummode & enum_mode_all_vectors) slide=0;

        double* pruning_func = enum_share_double[7];
        double* pftemp = pruning_func - jj;
        for (i=0;i<bs;i++) pruning_func[i] = normal_pruning_func[i] * clim;
        if ((enummode& enum_mode_find_shortest) || (enummode & enum_mode_find_abort)) {
            for (i=0;i<0.25*(kk-jj);i++) pruning_func[i] = max(pruning_func[i],to_double(c[jj+i])); 
        }
        //initialize memory for vector holding
        for (j=0;j<fnum;j++) {
            foundlen[j] = -1;
            foundvec[j][0] = jj + bs + 1;
            for (i=jj;i<=kk;i++) foundvec[j][i] = 0; 
        }

        //transmu
        double** transmu = enum_share_transmu;

        double* cd = enum_share_double[5];
        for (i=0;i<bs;i++) {
            for (j=i;j<bs;j++) {
                transmu[i][j] = to_double(mu[j+jj][i+jj]);
            }
        }
        //cd[i]=|b*_i|
        for (i=0;i<bs;i++) cd[i+jj] = sqrt(to_double(c[i+jj]));

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


        unsigned int lmct;  //local enum count
        long long int gmct = 0;   //global enum count
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

                while (t <= jjs-1) {
                    
                    dtemp = (yvec[t]+utildavec[t]) * cd[t];
                    dtemp *= dtemp;
                    ctilda[t] = ctilda[t+1] + dtemp;
                    cflag = 0;

                    if (ctilda[t] <= pftemp[t]) {
                        if (t<=jj+slide) {
                            //a vector s.t. pi_{i}|v|<|b*_i| (i=jj+slide+4) is found
                            if ((enummode & enum_mode_find_shortest) || (enummode & enum_mode_find_abort)) {
                                if (ctilda[t]<cd[t]*cd[t]*0.9999) {
                                    //Store the found vector to local memory 
                                    #pragma omp critical
                                    {
                                        mlen = storetomemory(foundlen, foundvec, fnum, utildavec, ctilda[t], t, kk,yvec,mu,jj,deltavec,enummode);
                                    }
                                    if (mlen >= 0) {
                                        if (t-jj < slide) {
                                            slide = t-jj;
                                        }
                                    }
                                    if ((mlen > 0 ) && (mlen < clim) && (slide==0)) {
                                        //update pruning function if the pruning radius is updated
                                        #pragma omp critical
                                        {
                                            clim = min(clim,mlen);
                                            for (i=0;i<kk-jj+1;i++) pruning_func[i] = normal_pruning_func[i] * clim;
                                            for (i=0;i<0.25*(kk-jj);i++) pruning_func[i] = max(pruning_func[i],to_double(c[jj+i]));
                                            if (enummode & enum_mode_find_abort) {
                                                endflag = kk;
                                            }
                                        }
                                    }
                                }
                            }
                            if (t==jj) {
                                if (enummode & enum_mode_find_abort) {
                                    #pragma omp critical
                                    {
                                        mlen = storetomemory(foundlen, foundvec, fnum, utildavec, ctilda[t], t, kk,yvec,mu,jj,deltavec,enummode);
                                    }
                                    endflag = kk;
                                }
                                if (enummode & enum_mode_all_vectors) {
                                //holding
                                    #pragma omp critical
                                    {
                                        storetomemory_allmode(foundlen,foundvec,fnum,findex,utildavec,ctilda[t],jj,kk);
                                    }
                                    //store v + k*b_1 with small norm
                                    double rr = yvec[t]+utildavec[t];
                                    double sign = +1;
                                    double utback = utildavec[t];
                                    if ((rr>0) && (s>jj)) sign=-1;

                                    i=2;
                                    while (1) {
                                        dtemp = yvec[t]+utback + sign*(int)(i/2);
                                        dtemp *= cd[t];
                                        dtemp *= dtemp;
                                        ctilda[t] = ctilda[t+1]+dtemp;
                                        if (ctilda[t]< clim) {
                                            utildavec[t] = utback +  sign*(int)(i/2);
                                            #pragma omp critical
                                            {
                                                storetomemory_allmode(foundlen,foundvec,fnum,findex,utildavec,ctilda[t],jj,kk);
                                            }
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
                            __m128d tmp;
                            tmp = _mm_load_sd(&t1d);
                            tmp = _mm_round_sd(tmp, tmp, 0);
                            _mm_store_sd(&t1d, tmp);
                            t1d = -t1d;
    #else
                            t1d = -round(t1d);
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
                            dtemp = (qyvec[t]+qutildavec[t]) * cd[t];
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
                                    t1d = -round(t1d);   
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
                                        if (qctilda[t+1]!=0) {
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
                            dtemp = (pyvec[t]+putildavec[t]) * cd[t];
                            dtemp *= dtemp;
                            pctilda[t] = pctilda[t+1] + dtemp;

                            if (pctilda[t] <= pftemp[t]) {
                                if (t<=jj+slide) {
                                    //a vector s.t. pi_{i}|v|<|b*_i| (i=jj+slide+4) is found
                                    if ((enummode &  enum_mode_find_shortest) || (enummode & enum_mode_find_abort)) {
                                        if (pctilda[t]<cd[t]*cd[t]*0.9999) {
                                            //Store the found vector to local memory 
                                            #pragma omp critical
                                            {
                                                mlen = storetomemory(foundlen, foundvec, fnum, putildavec, pctilda[t], t, kk,yvec,mu,jj,pdeltavec,enummode);
                                            }
                                            if (mlen >= 0) {
                                                if (t-jj < slide) {
                                                    slide = t-jj;
                                                }
                                            }
                                            if ((mlen > 0 ) && (mlen < clim) && (slide==0)) {
                                                //update pruning function if the pruning radius is updated
                                                #pragma omp critical
                                                {
                                                    clim = min(clim,mlen);
                                                    for (i=0;i<kk-jj+1;i++) pruning_func[i] = normal_pruning_func[i] * clim;
                                                    for (i=0;i<0.25*(kk-jj);i++) pruning_func[i] = max(pruning_func[i],to_double(c[jj+i]));
                                                    if (enummode & enum_mode_find_abort) {
                                                        endflag = kk;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    if (t==jj) {
                                        if (enummode & enum_mode_find_abort) {
                                            #pragma omp critical
                                            {
                                                mlen = storetomemory(foundlen, foundvec, fnum, putildavec, pctilda[t], t, kk,yvec,mu,jj,pdeltavec,enummode);
                                            }
                                            endflag = kk;
                                        }
                                        if (enummode &  enum_mode_all_vectors) {
                                        //holding
                                            #pragma omp critical
                                            {
                                                storetomemory_allmode(foundlen,foundvec,fnum,findex,putildavec,pctilda[t],jj,kk);
                                            }
                                            //store v + k*b_1 with small norm
                                            double rr = pyvec[t]+putildavec[t];
                                            double sign = +1;
                                            double utback = putildavec[t];
                                            if ((rr>0) && (s>jj)) sign=-1;
                                            i=2;
                                            while (1) {
                                                dtemp = pyvec[t]+utback + sign*(int)(i/2);
                                                dtemp *= cd[t];
                                                dtemp *= dtemp;
                                                pctilda[t] = pctilda[t+1]+dtemp;
                                                if (pctilda[t]< clim) {
                                                    putildavec[t] = utback +  sign*(int)(i/2);
                                                    #pragma omp critical
                                                    {
                                                        storetomemory_allmode(foundlen,foundvec,fnum,findex,putildavec,pctilda[t],jj,kk);
                                                    }
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
                                    t1d = -round(t1d);
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
                                        cout << " speed=" << 1.0*gmct / (gettimeofday_sec() - es) << "M/s slide=" << (int)slide;
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

        //Count # found vector
        findex = 0; 
        for (i=0;i<fnum;i++) {
           if (foundlen[i] != -1) {
               findex++;
           }
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


    int ENUM(mat_ZZ& L,double* foundlen,int** foundvec,int buff,double prob,double alpha,int elim,int vl,int opt,int parallel=1,char enummode=enum_mode_all_vectors,char finishmode=finish_mode_nowaste) {
        //return value is # found vectors
        //foundvec[i][j] is the coeff. of b_j in i-th vector

        double clim;
        double* normal_pruning_func;
        quad_float** mu;
        quad_float* c;

        int i;
        int m = L.NumRows();
        int findex = 0;

        normal_pruning_func = new double[m+1];
        mu = new quad_float*[m+1];
        for (i=0;i<m+1;i++) mu[i] = new quad_float[m+1];
        c = new quad_float[m+1];

         if (vl>=1) {
             cout << "Start ENUM-" << m << " prob=" << prob << " memory=" << buff << endl;
         }

        int slide = 0;
        GramSchmidt(L,mu,c,m);
        clim = alpha * lattice_tools::LatticeGH(c-1,m,INPUT_SQUARED);
        if (vl>=1) {
            cout << "Init clim=" << clim << endl;
        }
        set_approximate_extreme_pruning(normal_pruning_func,c,m,prob,0);

        //cout << "enum-start" << endl;
        if (parallel>=2) {
            //int parallel = omp_get_max_threads();
            ENUMCoremt(foundvec,foundlen,buff,findex,slide,0,m,elim,clim*clim,normal_pruning_func,mu,c,vl,parallel,enum_mode_all_vectors,finish_mode_exact);
        } else {
            //outputs are foundvec and foundlen
            ENUMCore(foundvec,foundlen,buff,findex,slide,0,m,elim,clim*clim,normal_pruning_func,mu,c,vl,enum_mode_all_vectors);
        }
        //cout << "enum-end" << endl;
        for (i=0;i<m+1;i++) delete [] mu[i];
        delete [] c;
        delete [] mu;
        return findex;
    }

    int ENUM(mat_ZZ& L,mat_ZZ& foundmat,double prob,double alpha,int elim,int vl,int opt,int parallel=1,char enummode=enum_mode_all_vectors,char finishmode=finish_mode_nowaste) 
    {
        initialize();
        int findex=0; //# found vectors
        int fnum=maxfoundvec;       //maxfoundvec is defined in the header file
        double* foundlen = new double[fnum];
        int** foundvec = enum_share_foundvec;
        findex = ENUM(L,foundlen,foundvec,fnum,prob,alpha,elim,VL3,0,parallel,enummode,finishmode);

        convert_intcoefftomat_zz(foundmat,L,foundvec,findex);
        delete [] foundlen;
        return findex;
    }
    //#define debug_skipdetectspeed
    double enum_speed_bench(int numthreads) {

        //returned speed/sec in Mnodes
        
        
        initialize();
        
        
        int ns = numthreads;
        if (ns<1) ns = 1;
        
        if (enum_speed[ns]>0) return enum_speed[ns];

        //Read from cache
        std::string speedstring;
        speedstring = ReadConf("enumspeed.cache","enum_th"+to_string(ns));
        if (speedstring!="")  {
            enum_speed[ns] = atof(speedstring.c_str());
            if (enum_speed>0) return enum_speed[ns];
        }
        cout << "Detecting lattice enumeration speed: # threads=" << numthreads << endl;
        
        mat_ZZ L;
        int n = 120;
        int q = 211;
        int i,j;
        L.SetDims(n,n);

        //test lattice
        for (i=0;i<n/2;i++) {
            L[i][i] = q;
            L[i+n/2][i+n/2] = 1;
            for (j=0;j<n/2;j++) {
                L[i+n/2][j] = rand() % q;
            }
        }

        
        int backv = verbose;
        verbose = 0;
        LLL_QP(L,0.999,0,0,0);    
        verbose = backv;
        
        int buff=200;
        double* foundlen;
        int** foundvec;
        foundlen = new double[buff+1];
        foundvec = new int*[buff+1];
        for (i=0;i<buff+1;i++) foundvec[i] = new int[n+1];

    #ifndef debug_skipdetectspeed
        if (numthreads>1) {
             omp_set_num_threads(numthreads); //setting number of threads        
             ENUM(L,foundlen,foundvec,buff,0.02*numthreads,1.00,1*numthreads,VL0,OPT_MULTITHREAD,numthreads);
        } else {
             ENUM(L,foundlen,foundvec,buff,0.02,1.00,1,VL0,0);
        }
        enum_speed[ns] = current_enum_speed;
    #else 
        enum_speed[ns] = max(10 * numthreads,10);
    #endif
        cout << "speed=" << current_enum_speed << "Mnodes/sec" << endl;
        for (i=0;i<buff+1;i++) delete [] foundvec[i];
        delete [] foundlen;
        delete [] foundvec;

        WriteConf("enumspeed.cache","enum_th"+to_string(ns),to_string(enum_speed[ns]));

        return enum_speed[ns];

    }
} // end of namespace lattice_enum



#endif
