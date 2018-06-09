#ifndef _inc_sampling_tools
#define _inc_sampling_tools

//10 for fast implementation for BKZ
//10000 for a near-exact estimation
#define SampleMult 10

//Sampling from objects
namespace sampling_tools{
    
    int init=0;
    double** sample_share_double;
    unsigned int* rstates;
    
    void initialize() {
        if (init!=0) return;
        
        int i;
        int numthread = omp_get_num_procs();
        sample_share_double = new double*[4*numthread+15];
        for (i=0;i<4*numthread+15;i++) sample_share_double[i] = new double[200*numthread];
        //shared area for generating random numbers
        rstates = new unsigned int[5*numthread];
        for (i=0;i<5*numthread;i++) {
            rstates[i] = (int(time(NULL))  + i*100)%0x40000000;
            //cout << "rs=" << rstates[i] << endl;
        }
        init =1;
    }

    
    //generating pseudo random real number from [0...1] 
    //very lazy random number generator
    double drand(unsigned int* th) {
        th[1] = th[1]*0x112233 + 0x9bcd231 + th[2];
        th[2] = th[2]*0x332211 + 0xf7abc10  +th[1];
        th[1] &= 0x3fffffff;
        return (double)(th[1]/1073741824.0);
    }

    //generating random real number from [a...b]
    double unique(double a,double b,unsigned int* th) {
        return a+(b-a)*drand(th);
    }

    //Generating normal distribution 
    double normal_dist(double mean,double sigma,unsigned int* th) {
        th[1] = th[1]*0x112233 + 0x9bcd231 + th[2];
        th[2] = th[2]*0x332211 + 0xf7abc10  +th[1];
        th[1] &= 0x3fffffff;
        th[2] &= 0x3fffffff;
        return mean + sigma * sqrt(-2*log(th[1]/1073741824.0)) * cos(2*3.141592653589*(th[2]/1073741824.0));
    }
    
    
    bool isContainedDKLR(double* y,double* L,double* R,int n,int l) {
        //Check: y[k] in [0,1] for k=1,...,l and
        //          L[k] < \sum y[i] < R[k] for k=1,...,l
        int i;
        double sum=0;
        for (i=1;i<=l;i++) {
            //cout << i << " " << y[i];
            if (y[i]<0) return false;
            if (y[i]>1) return false;
            sum += y[i];
            //cout << " " << sum << " " << R[i] << " " << L[i] << endl;
            if (sum>R[i]) return false;
            if (sum<L[i]) return false;
        }
        for (i=l+1;i<=n;i++) {
            if (y[i]<0) return false;
            if (y[i]>1) return false;
        }
        return true;
    }


    void SamplefromDKLR(double* y,double*L,double* R,int n,int l,unsigned int* th) {

        //Uniform sampling from the object (y[1],...,y[n]) s.t. 
        //  y[k] in [0,1] for k=1,...,n and L[k] < \sum_{i=1,....,k} y[i] < R[k] for k=1,...,n
        // by the hit-and-run algorithm
        int i,j=0;
        double* d = y + n + 1;
        double tl,tu;
        double s,ys;
        double rr=0;

        s = 0;
        //setting initial point in the object
        while (isContainedDKLR(y,L,R,n,l) == false) {
            y[1] = unique(L[1],R[1],th);
            s = y[1];
            for (i=2;i<=l;i++) {
                y[i] = unique(max(L[i]-s,0.0),max(0.95*R[i]-s,0.0),th);
                s += y[i];
            }
            for (i=l+1;i<=n;i++) y[i] = unique(0,1,th);
            j++;
            if (j>10) {
                j=100;
            }
        }

        do {
            rr += 1.0;
            if (rr>3) {
                rr = 0.01;
                do {
                    y[1] = unique(L[1],R[1],th);
                    s = y[1];
                    for (i=2;i<=l;i++) {
                        y[i] = unique(max(L[i]-s,0.0),0.9*R[i]-s,th);
                        s += y[i];
                    }
                    for (i=l+1;i<=n;i++) y[i] = unique(0,1,th);
                    j++;
                    if (j>10) {
                        cout << "j error?" << endl; 
                        for (i=1;i<=n;i++) {
                            cout << "RL[" << i << "]=" << R[i] << " " << L[i] << endl;
                        }
                        
                    }
                } while (isContainedDKLR(y,L,R,n,l) == false);
                rr = 0;
            }

            //random direction
            s = 0;
            for (i=1;i<=n;i++) {
                d[i] = normal_dist(0,1,th);
                s += d[i]*d[i];
            }
            s = sqrt(s);
            for (i=1;i<=n;i++) {
                d[i] /= s;
            }

            tl = -999999;
            tu = 999999;

            s = 0;
            ys = 0;
            for (i=1;i<=l;i++) {
                if (d[i]<0) {
                    tu = min(tu,-y[i]/d[i]);
                }
                if (d[i]>0) {
                    tl = max(tl,-y[i]/d[i]);
                }
                s += d[i];
                ys += y[i];
                if (s<0) {
                    tl = max(tl,(R[i]-ys)/s);
                    tu = min(tu,(L[i]-ys)/s);
                }
                if (s>0) {
                    tu = min(tu,(R[i]-ys)/s);
                    tl =  max(tl,(L[i]-ys)/s);
                }
            }

            for (i=l+1;i<=n;i++) {
                if (d[i]<0) {
                    tu = min(tu,-y[i]/d[i]);
                }
                if (d[i]>0) {
                    tl = max(tl,-y[i]/d[i]);
                }
                if (d[i]<0) {
                    tl = max(tl,(1-y[i])/d[i]);
                }
                if (d[i]>0) {
                    tu = min(tu,(1-y[i])/d[i]);
                }
            }
            tl = unique(tl,tu,th);
            for (i=1;i<=n;i++) y[i] += tl * d[i]; 
        } while (isContainedDKLR(y,L,R,n,l) == false);

    }

    void SamplefromCpKLR_plane(double* yD,double* yC,double* L,double* R,int n,unsigned int* th) {
        int i;
        double s,t = 0,u;
        SamplefromDKLR(yD,L,R,n-1,n-1,th);
        u = R[n];
        for (i=1;i<=n-1;i++) {
            t = unique(0,2*3.1415926535897932,th);
            s = sqrt(yD[i]);
            u -= yD[i];
            yC[2*i] = s * cos(t);
            yC[2*i-1] = s * sin(t);
        }
        s = sqrt(u);
        yC[2*i] = s * cos(t);
        yC[2*i-1] = s * sin(t);
    }

    bool checkR(double* a,double* R,int n) {
        double s;
        int i;
        s = a[1]*a[1];
        if (s>R[1]) return false;
        for (i=2;i<=n;i++) {
            s += a[i]*a[i];
            if (s>R[i]) return false;
        }
        return true;
    }

    double Ratio_prob(double* pf,int n) {
        initialize();
        int N = SampleMult * max(10000,n*n); 
        int i,k,k2;
        double* y,*yc;

        
        long long int ktotal=0;
        long long int ktotal2=0;
        int endflag=0;
        
        int numthread = omp_get_max_threads(); 

        //cout << "numthreads=" << numthread << endl;
        //cout << "N=" << N << endl;
        //numthread = 1;
        
        #pragma omp parallel  num_threads(numthread) private(i,y,yc,k,k2) 
        {
            //cout << omp_get_thread_num() << " start" << endl;
            
            srand(int(time(NULL)) ^ omp_get_thread_num());
            y = sample_share_double[4*omp_get_thread_num()];    
            yc = sample_share_double[4*omp_get_thread_num()+1];

            //cout << omp_get_thread_num() << " alloc" << endl;

            unsigned int* th = rstates + 3*omp_get_thread_num();
            th[1] += rand();
            th[2] ^= (rand() << 12);
            i = 0;
            double* eR = sample_share_double[4*omp_get_thread_num()+2];
            double* eL = sample_share_double[4*omp_get_thread_num()+3];
            for (i=1;i<=n/2;i++) eR[i] = pf[2*i];
            for (i=1;i<=n/2;i++) eL[i] = 0;

            //Initial point
            for (i=1;i<=n/2;i++) y[i] = eR[i] / 100.0;

            //Shuffle
            i = 0;
            while ((i<100*n) && (endflag==0)) {
                    SamplefromDKLR(y,eL,eR,n/2,n/2,th);
                    i++;
                    //if (i%100==0) cout << omp_get_thread_num() << " loop"  << i << endl;
            }
            k=0;
            k2 = 0;
            while (endflag==0) {
                SamplefromCpKLR_plane(y,yc,eL,eR,n/2,th);
                if (checkR(yc,pf,n-1)==true) k++;
                k2++;
                if ((k2 & 0x100) !=0) {
                #pragma omp critical
                    {
                        ktotal2 += k2;
                        k2 = 0;
                        ktotal += k;
                        k = 0;
                        if (ktotal2>N) endflag = 1;
                        //cout << omp_get_thread_num() << " " << ktotal << " " << ktotal2 << endl;
                    }                    
                }
            }
            #pragma omp atomic  
            ktotal += k;
            #pragma omp atomic  
            ktotal2 += k2;
            //cout << "rplast_block" << omp_get_thread_num() << " " << 1.0 * ktotal / ktotal2 << endl;
        } // end of omp parallel
        //cout << "rplast_ret=" << 1.0 * ktotal / ktotal2 << endl;
        return 1.0 * ktotal / ktotal2; 
    }
    
    double Ratio_prob_odd(double* pf,int n) {
        //assume n is odd
        initialize();
        int N = SampleMult * max(10000,2*n*n); 
        int i,k,k2;
        double* y,*yc;
        
        long long int ktotal=0;
        long long int ktotal2=0;
        int endflag=0;
        int error=0;
/*        
        cout << "pf: ";
        for (i=0;i<=n+2;i++) cout << pf[i] << " ";
        cout << endl;
  */      
        int numthread = omp_get_max_threads(); 
        //numthread = 1;
        int j;
/*        
        int* a1 = new int[numthread];
        int* a2= new int[numthread];
        for (i=0;i<numthread;i++) {
            a1[i] = 0;
            a2[i] = 0;
        }
  */      
        #pragma omp parallel  num_threads(numthread) private(i,j,y,yc,k,k2) 
        {
            srand(int(time(NULL)) ^ omp_get_thread_num());
            y = sample_share_double[4*omp_get_thread_num()];    
            yc = sample_share_double[4*omp_get_thread_num()+1];

            unsigned int* th = rstates + 3*omp_get_thread_num();
            th[1] += rand();
            th[2] ^= (rand() << 12);

            i = 0;
            double* eR = sample_share_double[4*omp_get_thread_num()+2];
            double* eL = sample_share_double[4*omp_get_thread_num()+3];
            for (i=1;i<=(n+1)/2;i++) eR[i] = pf[2*i];
            for (i=1;i<=(n+1)/2;i++) eL[i] = 0;

            //Initial point
            for (i=1;i<=(n+1)/2;i++) y[i] = eR[i] / 100.0;
            

            
            //Shuffle
            i = 0;
            while ((i<10*n) && (endflag==0)) {
                SamplefromDKLR(y,eL,eR,(n+1)/2,(n+1)/2,th);
                i++;
                //for (j=1;j<=n/2;j++) cout << y[j] << " ";
                //cout << endl;
            }
            k=0;
            k2 = 0;
            while (endflag==0) {
                for (j=0;j<1;j++) SamplefromCpKLR_plane(y,yc,eL,eR,(n+1)/2,th);
                for (j=1;j<n;j++) {
                    yc[j] /= sqrt(1.0-yc[(n+1)]*yc[(n+1)]);
                }
                if (std::isnan(yc[1])==true) {
                    error=1;
                    endflag=1;
                }
                if (checkR(yc,pf,n-1)==true) k++;
                k2++;
                if ((k2 & 0x100) !=0) {
                #pragma omp critical
                    {
                        //cout << omp_get_thread_num() << "/" << numthread <<  " " << k << " " << k2 << endl;
                        //a1[omp_get_thread_num()] += k;
                        //a2[omp_get_thread_num()] += k2;
                        ktotal2 += k2;
                        k2 = 0;
                        ktotal += k;
                        k = 0;
                        if (ktotal2>N) endflag = 1;
                        //cout << ktotal << " " << ktotal2 << endl;
                    }                    
                }

            }
            #pragma omp atomic  
            ktotal += k;
            #pragma omp atomic  
            ktotal2 += k2;
        }
        //cout << error << " " << ktotal << " / " << ktotal2 << endl;
/*
        double ratio=0;
        for (i=0;i<numthread;i++) {
            if (a2[i]>0) {
                cout << "thread-" <<  i << " " << a1[i] << " " << a2[i] << " " << (double)a1[i] / a2[i] << endl;
                ratio += (double)a1[i] / a2[i];
            }
        }
        ratio /= (double)numthread;
        cout << "ratio=" << ratio << " ratioB=" <<  1.0 * ktotal / ktotal2 << endl;
  */      
        if (error==0) {
                return 1.0 * ktotal / ktotal2;
        } else {
                return Ratio_prob_odd(pf,n);
        }
    }

    
    inline void Samplefromsphere(double* y,int n,double radius=1.0,unsigned int* th=0) {
        //returns y[1...n]
        if (th==0) {
            initialize();
            th = rstates;
        }
        
        int i;

        double* sum = y+n+1;
        *sum = 0;
        for (i=1;i<=n;i++) {
            y[i] = normal_dist(0,1,th);
            *sum += y[i]*y[i];
        }
        *sum=sqrt(*sum);
        for (i=1;i<=n;i++) {
            y[i] *= radius / (*sum) ;
            
        }

    }
    inline void Samplefromsphere(double* y,int n,unsigned int* th=0) {
        return Samplefromsphere(y,n,1.0,th);
    }
    
    double PruningProbDirect(double* R,int n) {

        sampling_tools::initialize();
        int N = 1 * n*n*omp_get_max_threads(); 
        //cout << "ppd: N=" << N << endl;
        
        int k,k2;
        double* y;

        long long int ktotal=0;
        long long int ktotal2=0;

        int eflag;
        int numthread = omp_get_max_threads(); 

        #pragma omp parallel  num_threads(numthread) private(y,k,k2) 
        {
            srand(int(time(NULL)) ^ omp_get_thread_num());
            y = sample_share_double[4*omp_get_thread_num()];    
            unsigned int* th = rstates + 3*omp_get_thread_num();
            th[1] += rand();
            th[2] ^= (rand() << 12);

            k=0;
            eflag = 0;
            k2 = 0;
            while (eflag==0) {
                    Samplefromsphere(y,n,1.0,th);
                    if (checkR(y,R,n-1)==true) k++;
                    k2++;
                    if ((k2 & 0x100) !=0) {
                        #pragma omp critical
                        {
                            ktotal2 += k2;
                            k2 = 0;
                            ktotal += k;
                            k = 0;
                            if (ktotal2>N) eflag = 1;

                        }                    
                    }
            }
            #pragma omp atomic  
            ktotal += k;
            #pragma omp atomic  
            ktotal2 += k2;
        }
        return 1.0*ktotal/ktotal2;
    }
    
};

#endif
