#ifndef _inc_compute_even_simplex
#define _inc_compute_even_simplex


template <typename T> T EvalAt(T* F,int n,double a) {
    //Assume a be a small number
    T ret,t;
    ret = 0;

    int i;
    t = 1;
    for (i=0;i<n;i++) {
        ret += t * F[i];
        t *= T(a);
    }
    return ret;
}

template <typename T> void Integ(T* ret, int n,T* F,double low) {    
    int i;
    for (i=0;i<n;i++) {
        ret[i+1] = F[i] / T(i+1.0);
    }
    ret[0] = -EvalAt(ret,n,low);
    return;
}
 
#define opt_volume 0x01
#define opt_volume_prob 0x02
#define opt_surface 0x03
#define opt_surface_prob 0x04

namespace EvenSimplex {

    int initialized=0;
    double*** F;
    quad_float*** QFF;
    double* es_share_double;
    int max_n;
    
    void Initialize () {
        if (initialized!=0) return;
        
        int n = 300;
        int i,j;
        F = new double**[n+2];
        for (i=0;i<=n;i++) {
            F[i] = new double*[n+2];
            for (j=0;j<n+2;j++) {
                F[i][j] = new double [n+2];
                //for (int k=0;k<n+2;k++) F[i][j][k] = 0;
            }
        }
        QFF = new quad_float**[n+2];
        for (i=0;i<=n;i++) {
            QFF[i] = new quad_float*[n+2];
            for (j=0;j<n+2;j++) {
                QFF[i][j] = new quad_float [n+2];
                //for (int k=0;k<n+2;k++) F[i][j][k] = 0;
            }
        }
        es_share_double = new double[n];
        
        initialized=1;
    }    

    //double precision does not return correct value for n>75
    double EvenSimplexVolumedouble(int n,double* R2) {
        int i,j;
        int k;
        Initialize();
        
        for (i=0;i<=n;i++) {
            for (j=0;j<=n+2;j++) {
                for (k=0;k<=n+2;k++) {
                    F[i][j][k] = 0;
                }
            }
        }

        
        F[1][0][1] = 1; //F[1,0]=y
        F[1][1][0] = R2[1]; //F[1,0]=R_1

        for (i=2;i<=n;i++) {
            Integ(F[i][0], i+1,F[i-1][0],to_double(0));
            for (j=1;j<=i;j++) {
                Integ(F[i][j],i,F[i-1][j],R2[j]);
                F[i][j][0] += EvalAt(F[i][j-1],i+1,R2[j]); 
            }
        }
        double pb = EvalAt(F[n-1][n-1],n+1,R2[n-1]);

        for (i=2;i<=n-1;i++) pb *= i;
        return pb;
    }
    double EvenSimplexVolumeqf(int n,double* R2,int option) {
        
        int i,j;
        int k;
        if (initialized==0) Initialize();
        
        if (n==1) return R2[1];
        
        for (i=0;i<=n;i++) {
            for (j=0;j<=n+2;j++) {
                for (k=0;k<=n+2;k++) {
                    QFF[i][j][k] = 0;
                }
            }
        }
        
        QFF[1][0][1] = 1; //F[1,0]=y
        QFF[1][1][0] = R2[1]; //F[1,0]=R_1

        for (i=2;i<=n;i++) {
            Integ(QFF[i][0], i+1,QFF[i-1][0],0);
            for (j=1;j<=i;j++) {
                Integ(QFF[i][j],i,QFF[i-1][j],R2[j]);
                QFF[i][j][0] += EvalAt(QFF[i][j-1],i+1,R2[j]); 
                //cout << "i=" << i << " j=" << j << " " << QFF[i][j][0] << endl;
            }
        }

        double ret = 0;
        switch(option){
        case opt_volume: return to_double(QFF[n][n][0]);
        case opt_surface:{
          quad_float pb = EvalAt(QFF[n-1][n-1], n + 1, R2[n-1]);
          ret = to_double(pb);
          break;
        }
        case opt_surface_prob:{
          quad_float pb = EvalAt(QFF[n-1][n-1],n+1,R2[n-1]);
          for(i = 2; i <= n - 1; ++i) pb *= i;
          ret = to_double(pb);
          break;
        }
        case opt_volume_prob:
          for(i = 2; i <= n; ++i) QFF[n][n][0] *= i;
          ret = to_double(QFF[n][n][0]);
        }
        return ret;
        
        // if (option==opt_volume) return to_double(QFF[n][n][0]);
        // //cout << QFF[n][n][0] << endl;

        // quad_float pb = EvalAt(QFF[n-1][n-1],n+1,R2[n-1]);
        // //cout << pb << endl;
        
        // if (option==opt_surface) return to_double(pb);
        // for (i=2;i<=n-1;i++) pb *= i;
        // //cout << pb << endl;
        // if (option==opt_surface_prob) return to_double(pb);

        // for (i=2;i<=n;i++) {
        //     //cout << "i=" << i << " "  << QFF[n][n][0] << endl;
        //     QFF[n][n][0] *= i;
        // }
        // //cout << QFF[n][n][0] << endl;
        // if (option==opt_volume_prob) return to_double(QFF[n][n][0]);
        
    }
    
    RR EvalAt(vec_RR& F,double a) {
        //Assume 0<=a<=1
        RR ret,t;
        ret = 0;

        int i,n;
        n = F.length();
        t = 1;
        for (i=0;i<n;i++) {
            ret += t * F[i];
            t *= to_RR(a);
        }
        return ret;;
    }

    vec_RR Integ(vec_RR& F,double low) {

        vec_RR ret;

        int i,n;
        n = F.length();
        ret.SetLength(n+1);
        ret[0] = 0;
        for (i=0;i<n;i++) {
            ret[i+1] = F[i] / to_RR(i+1);
        }
        ret[0] = -EvalAt(ret,low);

        return ret;
    }


    RR EvenSimplexVolume(int n,double* R2,int option) {

        //Computing the volume of n-dimensional trancated simplex
        //{ (x1,...,xn) : \sum[i=1,...,l] x_i < R_l for l=1,...,n}
        //F[n,n] is the volume
        //pb is the probability that Pr{x <- surface of full-simplex}[x is in trancated simplex]
        
        if (n<=60) {
            return to_RR(EvenSimplexVolumeqf(n,R2,option));
        }
        
        int bprec = RR::precision();
        RR::SetPrecision(max(2*n,120));
        
        int i,j;
        vec_RR** F;
        F = new vec_RR*[n+2];
        for (i=0;i<=n;i++) F[i] = new vec_RR[n+2];

        F[1][0].SetLength(2);
        F[1][0][1] = 1; //F[1,0]=y

        F[1][1].SetLength(2);
        F[1][1][0] = R2[1]; //F[1,0]=R_1

        for (i=2;i<=n;i++) {
            F[i][0] = Integ(F[i-1][0],0);
            for (j=1;j<=i;j++) {
                F[i][j] = Integ(F[i-1][j],R2[j]);
                F[i][j][0] += EvalAt(F[i][j-1],R2[j]); 
                //cout << "F[" << i << "," << j << "]=";
                //for (int k=i-j;k>=0;k--) cout << F[i][j][k] << "*x^" << k << " ";
                //cout << endl;
                if (i==j) {
                    RR v;
                    v = F[i][i][0];
                    for (int k=2;k<=i;k++) v *= k;
                    //cout << "frac_check(" << i << "): " << F[i][i] << " ... " << v << endl;
                }
            }
            
        }
    /*    
        RR pb = EvalAt(F[n-1][n-1],to_RR(R2[n-1]));
        for (i=2;i<=n-1;i++) pb *= i;

        cout << pb << endl;
        cout << F[n][n] << endl;
    */
        RR ret;

        if ((option==opt_volume) || (option==opt_volume_prob)) {
            ret = F[n][n][0];
        }
        //cout << "ret=" << ret << endl;

        if ((option==opt_surface) || (option==opt_surface_prob)) {
            ret = EvalAt(F[n-1][n-1],R2[n-1]);
        }
        //cout << ret << endl;

        if (option==opt_surface_prob) {
            for (i=2;i<=n-1;i++) ret *= i;
        }
        //cout << ret << endl;

        if (option==opt_volume_prob) {
            for (i=2;i<=n;i++) ret *= i;
            //cout << "last_check: " << ret << endl;
        }
        //cout << ret << endl;

        for (i=0;i<=n;i++) delete [] F[i];
        delete [] F;
        RR::SetPrecision(bprec);
        return ret;
    }
    
}



#endif
