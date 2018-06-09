#ifndef _inc_latticetools_cpp
#define _inc_latticetools_cpp

#define INPUT_SQUARED 0x40
#define INPUT_NONSQUARED 0x20


//return Euclidean length in double
double LengthOf(vec_ZZ& b) {
    ZZ ip;
    double r;
    InnerProduct(ip,b,b);
    conv(r,ip);
    return sqrt(r);
}

RR HalfGamma(int n) {
    //Gamma(n/2)
    int i;
    RR ret;
    ret = to_RR(1);
    if (n%2==0) {
        for (i=1;i<n/2;i++) ret *= i;
    } else {
        ret = sqrt(3.141592653589793238462);
        for (i=1;i<n;i+=2) {
            ret *= 0.5 * i;
        }
    }
    return ret;
}

//Returns the volume of k-ball of given radius
RR VolumeBall(int k,double radius) {
    RR ret;
    ret = exp(to_RR(k*log(radius) + 0.5*k*log(3.1415926535897932)));
    ret /= HalfGamma(k+2);
    return ret;
}

template <typename T> void GramSchmidt(mat_ZZ& L,T** mu,T* rlen,int n,int ioffset=0,int joffset=0) {

    //indexes are slided by ioffset and joffset

    int i,j,k;
    T** r;
    T** mur;
    r = new T*[n];
    for (i = 0;i<n;i++) r[i] = new T[n];

    mur = new T*[n];
    for (i = 0;i<n;i++) mur[i] = new T[n];

    ZZ ip;

    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) {
            InnerProduct(ip,L[i],L[j]);
            conv(r[i][j],ip);   //r[i][j] = L[i]*L[j];
            for (k=0;k<j;k++) {
                r[i][j] -= r[j][k] * mur[i][k];
            }
            if (i>j) mur[i][j] = r[i][j] / r[j][j];
        }
    }
    for (j=0;j<n;j++) {
        conv(rlen[j+ioffset],r[j][j]);    // = |b*_j|^2
    }
    
    for (i=0;i<n;i++) {
        for (j=0;j<i;j++) {
            conv(mu[i+ioffset][j+joffset],mur[i][j]);
        }
        mu[i+ioffset][i+joffset] = 1;
    }

    for (i = 0;i<n;i++) delete [] r[i];
    delete [] r;
    
    for (i = 0;i<n;i++) delete [] mur[i];
    delete [] mur;
}



namespace lattice_tools {
    
    int init=0;
    
    RR* vol_unit_ball;
    double* gaussian_heuristic_constant;
    double* smallghconst;
    double* simhkzconst;

    quad_float** mu;
    quad_float* c;
    
    void initialize() {
        if (init!=0) return;

        int oldprec = RR::precision();
        RR::SetPrecision(256);
        int i,j;
        vol_unit_ball = new RR[latticemaxdim+1];
        gaussian_heuristic_constant = new double[latticemaxdim+1];
        smallghconst = new double[51];
        simhkzconst = new double[51];
        for (i=1;i<=latticemaxdim;i++) {
            gaussian_heuristic_constant[i] = to_double( exp(log(VolumeBall(i,1.0))*(-1.0/i)));
            vol_unit_ball[i] = VolumeBall(i,1.0);
        }
        RR::SetPrecision(oldprec);
        
        mu = new quad_float*[latticemaxdim+1];
        for (i=0;i<=latticemaxdim;i++) mu[i] = new quad_float[latticemaxdim];
        c = new quad_float[latticemaxdim+1];
        

        double CNconstant[] = {0.593208,0.582161,0.561454,0.544344,0.522066,0.502545,0.479628,
    0.459819,0.438675,0.413708,0.392483,0.370717,0.344447,0.322574,
    0.297318,0.273761,0.249247,0.225483,0.199940,0.173832,0.147417,
    0.123425,0.100035,0.074487,0.043089,0.020321,-0.013844,-0.042863,
    -0.068204,-0.093892,-0.124345,-0.151097,-0.183912,-0.214122,
    -0.241654,-0.274612,-0.302966,-0.330965,-0.367514,-0.391956,
    -0.426507,-0.457813,-0.488113,-0.518525,-0.554184,-0.585479,
    -0.617705,-0.646749,-0.671864,-0.687300};    

        for (i=0;i<50;i++) {
            //cout << "CN[" << i << "]=" << CNconstant[i] << endl;
        }
        
        for (i=0;i<25;i++) swap(CNconstant[i],CNconstant[49-i]);
        for (i=0;i<=49;i++) {
            CNconstant[i] = exp(CNconstant[i]);
        }
        double det;
        for (i=1;i<=49;i++) {
            det = 0;
            for (j=0;j<=i;j++) det += log(CNconstant[j]);
            det = exp(1.0 * det/(i+1)) * lattice_tools::gaussian_heuristic_constant[i+1];      //gaussian heuristic
            smallghconst[i+1] = CNconstant[i] / det;
            //cout << "sgc[" << i+1 << "]=" << smallghconst[i+1] << endl; 
        }

        for (i=0;i<25;i++) swap(CNconstant[i],CNconstant[49-i]);
        det=0;
        for (j=0;j<50;j++) det += log(CNconstant[j]);
        det = exp(det / 1.0 / 50);   //det^(1/n)
        for (i=1;i<=50;i++) simhkzconst[i] = CNconstant[i-1] / det;

        init = 1;
    }
    
    void finish() {
        delete [] vol_unit_ball;
        delete [] gaussian_heuristic_constant;
        delete [] c;
        for (int i=0;i<=latticemaxdim;i++) delete [] mu[i];
        delete [] mu;
        init = 0;
    }
    
    
    template <typename T> double LatticeVolumeRoot(T* c,int dim,char opt) {
        initialize();
        //gaussian heuristic radius
        //computing from c[1...dim]
        T det;
        conv(det,0);
        int i;
        for (i=1;i<=dim;i++) {
                det += log(c[i]);
         } //det = log(det)
        if (opt==INPUT_SQUARED) det *= 0.5;
        
         det /= dim;
         det = exp(det);
         return to_double(det);
    }

    template <typename T> double LatticeGH(T* c,int dim,char opt=INPUT_NONSQUARED) {
        initialize();
        //gaussian heuristic radius
        //computing from c[1...dim]
        T det;
        conv(det,0);
        int i;
        for (i=1;i<=dim;i++) {
                det += log(c[i]);
         } //det = log(det)
        if (opt==INPUT_SQUARED) det *= 0.5;
        
         det /= dim;
         det = exp(det) * gaussian_heuristic_constant[dim];    //gaussian heuristic
         return to_double(det);
    }

    template <typename T> double LatticeApprox(T* c,int dim,char opt) {
        initialize();
        //approximation factor
        //computing from c[1...dim]
        return to_double(c[1]) / LatticeGH(c,dim,opt);
    }

    double LatticeApprox(mat_ZZ& L) {
        initialize();
        GramSchmidt(L,mu,c,L.NumRows());
        return LatticeApprox(c-1,L.NumRows(),INPUT_SQUARED);
    }
    
    double LatticeGH(vec_RR& c,int jj,int dim,char opt) {
        initialize();
        //gaussian heuristic radius
        //computing from c[1...dim]
        double det;
        conv(det,0);
        int i;
        for (i=jj;i<jj+dim;i++) {
                det += to_double(log(c(i)));
         } //det = log(det)
        if (opt==INPUT_SQUARED) det *= 0.5;
        
         det /= dim;
         det = exp(det) * gaussian_heuristic_constant[dim];    //gaussian heuristic
         return to_double(det);
    }
    double LatticeVolumeRoot(vec_RR& c,int jj,int dim,char opt) {
        initialize();
        //computing from c[1...dim]
        double det;
        conv(det,0);
        int i;
        for (i=jj;i<jj+dim;i++) {
                det += to_double(log(c(i)));
         } //det = log(det)
        if (opt==INPUT_SQUARED) det *= 0.5;
        
         det /= dim;
         det = exp(det);
         return to_double(det);
    }

    template <typename T> double LatticeVolume(T* c,int dim,char opt) {
        initialize();
        //gaussian heuristic radius
        //computing from c[1...dim]
        T det;
        conv(det,0);
        int i;
        for (i=1;i<=dim;i++) {
                det += log(c[i]);
         } //det = log(det)
         if (opt==INPUT_SQUARED) det *= 0.5;
        
         //det /= dim;
         det = exp(det);
         return to_double(det);
    }
    template <typename T> void Normalize(T* c,int dim,char opt=0) {
        initialize();
        //normalize c[1...n] s.t. (prod of c[i]) =1
        T det;
        conv(det,0);
        int i;
        for (i=1;i<=dim;i++) {
                det += log(c[i]);
         } //det = log(det)
         det /= dim;
         det = exp(det);
        for (i=1;i<=dim;i++) c[i] /= det;
    
    }

    double LatticeGH(double* c,int dim) { 
        initialize();
        return LatticeGH(c,dim,INPUT_NONSQUARED);
    }
    
    double LatticeGH(mat_ZZ& L) {
        initialize();
        GramSchmidt(L,mu,c,L.NumRows());
        return LatticeGH(c-1,L.NumRows(),INPUT_SQUARED);
    }
/*
    double LatticeFEC(mat_ZZ& L) {
        initialize();
        GramSchmidt(L,mu,c,L.NumRows());

        RR fec = FullENUMCost(c-1,L.NumRows(),1.0,INPUT_SQUARED);
        return to_double(log(fec));
    }
  */
    double LatticeVolumeRoot(mat_ZZ& L) {
        initialize();
        GramSchmidt(L,mu,c,L.NumRows());
        return LatticeVolumeRoot(c-1,L.NumRows(),INPUT_SQUARED);
    }

    double LatticeVolume(mat_ZZ L) {
        initialize();
        GramSchmidt(L,mu,c,L.NumRows());
        return LatticeVolume(c-1,L.NumRows(),INPUT_SQUARED);
    }

    quad_float LatticeHaflVolume(mat_ZZ& L) {
        initialize();
        GramSchmidt(L,mu,c,L.NumCols());
        quad_float ret;
        int i;
        int n = L.NumRows();

        ret = 1;
        for (i=0;i<n/2;i++) {
            ret *= sqrt(c[i]);
        }
        for (i=n/2;i<n;i++) ret /= sqrt(c[i]);
        return ret;
    }
    

    
    
    template <typename T> void alpha_decompose(T* dip,vec_ZZ& d,mat_ZZ& li,quad_float **mu,int offset=0) {

        //d[0...n-1]: vector
        //li[0...n-1][0..n-1]: lattice
        //mu[0...n-1][0..n--1]: GS coeff
        //compute dip[i]=<v,b*i> 
        //This is equal to: decomposing v=\sum alpha[i].b*[i] <=> dip[i]=alpha[i].|b*i|^2
        int i,k;
        int n = li.NumRows();
        for (i=0;i<n;i++) {
            dip[i] = 0;
            for (k=0;k<n;k++) dip[i] += conv<T>(d[k] * li[i][k]);        //<d[j].b[i]>
            for (k=0;k<i;k++) dip[i] -= conv<T>(mu[i+offset][k+offset]) * dip[k];     //<d[j],b*[i]>
        }
    }

    template <typename T> void alpha_decompose(T* alpha,vec_ZZ& v,mat_ZZ& L,int opt=0) {
        initialize();
        GramSchmidt(L,mu,c,L.NumRows());
        alpha_decompose(alpha,v,L,mu);
        if (opt==1) {
            for (int i=0;i<L.NumCols();i++) {
                //to adjust v=\sum alpha[i].b*[i] <=> dip[i]=alpha[i]
                alpha[i] /= conv<T>(c[i]);    
            }
        }
        if (opt==2) {
            for (int i=0;i<L.NumCols();i++) {
                //to adjust v=\sum alpha[i].b*[i] <=> dip[i]=alpha[i]*|b*[i]|
                alpha[i] /= sqrt(conv<T>(c[i]));    
            }
        }
    }
    
    RR FullENUMCost(quad_float* c,int m,double radius,int opt=INPUT_SQUARED) {

        //Assume that the input is squared value

        int i;
        RR ret = to_RR(0);
        double gh=radius;

        //cout << "gh_org=" << gh << endl;

        RR t = to_RR(1);

        for (i=1;i<=m;i++) {
            if (opt==INPUT_SQUARED) {
                t *= to_RR(gh / sqrt(c[m-i+1]));
                ret += t * vol_unit_ball[i];
            } 
            if (opt==INPUT_NONSQUARED) {
                t *= to_RR(gh / c[m-i+1]);
                ret += t * vol_unit_ball[i];
            } 
            //cout << i << " " << ret << " " << c[m-i+1] << endl;
        }    

        //cout << "fec=" << ret << endl;

        return ret;
    }

    RR FullENUMCost(quad_float* c,int m,int opt=INPUT_SQUARED,double alpha=1.0) {

        double gh = alpha * LatticeGH(c,m,opt);
        return FullENUMCost(c,m,gh,opt);

    }

    RR FullENUMCost(vec_RR& c,int jj,int m,int opt=INPUT_SQUARED) {

        //Assume that the input is squared value

        int i;
        RR ret = to_RR(0);
        double gh = LatticeGH(c,jj,m,opt);

        //cout << "gh_org=" << gh << endl;

        RR t = to_RR(1);

        for (i=1;i<=m;i++) {
            if (opt==INPUT_SQUARED) {
                t *= to_RR(gh / sqrt(c(jj+m-i)));
                ret += t * vol_unit_ball[i];
            }
            if (opt==INPUT_NONSQUARED) {
                t *= to_RR(gh / c(jj+m-i));
                ret += t * vol_unit_ball[i];
            }
        }

        return ret;
    }

    RR FullENUMCost(mat_ZZ& L) {
        initialize();
        GramSchmidt(L,mu,c,L.NumCols());
        return FullENUMCost(c-1,L.NumCols(),INPUT_SQUARED);
    }

    void LatticeGSLengths(mat_ZZ L,quad_float* c) {
        //output is c[1..n]
        initialize();
        GramSchmidt(L,mu,c+1,L.NumRows());
    }
    
    
}



#endif
