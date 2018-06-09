#ifndef _inc_latticemisc_cpp
#define _inc_latticemisc_cpp

template <typename T,typename T2> void vec_copy(std::vector<T>& a,std::vector<T2>& b) {
    a.resize(b.size());
    for (int i=0;i<b.size();i++) a[i] = b[i];
}

template <typename T,typename T2> void vec_copy(std::vector<std::vector<T> >& a,std::vector<std::vector<T2> >& b) {
    a.resize(b.size());
    for (int i=0;i<b.size();i++) vec_copy(a[i],b[i]);
}

int to_int(double b) {return (int)b;};

template <typename T> int to_int(T b) {return boost::lexical_cast<int>(b);};

void erasepoint(std::string& xstr) {
    for (int i=0;i<xstr.length();i++) {
        if (xstr[i]=='.') {
            xstr = xstr.substr(0,i);
            return;
        }
    }
    return;
}

template <typename T> std::string to_intstring(T& a) {
    std::string xstr = a.str(0,std::ios_base::fixed);
    erasepoint(xstr);
    return xstr;
}

template <> std::string to_intstring(double& d) {
    float10 a = d;
    std::string xstr = a.str(0,std::ios_base::fixed);
    erasepoint(xstr);
    return xstr;
}

template <> std::string to_intstring(long double& d) {
    float10 a = d;
    std::string xstr = a.str(0,std::ios_base::fixed);
    erasepoint(xstr);
    return xstr;
}

void RowTransformwrap(vec_ZZ& A, vec_ZZ& B, const ZZ& MU1) {
    //A=A-MU*B
    if (MU1==0) return;
    if (MU1==1) {
        A -= B;
        return;
    }
    if (MU1==-1) {
        A += B;
        return;
    }
    A -= B*MU1;
    return;
    
}

bool conv(mpfr_float& a,ZZ& b) {    a = (mpfr_float)to_stdstring(b); }

bool conv(float50& a,ZZ& b) {    a = (float50)to_stdstring(b); }
bool conv(float40& a,ZZ& b) {    a = (float40)to_stdstring(b); }
bool conv(float30& a,ZZ& b) {    a = (float30)to_stdstring(b); }
bool conv(float20& a,ZZ& b) {    a = (float20)to_stdstring(b); }
bool conv(float15& a,ZZ& b) {    a = (float15)to_stdstring(b); }
bool conv(float10& a,ZZ& b) {    a = (float10)to_stdstring(b); }


bool conv(ZZ& a,bkzfloat& b) {
    conv(a,to_stdstring(b).c_str());
}

bool conv(quad_float& a,ZZ& b) {
    a = to_quad_float(b);
}

bool conv(double& a,ZZ& b) { a = to_double(b); }
bool conv(long double& a,ZZ& b) { a = boost::lexical_cast<long double>(b); }



template <typename T> void copyvec(std::vector<T>& a,std::vector<T>& b) {
    a.resize(b.size());
    for (int i=0;i<b.size();i++) a[i] = b[i];
}

template <typename T> inline void addvector(T* a,int c,T*b,int n) {
    for (int i=0;i<n;i++) a[i] += b[i] * c;
}

template <typename T> inline void display(T* a,int n) {
    cout << "[";
    for (int i=0;i<n;i++) cout << a[i] << " ";
    cout << "]" << endl;
}

template <typename T> void conv_to_ZZ(ZZ& a,T& b) {
    if (abs(b) < 2e+9) {
        conv(a,(long int)b);
        return;
    }
    
    //Todo: better method
    double logdet = boost::lexical_cast<double>(log(abs(b)));
    int prec = logdet / 2.302585092;
    prec += 10;
    

#pragma omp critical
 {
    int precisionback = mpfr_float::default_precision() ;
    mpfr_float::default_precision(prec);
    mpfr_float detmp = static_cast<mpfr_float>(b);
    detmp = round(detmp);
    std::string t = detmp.convert_to<std::string>();
    a = to_ZZ(t.c_str());
    mpfr_float::default_precision(precisionback);
 }
    return;
}




#endif