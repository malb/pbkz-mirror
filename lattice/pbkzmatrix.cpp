#ifndef _inc_pbkz_matrix_define
#define _inc_pbkz_matrix_define

//vector of small elements

template <typename T> class pbkzvector :  public std::vector<T> {
    private:
    public:

    int length() { return this->size(); };
};

template <typename T> inline void addvector(pbkzvector<T>& a,int c,pbkzvector<T>& b) {
    int n = b.size();
    for (int i=0;i<n;i++) a[i] += b[i] * c;
}

template <typename T,typename T2> inline void addvector(pbkzvector<T>& a,T2 c,pbkzvector<T>& b) {
    int n = b.size();
    for (int i=0;i<n;i++) a[i] += b[i] * c;
}

template <typename T> inline T squared_norm(pbkzvector<T>& a) {
    T ret = 0;
    for (int i=0;i<a.size();i++) {
        ret += a[i]*a[i];
    }
    return ret;
}

template <typename T> inline T LengthOf(pbkzvector<T>& a) {
    T ret = squared_norm(a);
    ret = sqrt(ret);
    return ret;
}

template <typename T> std::ostream& operator <<(std::ostream &stdout, pbkzvector<T>  &arg) {

    stdout << "[";
    for (int i=0;i<arg.size();i++) {
        stdout << arg[i];
        if (i+1<arg.size()) stdout << " ";
    }
    stdout << "]";
    return stdout;
}

template <typename T> std::istream& operator >>(std::istream &stdin, pbkzvector<T>  &arg) {
    int left=0;
    int right=0;
    int i=0;
    std::string temp;
    std::string line;
    while (stdin) {
        stdin >> temp;
        //cout << "t1=" << temp << endl;
        while (temp.length()>0) {
            if (temp[0]=='[') {
                left++;
                temp = temp.substr(1,temp.length()-1);
            } else 
            if (temp[temp.length()-1]==']') {
                right++;
                temp = temp.substr(0,temp.length()-1);
            } else {
                line += temp + " ";
                temp="";
            }
        }
        if (left-right==1) {
            if (i>=arg.size()) arg.resize(i+1);
            if (right>=1) line = "[" + line + "]";
            line >> arg[i]; 
            i++;
            line = "";
        }
    }
    return stdin;
}

//Matrix with small elements
//This is a high-speed compatible class with mat_ZZ in NTL
template <typename T> class pbkzmatrix : public pbkzvector<pbkzvector<T> > {

    private:
        
    public:
    
    //For compatibility with NTL library    
    void SetDims(int rows,int cols) {
        this->resize(rows);
        for (int i=0;i<rows;i++) {
            (*this)[i].resize(cols);
        }
    }
    
    int NumRows() {
        return this->size();
    }

    int NumCols() {
        if (NumRows()==0) return 0;
        return (*this)[0].size();
    }

    pbkzvector<T>& operator [] (int i) {
        if ((i < 0) || (i >= this->size()))  {
            cout << "pbkzmatrix: range error " << i << endl;
            exit(0);
        }
        pbkzvector<T>* p = this->data();
        return p[i];
    }
    
};

template <typename T> void SaveLattice(pbkzmatrix<T>& L,std::string fname) {
    ofstream logstream;
    logstream.open(fname.c_str(),ios_base::trunc);
    logstream.precision(12);
    logstream << L;
    logstream.close();
}

template <typename T> void LoadLattice(pbkzmatrix<T>& L,std::string fname) {

    ifstream logstream;
    logstream.open(fname.c_str(),ios_base::in);
    try {
        logstream >> L;
    }
    catch (std::exception e) {
        L.SetDims(1,1);
    }
    logstream.close();
}

#endif