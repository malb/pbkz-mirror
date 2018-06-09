#ifndef _inc_enum_wrapper_cpp
#define _inc_enum_wrapper_cpp

//#include "bkzlibwrapper.hpp"

void coeff_to_zz(vec_ZZ& v,std::vector<int>& coeffs,mat_ZZ& L,int offset) {

    v = coeffs[0] * L[offset];
    int i = 0;
    for (int i=1;i<coeffs.size();i++) {
        if (coeffs[i]!=0) {
            v += coeffs[i] * L[i+offset];
        }
    }
}


struct FoundVectorInfo {
    double norm = -1;
    int projection = -1;
    std::vector<int> coeffs;
    
    void display() {
        cout << "[" << norm << " " << " " << projection << ": ";
        for (int i=0;i<coeffs.size();i++) cout << coeffs[i] << " ";
        cout << "]" << endl;
    }
    
};

bool operator < (FoundVectorInfo& a,FoundVectorInfo& b) {
    if (a.projection < b.projection) return true;
    if (a.projection == b.projection) {
        if (a.norm < b.norm) return true;
    }
    return false;
}

void clear(FoundVectorInfo& a) {
    a.norm = -1;
    a.projection = -1;
    a.coeffs.resize(0);
    
}

struct FoundENUMVectors {
    int maxhold = 10000;
    std::vector<FoundVectorInfo> data;
    int foundnum;
    long long int totalnodes;
    double etime;
    
    void clearmemory() {
        data.clear();
        data.resize(maxhold);
        foundnum = 0;
        totalnodes = 0;
    }
    void clearmemory(int hold) {
        maxhold = hold;
        clearmemory();
    }
    
    int countnonzero() {
        int j = 0;
        for (int i=0;i<data.size();i++) {
            if (data[i].norm>0) j++;
        }        
        return j;
    }

    int getminproj() {
        int j = -1;
        for (int i=0;i<data.size();i++) {
            if (data[i].norm>0) {
                if (j==-1) {
                    j = data[i].projection;
                }
                j = min(j,data[i].projection);
            }
        }        
        return j;
    }
    
    int shiftproj() {
        int j = getminproj();
        for (int i=0;i<data.size();i++) {
            if (data[i].norm>0) {
                data[i].projection -=j;
            }
        }
        return j;
    }

    int getemptyslot() {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm==-1) return i;
        }
        return -1;
    }

    int getminslot() {
        int ms = 0;
        for (int i=0;i<data.size();i++) {
            if (data[i].norm==-1) return i;
            if (data[ms] < data[i]) ms = i;
        }
        return ms;
    }
    
    void storeelement(FoundVectorInfo& newFV) {
        int ms = getemptyslot();
        if (ms<0) {
            ms = getminslot();
        }
        
        if (( data[ms].norm==-1 ) ||
             (data[ms].projection > newFV.projection) ||
             ((data[ms].projection == newFV.projection) && ( data[ms].norm > newFV.norm )) ) {
                //replace
            data[ms].norm = newFV.norm;
            data[ms].projection = newFV.projection;
            data[ms].coeffs.resize(newFV.coeffs.size());
            for (int i=0;i<newFV.coeffs.size();i++) data[ms].coeffs[i] = newFV.coeffs[i];
        }
    }
    
    void display() {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm>0) {
                cout << "slot=" << i << " proj=" << data[i].projection << " norm=" << data[i].norm << "(" << sqrt(data[i].norm) << "): ";
                for (int j=0;j<data[i].coeffs.size();j++) cout << data[i].coeffs[j] << " ";
                cout << endl;
            }
        }        
    }

    template <typename T> mat_ZZ get_zz_rep(LatticeBasis<T>& B,int offset=1) {
        int j = countnonzero();
        mat_ZZ ret;
        ret.SetDims(j,B.L.NumCols());
        j = 0;
        for (int i=0;i<data.size();i++) {
            if (data[i].norm > 0) {
                coeff_to_zz(ret[j],data[i].coeffs,B.L,data[i].projection-1);
                j++;
            }
        }
        return ret;
    }
    
    void sort() {
        std::sort(data.begin(),data.end());
    }

    void simplify(FoundVectorInfo& a) {
        //this breaks norm information
        int frontzero=0;
        int backzero=a.coeffs.size()-1;

        while ((backzero > 1) && (a.coeffs[backzero]==0)) backzero--;
        a.coeffs.resize(backzero+1);

        while ((frontzero < a.coeffs.size()) && (a.coeffs[frontzero]==0)) frontzero++;
        for (int i=0;i<a.coeffs.size()-frontzero;i++) a.coeffs[i] = a.coeffs[i+frontzero];
        a.coeffs.resize(a.coeffs.size()-frontzero);
        a.projection += frontzero;
    }

    void cutbylowerbound(bkzfloat radius) {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm <= radius) {
                data[i].projection = -1;
            }
        }
    }
    void duplicate() {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm > 0) {
                simplify(data[i]);
            }
        }
        
        for (int i=0;i<data.size();i++) {
            if (data[i].norm <= 0) {
                data[i].projection = -1;
            }
        }
        sort();
        for (int i=1;i<data.size();i++) {
            if ( data[i-1].projection== data[i].projection) {
                if (data[i-1].coeffs.size()==data[i].coeffs.size()) {
                    char eqflag=1;
                    for (int j=0;j<data[i].coeffs.size();j++) {
                        if (data[i].coeffs[j] != data[i-1].coeffs[j] ) {
                            eqflag = 0;
                            break;
                        }
                    }
                    if (eqflag == 1) {
                        clear(data[i-1]);
                    }
                }
            }
        }
    }
    void erasesingle() {
        for (int i=0;i<data.size();i++) {
            if (data[i].norm > 0) {
                int nonzerocount=0;
                for (int j=0;j<data[i].coeffs.size();j++) {
                    if (data[i].coeffs[j] != 0) nonzerocount++;
                }
                if (nonzerocount==1)  data[i].norm = -1;
            }
        }
    }
};


#include "vectorenumeration_boost.cpp"

template <typename T> void GSrepresentation(std::vector<T>& gsrept,vec_ZZ& vv,mat_ZZ& L,std::vector<std::vector<T> >& mu,std::vector<T>& cd) {
    int n = L.NumCols();
    int m = L.NumRows();
    gsrept.resize(m+1);
    
    //vv[0...n-1]: input vector
    //L[0...n-1][0..m-1]: lattice
    //mu[1...n][1..n]: GS coeff
    //cd[1...n]: |b*i|^2
    //compute gsrept[i]=<v,b*i> 
    //This is equal to: decomposing v=\sum alpha[i].b*[i] <=> dip[i]=alpha[i].|b*i|^2
    int i,k;
    for (i=1;i<=m;i++) {
        gsrept[i] = 0;
        for (k=0;k<n;k++) gsrept[i] += boost::lexical_cast<T>(vv[k] * L[i-1][k]);        //<d[j].b[i]>
        for (k=1;k<i;k++) gsrept[i] -= (T)(mu[i][k]) * gsrept[k];     //<d[j],b*[i]>
    }
    for (i=1;i<=m;i++) gsrept[i] /= cd[i];
}

template <typename T> void Coeffrepresentation(std::vector<T>& crept,vec_ZZ& vv,mat_ZZ& L,std::vector<std::vector<T> >& mu,std::vector<T>& cd) {
    GSrepresentation(crept,vv,L,mu,cd);
    int m = crept.size()-1;
    for (int i=m;i>=1;i--) {
        crept[i] = round(crept[i]);
        for (int j=i-1;j>=1;j--) {
            crept[j] = crept[j] - crept[i] * mu[i][j];
        }
    }
}

template <typename T> void GSrepresentation(std::vector<T>& gsrept,vec_ZZ& vv,LatticeBasis<T>& B) {
    B.updateGSBasis();
    GSrepresentation(gsrept,vv,B.L,B.gs.mu,B.gs.cd);
}

template <typename T> void Coeffrepresentation(std::vector<T>& gsrept,vec_ZZ& vv,LatticeBasis<T>& B) {
    B.updateGSBasis();
    Coeffrepresentation(gsrept,vv,B.L,B.gs.mu,B.gs.cd);
}

template <typename T> mat_ZZ ENUMbackend(LatticeBasis<T>& B,vec_ZZ& v,PruningFunction& PF,bkzfloat radius,int enumoption,int vl=0,std::string stringoptions="",int modeoption=modesvp) {
    FoundENUMVectors EV; 
    int slide = 0; 
    
    int elim = 10000;
    
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

    bkzfloat lowerradius=0;
    if (options["lowerbound"]!="") {
        lowerradius = boost::lexical_cast<bkzfloat>(options["lowerbound"]);
    }

    if (options["vl"]!="") {
        vl = atoi(options["vl"].c_str());
    }
    
    int parallel = 1;
    if (options["parallel"]!="") {
        parallel = atoi(options["parallel"].c_str());
    }

    bool sizereduceflag = false;
    if (options["sizereduce"]=="true") {
        sizereduceflag=true;
    }
    
    B.updateGSBasis();


    int finishmode = finish_mode_nowaste;
    if (options["finishmode"]=="exact") {
        finishmode = finish_mode_exact;
    }
    if ((enumoption & enum_mode_all_vectors) != 0) {
        finishmode = finish_mode_exact;
    }
    //end of parameter 

    
    if (modeoption==modesvp) {

        if (parallel<=1) {
            if (iend-istart+1 < 150) {
                //lattice_enum::ENUMCoreboost<double,double,T>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
                lattice_enum::ENUMCoreboost_double<T>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
            } else {
                lattice_enum::ENUMCoreboost<long double,long double,T>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption);
            }
        } else {
            if (iend-istart+1 < 150) {
                lattice_enum::ENUMCoreboostmt_double<T>(EV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
            } else {
                lattice_enum::ENUMCoreboostmt<long double,long double,T>(EV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,vl,enumoption,finishmode);
            }
        }
    }     

    if (modeoption==modecvp) {
        //Todo: getgspresentation
        std::vector<T> gsrept;
        GSrepresentation(gsrept,v,B);

        if (parallel<=1) {
            if (iend-istart+1 < 150) {
                lattice_enum::ENUMCVCoreboost_double<T>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,gsrept,vl,enumoption);
            } else {
                lattice_enum::ENUMCVCoreboost<long double,long double,T>(EV,slide,istart,iend-istart+1,elim,radius*radius,PF.rd,B.gs.mu,B.gs.cd,gsrept,vl,enumoption);
            }
        } else {
            if (iend-istart+1 < 150) {
                lattice_enum::ENUMCVCoreboostmt_double<T>(EV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,gsrept,vl,enumoption,finishmode);
            } else {
                lattice_enum::ENUMCVCoreboostmt<long double,long double,T>(EV,slide,istart,iend-istart+1,elim,parallel,radius*radius,PF.rd,B.gs.mu,B.gs.cd,gsrept,vl,enumoption,finishmode);
            }
        }
    }     

    if (vl>=3) EV.display();
    EV.duplicate();
    if (lowerradius>0) EV.cutbylowerbound(lowerradius);

    mat_ZZ ret;
    ret = EV.get_zz_rep(B);
    
    if ((istart>1) && (sizereduceflag==true)) {
        //Size-reduction
        LatticeBasis<T> B2;
        B2.resize(istart+1,B.L.NumCols());
        for (int i=0;i<istart;i++) {
            B2[i] = B[i];
        }
        B2.updateGSBasis(1,istart-1);
        for (int i=0;i<ret.NumRows();i++) {
            B2[istart] = ret[i];
            B2.GScomputed = istart-1;
            SizeReduce(B2,0,istart);
            ret[i] = B2[istart];
        }
    }
    return ret;
}

template <typename T> mat_ZZ ENUM(LatticeBasis<T>& B,bkzfloat radius,bkzfloat tprob,int enumoption=0,int pfoption=0,int vl=0,std::string otheroptions="") {
    PruningFunction PF;
    bkzfloat prob = pruning_func::SetPruningFunction(B,PF, radius , tprob,pfoption,vl);
    vec_ZZ zerovec;
    return ENUMbackend<T>(B,zerovec,PF,radius,enumoption,vl,otheroptions,modesvp);
}

template <typename T> mat_ZZ ENUMCV(LatticeBasis<T>& B,vec_ZZ& targetv,bkzfloat radius,bkzfloat tprob,int enumoption=0,int pfoption=0,int vl=0,std::string otheroptions="") {
    PruningFunction PF;
    bkzfloat prob = pruning_func::SetPruningFunction(B,PF, radius , tprob,pfoption);
    return ENUMbackend<T>(B,targetv,PF,radius,enumoption,vl,otheroptions,modecvp);
}


#endif