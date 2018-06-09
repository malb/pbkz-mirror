
#include "bkzlibwrapper.hpp"

#ifndef _bkzpreprocessboostmain_common
#define  _bkzpreprocessboostmain_common
    
//Heuristic preprocessing subroutine to find a basis of low enumeration cost
template <typename PFLOAT,typename CFLOAT> bkzfloat ProgressiveBKZPreprocessCore(LatticeBasis<CFLOAT>&B,mat_ZZ* U,int istart,int iend,bkzfloat costlimit,bkzfloat alphalimit,double ualpha,bkzfloat uprob,int vl=0) {

    bkzfloat current_cost_lb;
    //Greedy BKZ
    int endflag = 0;
    int bstart;
    int bsize;
    int lowest_bsize = 20;
    bsize = lowest_bsize;
    int bmax = bsize;
    
    int loop=0;
    char direction = +1;    //-1=backward +1=forward

    if (direction == 1) bstart = istart;
    if (direction == -1) bstart =  istart + (iend-istart+1-bsize);
    
    int SRcomputed = 0;

    double radius = ualpha * lattice_tools::LatticeGH(B.gs.c,istart,iend);
    std::vector<CFLOAT> lbt;
    lbt.resize(iend+1);
    ENUMCostLB_maketable(lbt,istart,iend,uprob);

    
    while (endflag == 0) {
        //search suitable (bstart,bend) s.t. large |b*i| and small cost
        int localbs = min(bsize,iend-bstart+1);
        
        bkzfloat ghratio = B.gs.c[bstart]/lattice_tools::LatticeGH(B.gs.c,bstart,bstart+localbs-1);
        char executeflag = 0;
        bkzfloat blockcost;
        if (ghratio > 1.0 + 1.0/localbs) {
            //check cost
             blockcost =  bkzconstants::target_cost_lb_basis(B.gs.c,bstart,bstart+localbs-1);
             if (blockcost < 2 * bkzconstants::target_cost_lb(localbs)) {
                executeflag = 1;
            }
        }
        
        if (bsize==lowest_bsize) executeflag = 1;

        if (executeflag!=0) {
            double alpha,r;
            bkzfloat prob;
            PBKZParam(alpha,r,prob,localbs);
            int result;
            for (int i=SRcomputed+1;i<=bstart+localbs-1;i++) {
                SizeReduce(B,U,i);
            }
            double radius = alpha * (double)lattice_tools::LatticeGH(B.gs.c,bstart,bstart+localbs-1);
            result = EnumthenUpdate<CFLOAT>(B,U,bstart,bstart+localbs-1,radius,prob,0,vl,"");
            SRcomputed = bstart-1;    //Output of EnumthenUpdate do not guarantee size reduction

            if (result!=0) {
                if (vl>=3) {
                    cout << "check bstart=" << bstart << " bsize=" << localbs << "           \r";
                    cout.flush();
                }
                bmax = max(bmax,localbs);
                if (++loop%100==0) {
                    int swap;
                    local_LLL(B,U,0.999,istart,iend,stnormal,VL0,swap);
                }

            }
        }
        bstart += direction;
        
        int trigger_index = iend-bsize+2;
        if (direction == -1) trigger_index = istart-1;

        if (bstart==trigger_index) {
            
            current_cost_lb = ENUMCostLBwithtable(B.gs.c,lbt,istart,iend,radius); 
            
            bkzfloat alpha = B.gs.c[istart] / lattice_tools::LatticeGH(B.gs.c,istart,iend);

            if (vl>=1) {
                cout << "preprocess block=(" << istart << "," << iend << ") cost_lb=" << current_cost_lb << "/" << costlimit;
                if (alphalimit>0) cout << " alpha=" << alpha << "/" << alphalimit << " ";
                cout << "  \r";
                if (vl>=5) cout << endl;
                cout.flush();
            }
            if (current_cost_lb < costlimit) break;
            if ((alphalimit >0) && (alpha < alphalimit)) break;
            
            bsize--;
            if (bsize<lowest_bsize) {
                int shift = 10;
                if (bmax > 50) shift = 5;
                bsize = min(bmax+shift,(iend-istart+1)-20);
                if (vl>=3) {
                    B.gs.displayc(istart,iend);
                }
            }
            
            if (direction == 1) bstart = istart;
            if (direction == -1) bstart =  istart + (iend-istart+1-bsize);
        }
    }
    return current_cost_lb;
}


#endif



#ifdef __bkzpreprocessbig
template <typename PFLOAT,typename CFLOAT> void ExtractSubBlock(LatticeBasis<CFLOAT>&  LB,LatticeBasis<PFLOAT>& B,int istart,int iend) {
#endif

#ifdef __bkzpreprocesssmall
template <typename PFLOAT,typename CFLOAT> void ExtractSubBlock(SmallLatticeBasis<CFLOAT,CFLOAT>&  LB,LatticeBasis<PFLOAT>& B,int istart,int iend) {
#endif
    
    //Extracting sublattice B[istart..iend] to LB using Gram-Schmidt data
    B.updateGSBasis(1,iend);
    
    //Find factor
    PFLOAT mings = B.gs.c[istart];
    for (int i=istart+1;i<=iend;i++) {
        mings = min(mings,B.gs.c[i]);
    }    
    mings = mings/10000;
    int dim = iend - istart + 1;
    
    //diagonal
    for (int i=0;i<dim;i++) {
        for (int j=0;j<dim;j++) {
            PFLOAT ll = B.gs.mu[i+istart][j+istart] * B.gs.c[j+istart]/mings;
#ifdef __bkzpreprocesssmall
            LB.L[i][j] = (CFLOAT)ll;
#endif
#ifdef __bkzpreprocessbig
            //need to convert NTL_ZZ
            conv_to_ZZ(LB.L[i][j] ,ll);
#endif
        }
        PFLOAT ll = B.gs.c[i+istart]/mings;
#ifdef __bkzpreprocesssmall
        LB.L[i][i] = (CFLOAT)ll;
#endif
#ifdef __bkzpreprocessbig
            conv_to_ZZ(LB.L[i][i] ,ll);
#endif
    }
}   //End of ExtractSubBlock


#ifdef __bkzpreprocessbig
    void SetIdentity(mat_ZZ&  U) {
#endif

#ifdef __bkzpreprocesssmall
    template <typename T> void SetIdentity(pbkzmatrix<T>& U) {
#endif
        for (int i=0;i<U.NumCols();i++) {
            for (int j=0;j<U.NumRows();j++) {
                if (i==j) {
                    U[i][j] = 1;
                } else {
                    U[i][j] = 0;
                }
            }
        }
}
    
    
    