#ifndef _inc_pbkzmain_cpp
#define _inc_pbkzmain_cpp

#include "pbkzproperty.cpp"



namespace progressive_bkz {


#define insert_all 0x01
#define insert_gauss 0x02

    int BasisUpdate(mat_ZZ& B,mat_ZZ* U,quad_float** mu,quad_float* c,quad_float* b,quad_float** B1,int jj,int bs,int slide,int** foundvec,double* foundlen,int fnum,int m,int strategy) {

        int i,j,k;
        int kk = jj + bs - 1;
        int n = B.NumCols();
        quad_float t1;
        quad_float *tp;
        quad_float delta;
        delta = to_quad_float(0.99);
        long quit = 0;

        ZZ MU;

    #ifdef debug_clength_check
        checkapprox(B,B1);
    #endif

        //Erase zero vectors
        for (i=0;i<fnum;i++) {
            int ck=0;
            for (j=jj;j<=kk;j++) {
                if (foundvec[i][j]!=0) {
                    ck=1;
                    break;
                }
            }
            if (ck==0) {
                for (k=0;k<=kk;k++) foundvec[i][k] = foundvec[fnum-1][k];
                foundlen[i] = foundlen[fnum-1];
                i--;
                fnum--;
            }
        }    

        //Erase the duplicates
        for (i=0;i<fnum;i++) {
            int ck=0;
            int k;
            for (j=i+1;j<fnum;j++) {
                k = jj;
                ck = 0;
                while (k<=kk) {
                    if (foundvec[i][k]!=foundvec[j][k]) {
                        ck++;
                        break;
                    }
                    k++;
                }
                if (ck==0) {
                    foundvec[i][0] = min(foundvec[i][0],foundvec[j][0]);
                    foundlen[i] = max(foundlen[i],foundlen[j]);
                    for (k=0;k<=kk;k++) foundvec[j][k] = foundvec[fnum-1][k];
                    foundlen[j] = foundlen[fnum-1];
                    j--;
                    fnum--;
                }
            }
        }

        //counting non-zero vectors
        int fct=0;
        for (int fi=0;fi<fnum;fi++) {
           if (foundlen[fi]>0) {
                   fct++;
           }
        }



        if ((fct>0) && (strategy==insert_all)) {
            //sorting by (index,length)
            //Rewrite to the marge sort
            for (i=0;i<fct;i++) {
                for (j=i+1;j<fct;j++) {
                  if (((foundvec[i][0] > foundvec[j][0]) || (foundvec[i][0] == foundvec[j][0])) && (foundlen[i] > foundlen[j])) {
                       //swap
                       swap(foundvec[i],foundvec[j]);
                       swap(foundlen[i] , foundlen[j]);                      
                    }
                }
            }

            jj+= slide;
            //make zero buffer from jj,...,jj+fct-1
            for (i = m; i >= jj; i--) {
            // swap i, i+fct
            swap(B(i+fct), B(i));
            if (U) swap((*U)(i+fct), (*U)(i));
            swap(B1[i+fct],B1[i]);
            swap( b[i+fct], b[i]);
            }
            for (j = jj; j < jj+fct; j++) {
                for (i = 1; i <= n; i++) conv(B(j, i), 0);
                if (U) {
                    for (i = 1; i <= m; i++)
                        conv((*U)(j, i), 0);
                }
            }

            //Storing found vectors to the buffer
            for (int fi=0;fi<fct;fi++) {
                int* uvtemp = foundvec[fi];
                for (i = jj; i <= kk; i++) {
                    if (uvtemp[i] == 0) continue;
                    conv(MU, uvtemp[i]);
                    RowTransform2(B(jj+fi), B(i+fct), MU);
                    if (U) RowTransform2((*U)(jj+fi), (*U)(i+fct), MU);
                }
                for (i = 1; i <= n; i++) conv(B1[jj+fi][i], B(jj+fi, i));
                b[jj+fi] = InnerProduct(B1[jj+fi], B1[jj+fi], n);   
            }


            //Sort by the LLL
            verbose = 0;
            //cout << endl;
            //cout << "B=[" << endl;
            //cout << B << endl;
            ll_LLL_QP(B, U, delta, 0, 0, B1, mu, b, c, kk+fct, jj-slide, quit);

            //Erasing zero vectors
            for (i = kk+1; i <= m; i++) {
                // swap i, i+fct
                swap(B(i+fct), B(i));
                if (U) swap((*U)(i+fct), (*U)(i));
                tp = B1[i+fct]; B1[i+fct] = B1[i]; B1[i] = tp;
                t1 = b[i+fct]; b[i+fct] = b[i]; b[i] = t1;
            }
            jj-=slide;
        }

        //compute c for next index
        ll_LLL_QP(B, U, delta, 0, 0,B1, mu, b, c, min(kk+1,m),min(kk,m), quit);

    #ifdef debug_clength_check
        for (i=1;i<=min(kk+1,m);i++) {
            if (c[i]==0) {
                for (i=1;i<=m;i++) {
                    cout << "c[" << i << "]=" << c[i] << endl;
                }
                exit(0);
            }
        }
    #endif    
        return fct;
    }



    int BKZCore(mat_ZZ& B,mat_ZZ* U,quad_float** mu,quad_float* c,quad_float* b,quad_float** B1,int& lcindex,int jj,int bs,int fnum,double probs,double radius,char mode,int elim,int m,int parallel,int enumopt,int vl,double delta,bool optimizepf=false) {

        sharememalloc();

        int kk = min(jj + bs - 1,m);
        int slide;
        int findex; //# found vectors
        double* foundlen = bkz_share_double[8];
        int** foundvec = bkz_share_foundvec;

    #ifdef debug_clength_check
        quad_float* cback = BKZsharememory::bkz_share_quad_float2[0];   //mempry for volume factor
        for (int i=1;i<jj+bs;i++) cback[i] = c[i];
    #endif

        if (lcindex < kk) {
            //If LLL basis was not computed in enough index
            long quit = 0;
            if (vl>=0) cout << "Compute LLL[" << lcindex << "," << kk << "]" << endl;
            ll_LLL_QP(B, U, to_quad_float(0.99), 0, 0,B1, mu, b, c, kk,lcindex, quit);
        }
        //Setting pruning radius and function
        double clim = InitRadius(c,jj,kk-jj+1,radius,mode);
        double* normal_pruning_func = bkz_share_double[6];
        set_approximate_extreme_pruning(normal_pruning_func,c+jj,kk-jj+1,probs,0);
        if ((optimizepf==true) && (kk-jj+1>=40)) {
            pruning_func::optimize_pruning_function(normal_pruning_func-1,sqrt(clim),c+jj-1,kk-jj+1,probs,vl);
        }

        if (parallel==1) {
            //outputs are foundvec and foundlen
            lattice_enum::ENUMCore(foundvec,foundlen,fnum,findex,slide,jj,kk-jj+1,elim,clim,normal_pruning_func,mu,c,vl,enumopt );
        } else {
            lattice_enum::ENUMCoremt(foundvec,foundlen,fnum,findex,slide,jj,kk-jj+1,elim,clim,normal_pruning_func,mu,c,vl,parallel,enumopt,finish_mode_nowaste);
        }
        //Update the basis
        lcindex = min(kk+2,m);
        int ret = BasisUpdate(B,U,mu,c,b,B1,jj,bs,slide,foundvec,foundlen,fnum,m,insert_all);

    #ifdef debug_clength_check
        int i = 1;
        while (i<jj+bs) {
            if (cback[i] > c[i]) break;
            if (cback[i] > c[i]) {
                cout << "|b*i| check error" << endl;
                break;
            }
            i++;
        }
    #endif
        return ret;
    }


    void BKZpreprocessOPT(mat_ZZ& B,mat_ZZ* U,quad_float** mu,quad_float* c,quad_float* b,quad_float** B1,int& lcindex,int jj,int bs,int fnum,double probs,double radius,double enumpersec, int limitsec,int m,int parallel,int vl,int rl) {
        //Memo: vl=verbose level
        //      rl=recursive_level

        quad_float **minB1 = progressive_bkz::bkz_share_quad_float3[4+rl*2];
        quad_float **minmu = progressive_bkz::bkz_share_quad_float3[5+rl*2];
        quad_float *minc = progressive_bkz::bkz_share_quad_float2[11+rl*2];
        quad_float *minb = bkz_share_quad_float2[12+rl*2];

        int i,j;
        int n = B.NumCols();
        mat_ZZ minB;    //Backup of basis that achieves minumum cost
        mat_ZZ minU;    //Backup of Unimodular that achieves minumum cost
        //Backup data
        for (i=0;i<m;i++) minc[i] = 0;
        if (U) minU = *U;
        minB = B;
        for (i=0;i<=m;i++) {
            for (j=0;j<=n;j++) minB1[i][j] = B1[i][j];
            for (j=0;j<i;j++) minmu[i][j] = mu[i][j];
            minb[i] = b[i];
            minc[i] = c[i];
        }

        quad_float mincost;
        quad_float totalmincost;

        quad_float ccost;
        quad_float ucost;
        quad_float firstcost;

        int subbs;
        //For computing costs
        firstcost = mincost = ucost = 2*ENUMCost(c,jj,bs,probs,radius,'0');
        if (vl>=3) cout << "Preprocess init cost=" << mincost << " (" << mincost / enumpersec << "sec)" << endl;

        double r,alpha,prob;
        int localupdate;
        //initial values
        subbs = 20;

        int subtour = 0;
        double logfec;
        double targetlogfec;
        double startcputime = clock();
        double pftime;

        double totalminsec=to_double(mincost / enumpersec);
        double costsec;

        while (1) {
            while(1) {
                //Need to update blocksuze?
                logfec = to_double(log(lattice_tools::FullENUMCost(c+jj-1,bs,INPUT_SQUARED)));
                targetlogfec = ComputeLogTargetFEC(bs,subbs);
                if ((logfec < targetlogfec) && (subtour>5)) {
                    subbs++;
                    subtour = 0;
                } else {
                    break;
                }
                break;
            };


            //Compute params
            r= compute_target_r(subbs);
            alpha = 1.0 / lattice_tools::gaussian_heuristic_constant[subbs] * exp(-0.25*(subbs-1.0)*log(r)) * ((1.0+subbs)/subbs);
            if (alpha<1.00) alpha = 1.0;
            prob = 2.0 * exp(-1.0*subbs*log(alpha));
            if (prob>0.5) prob=0.5;
            subtour++;
            //old param setting
            //alpha = 1.05;
            //prob = exp(-1.0*subbs*log(1.05));

            //applyBKZ
            localupdate = 0;
            for (i=jj;i<=jj+bs-1;i++) {
                int localbs = min(subbs,jj+bs-i);
                if (BKZCore(B,U,mu,c,b,B1,lcindex,i,localbs,fnum,prob,alpha,'G',10,m,1,enum_mode_find_shortest,vl-2,0)!=0) {
                    localupdate=1;
    #ifdef _inc_vectorpool_cpp
                    if (recp!=0) {
                        for (j=jj-1;j<jj+localbs;j++) {
                            Vectorpool::Addvector(*recp,B[j]);
                        }
                    }
    #endif
                }
            }
            double cputime = (clock() - startcputime) / CLOCKS_PER_SEC; //elapsed CPU time

            //Backup?
            pftime = clock();
            ucost = 2*ENUMCost(c,jj,bs,probs,radius,'0');        
            startcputime += (clock()-pftime); 
            if (mincost > ucost) {
                mincost = ucost;
                minB = B;
                for (i=0;i<=m;i++) {
                    if (minB1[i]!=B1[i]) {
                        for (j=0;j<=n;j++) minB1[i][j] = B1[i][j];
                        for (j=0;j<i;j++) minmu[i][j] = mu[i][j];
                        minb[i] = b[i];
                        minc[i] = c[i];
                        if(U) minU[i] = (*U)[i];
                    }
                }
                totalminsec = min(totalminsec,to_double(mincost / enumpersec) + cputime);
            }

            //Break the loop?
            costsec = cputime + to_double(mincost / enumpersec);
            if (localupdate==0) {
                subbs++;
                subtour = 0;
            }


            //outputinfo
            if (vl>=3) {
                cout << "preprocess[" << jj << ":" << jj+bs-1 << "] bs=" << subbs << " FEC=" << logfec << "/tg" << targetlogfec;
                cout << " cost[sec]=" << ucost/enumpersec << "/" << mincost/enumpersec << "/" << cputime;
                if (limitsec <~0) {
                    cout << "(R" <<  1.05 * totalminsec - costsec << ")\r";
                } else {
                    cout << "(R" <<  limitsec - cputime << ")\r";
                }
                cout.flush();
            }
            if (costsec > 1.05 * totalminsec) break;
            if (cputime >= 3600 * 10) break;
            if (limitsec >0) {
                if (cputime >= limitsec) break;
            }
        }
        if (vl>=3) cout << endl;

        //Readfrombackup
        B = minB;
        for (i=0;i<=m;i++) {
            if (minB1[i]!=B1[i]) {
                for (j=0;j<=n;j++) B1[i][j] = minB1[i][j];
                for (j=0;j<i;j++) mu[i][j] = minmu[i][j];
                b[i] = minb[i];
                c[i] = minc[i];
                if (U) (*U)[i] = minU[i];
            }
        }
        if (vl>=3) {
            cputime = clock();
            cputime = (cputime - startcputime) / CLOCKS_PER_SEC;
            cout << "end preprocess at level" << rl << " first=" << firstcost / enumpersec << " current=" << mincost/enumpersec << " cputime=" << cputime << "   " << endl;
        }
        //cout << "mincost=" << mincost << endl;
    }

    long BKZmain(mat_ZZ& BB, mat_ZZ* UU,BKZproperty& BP) {

        RR delta = to_RR(0.99); // delta of LLL parameter
        double bdelta = 0;      //delta of Hermite factor
        int* beta = BP.beta;
        int* enumlim = BP.enumlim;

        long m = BB.NumRows();
        long n = BB.NumCols();
        int isapproxcomputed = 0;

        int nolargeextract=0;
        if (m<400) nolargeextract = 1;
        if (BP.endindex>0) {
            if (BP.endindex < 400) nolargeextract = 1;
        }
        

        //option heuristic break
        double prevlogfec=-1;

        NumSwaps = 0;
        int lllverboseback = verbose; 
        verbose =0;

        long i, j;
        RR ucost,prevcost;

        int maxhold = BP.holdvecs;
        mat_ZZ B,B_temp;
        B = BB;
        B.SetDims(m+1+maxhold, n);

        mat_RR B1;
        mat_RR mu;
        vec_RR c; // squared lengths of Gramm-Schmidt basis vectors
        vec_RR b; // squared lengths of basis vectors
        if (nolargeextract==0) {
            B1.SetDims(n+maxhold,n+1);
            mu.SetDims(n+1,n+1);
            c.SetLength(n+maxhold);
            b.SetLength(n+maxhold);
        }    

        quad_float **qB1 = bkz_share_quad_float3[0];
        quad_float **qmu = bkz_share_quad_float3[1];
        quad_float *qc = bkz_share_quad_float2[1]; // squared lengths of Gramm-Schmidt basis vectors
        quad_float *qb = bkz_share_quad_float2[2]; // squared lengths of basis vectors

        //wang

        //Translation unitary matrix
        mat_ZZ Ulocal;
        mat_ZZ *U;
        if (UU) {
           Ulocal.SetDims(m+1+maxhold, m);
           for (i = 1; i <= m; i++)
             for (j = 1; j <= m; j++)
              conv(Ulocal(i, j), 0);

           for (i = 1; i <= m; i++)
              conv(Ulocal(i, i), 1);
           U = &Ulocal;
        }
        else
           U = 0;

        long quit;
        long jj=0, kk;
        int lcindex;    //LLL, |b*i|, and mu[i,j] of RR are computed till this index 
        int qlcindex;   //LLL, |b*i|, and mu[i,j] of "quad_float approximation" are computed till this index
        if (nolargeextract==0) {
            //useRR
            compute_approxRR(B,m,n,B1,mu,b,c,U);
            for (i=0;i<=n;i++) {
                //cout << "c[" << i << "]=" << sqrt(c[i]) << endl;
            }

        } else {
            compute_approx(B,m,n,qB1,qmu,qb,qc,U);
            qlcindex = m;
        }
        lcindex = m;

        int nofirst=0;
        init_red_fudge();
        RR tour_max_cost = to_RR(-1);
        BP.return_tour_max_cost=0;
	 double heuristic_mult = 1.0;

        jj = BP.startindex-1;
        int update=1;
        double probs = BP.pruning_prob;
        int parallel=1;
        double entiregh;
        double entiredet; 
        double startbkztime=clock();
        //total_enum_nodes = 0;       //shared variable

        if (nolargeextract==0) {
            entiregh = lattice_tools::LatticeGH(c,1,m,INPUT_SQUARED);
            entiredet = lattice_tools::LatticeVolumeRoot(c,1,m,INPUT_SQUARED);    //det^(1/n)
        } else {
            entiregh = lattice_tools::LatticeGH(qc,m,INPUT_SQUARED);
            entiredet =  lattice_tools::LatticeVolumeRoot(qc,m,INPUT_SQUARED);    //det^(1/n)
        }

        int istart,iend;
        istart = max(1,BP.startindex);
        iend = BP.endindex;
        if (istart >= m) istart = m;
        if (iend >= m) iend = m;
        if (iend<0) iend = m;
        if (istart > iend ) swap(istart,iend);

        
        //For finding short vector
        RR min_allenum_cost = to_RR(-1);
        mat_ZZ minB = B;


        //For logoutput
        ofstream logstream;
        if (BP.tlname!="") { 
            logstream.open(BP.tlname.c_str(),ios_base::app);
        }


        RR sim_min_cost = minimum_cost(m,1.0);   

        //quad_float approximation of RR values
        int qastart = -1;
        int qaend = -1; //store quad_float approx. of [qastart,...,qaend] to qB1,qmu,qc,qb
        int qamax = 250;    //max blocksize of quad_float approx 
        int qasize = 0;
        mat_ZZ localB;
        mat_ZZ localU;
        mat_ZZ localB2;
        mat_ZZ localB3;
        if (nolargeextract==0) {
            localB2.SetDims(qamax+50,m);
        } else {
            qastart = 1;
            qaend = m;
            qamax = m;
            qasize = m;
        }
        RR factor;

        //initial params
        jj = istart-1;
        nofirst = 1;    //number of processing tour
        BP.processed_tours = 0;
    #ifdef _inc_vectorpool_cpp
        //For recvec
        if (BP.recvec==true) {
            BP.recvecdata.clear();
            for (i=0;i<m;i++) {
                Vectorpool::Addvector(BP.recvecdata,B[i]);
            }
        }
    #endif    
        //start of tour
        while (1) { 


            //setting start and end index
            jj++;   //index
            kk = min(jj+beta[0]-1, iend);  //process the block [jj...jj+kk-1]
            if (jj-1==BP.breakindex) break; 

            //Need to update quad_float approx?
            if ((nolargeextract==0) && ((jj < qastart) || (qaend < min(kk,iend)) )) {
                if (BP.verboselevel>=2) cout << "Update quad_float approx block [" << qastart << ":" << qaend << "]" << endl; 
                int qflag = 0;
                if (qastart!=-1) {
                    qflag = 1;
                    //Update the previous data
                    //Update B and U using localU

                    for (i=0;i<qasize;i++) {
                        for (j=0;j<qasize;j++) {
                            if (j==0) {
                                localB2[i] = localU[i][j] * B[j+qastart-1];
                            } else {
                                localB2[i] += localU[i][j] * B[j+qastart-1];
                            }
                        }
                    }
                    for (i=0;i<qasize;i++) {
                        B[i+qastart-1] = localB2[i];
                    }
                    for (i = qastart; i <=qaend; i++) {
                        for (j = 1; j <= n; j++) {
                            conv(B1(i,j), B(i, j));
                        }
                    }
                    for (i = qastart; i <= qaend; i++) {
                        InnerProduct(b(i),B1(i), B1(i));
                    }

                    ll_LLL_RR(B, U, delta, 0, 0, B1, mu, b, c, qaend, qastart, quit);  //must not swap
                    lcindex = qaend;
                }

                //Make local qf approx
                qastart = max(1,jj-1);
                qaend = min(iend,qastart + qamax -1); 
                if (iend - qaend<30) qaend = iend;

                qasize = qaend - qastart + 1;

                if ((qflag==1) && (qaend > lcindex)) {
                    for (i = lcindex; i <=qaend; i++) {
                        for (j = 1; j <= n; j++) {
                            conv(B1(i,j), B(i, j));
                        }
                    }
                    for (i = lcindex; i <= qaend; i++) {
                        InnerProduct(b(i),B1(i), B1(i));
                    }
                    ll_LLL_RR(B, U, delta, 0, 0, B1, mu, b, c, qaend, lcindex, quit);
                    lcindex=qaend;
                }

                if (BP.verboselevel>=2) cout << "Make quad_float approx block [" << qastart << ":" << qaend << "]" << endl; 
                localB.SetDims(qasize+BP.holdvecs,qasize);
                localU.SetDims(qasize+BP.holdvecs,qasize);

                power2(factor,10+(int)(0.23*qasize));
                factor /= sqrt(c(qaend));
                isapproxcomputed = 1;
                for (i=0;i<qasize;i++) {
                    for (j=0;j<i;j++) {
                        localB[i][j] = to_ZZ(sqrt(c(qastart+j)) * mu(qastart+i,qastart+j)* factor) ;
                    }
                    localB[i][i] = to_ZZ(sqrt(c(qastart+i))* factor);
                    for (j=i+1;j<qasize;j++) {
                        localB[i][j] = 0;
                    }
                }
                for (i=qasize;i<qasize+BP.holdvecs;i++) {
                    for (j=0;j<qasize;j++) {
                        localB[i][j] = 0;
                    }
                }
                for (i=1;i<=qasize;i++) {
                    for (j=1;j<=qasize;j++) {
                        localU(i,j) = 0;
                        if (i==j) localU(i,j) = 1;
                    }
                }
                compute_approx(localB,qasize,qasize,qB1,qmu,qb,qc,&localU);
                qlcindex = qasize;
            }


            //Operations about logFEC: after processing last index <=> before processing first index
            double logfec=0;
            
            //Output information and logfile
            std::ostringstream info;


            if (BP.verboselevel>=1) {
                if (BP.verboselevel>=3) cout << "--------------------------" << endl;
                if ((jj==istart) || (BP.verboselevel>=2)) {

                    info << timestr() << " " << "BKZ-" << beta[0] << " Tour=" << nofirst << "/" << BP.tourlim;
                    info << " i=" << jj << " swap=" << NumSwaps; 
                    info << " prob=" << probs << " Time=" << ((double)clock() - first_init_time) / CLOCKS_PER_SEC;
                    info << " (pftime=" << total_pruningfunc_time <<  ")";
                    cout << info.str() << endl;

                    if (jj==istart) {
                        std::ostringstream info2;
                        if (nolargeextract==0) {
                                info2 << "|b1|="  << sqrt(c(1)) << "=" << sqrt(c(1))/entiregh << "GH=" << sqrt(c(1))/entiredet << "det^(1/n)";
                        } else {
                                info2 << "|b1|="  << sqrt(qc[1]) << "=" << sqrt(qc[1])/entiregh << "GH=" << sqrt(qc[1])/entiredet << "det^(1/n)";
                        }
                        cout << info2.str() << " ";

                        if (logfec==0) {
                            if (nolargeextract==1) {
                                logfec = to_double(log(lattice_tools::FullENUMCost(qc+istart-1,iend-istart+1)));
                            }  else {
                                logfec = to_double(log(lattice_tools::FullENUMCost(c,istart,iend-istart+1)));
                            } 
                        }
                        cout << "logfec=" << logfec << "(" << exp(logfec-log(sim_min_cost)) << " approx)" <<endl;
                    }
                }
            }

            if (BP.tlname!="") {
                //Output numerical table

                if ((jj==istart) || (BP.tllevel>=2)) {
                    if (logfec==0) {
                        if (nolargeextract==1) {
                            logfec = to_double(log(lattice_tools::FullENUMCost(qc+istart-1,iend-istart+1)));
                        }  else {
                            logfec = to_double(log(lattice_tools::FullENUMCost(c,istart,iend-istart+1)));
                        } 
                    }
                    double b1len;
                    if (nolargeextract==0) {
                        b1len = to_double(sqrt(c(1)));
                    } else {
                        b1len = to_double(sqrt(qc[1]));
                    }
                    logstream << nofirst << "\t";   // #tours
                    logstream << beta[0] << "\t";   //current blocksize (as input parameter)
                    logstream << jj << "\t";   //processing index
                    logstream << NumSwaps << "\t";  //# LLL swaps
                    logstream << ((double)clock() - first_init_time) / CLOCKS_PER_SEC << "\t";  //Total CPU times since progressive_bkz::sharememalloc() was called
                    logstream << total_pruningfunc_time << "\t";    //Total CPU time to generate pruning functions
                    logstream << b1len << "\t"; //|b1|
                    logstream << b1len/entiregh << "\t"; //|b1|/GH(L)
                    logstream << b1len/entiredet << "\t"; //|b1|/det(L)^{1/n}
                    logstream << exp(log(b1len/entiredet)*1.0/m)  << "\t"; // root Hermite factor
                    logstream << logfec << "\t"; //log(FullENUMCost)
                    logstream << BP.abortbyfec << "\t"; //Target log(FullENUMCost)
                    logstream <<  exp(logfec-log(sim_min_cost)) << "\t"; //FEC(L)/FEC(simulated-HKZ-basis)
                    logstream << endl;
                } 
                
            }
            
            if (jj==istart) { 
                if (nolargeextract==1) {
                    logfec = to_double(log(lattice_tools::FullENUMCost(qc+istart-1,iend-istart+1)));
                    if (BP.abortbyapprox>0) {
                        BP.current_approx_ratio = to_quad_float(lattice_tools::FullENUMCost(qc,m) / sim_min_cost); 
                    }
                } 
                if (nolargeextract==0) {
                        logfec = to_double(log(lattice_tools::FullENUMCost(c,istart,iend-istart+1)));
                } 
                if (BP.abortbyfec>0) {
                    if (BP.verboselevel>=1) cout << "logFEC=" << logfec << " target_logFEC=" << BP.abortbyfec << endl; ;
                    if (logfec< BP.abortbyfec) break;
                }

                if (nolargeextract==1) {
                    if (BP.abortbyapprox>0) {
                        if (BP.verboselevel>=1) cout << "entire-approx-ratio=" << BP.current_approx_ratio << endl;
                        if (BP.current_approx_ratio < BP.abortbyapprox) break;
                    }
                }

                if (prevlogfec >0) {
			if (logfec > prevlogfec) {
				heuristic_mult *= 2.0;
				//cout << "hm=" << heuristic_mult << endl;
                        	if (BP.abortbyfec>0) {
                			if (BP.heuristicbreak==true) {
                           			if (BP.verboselevel>=VL3) cout << "FEC increasing -> break" << endl;
                           			break;
					}
                    		}
                    }
                }
                prevlogfec = logfec;
                BP.processed_tours++;
            }

            //end of output information and logfile
            char supported=0;

            //need to extend blocksize?
            if ((nofirst==1) && (jj==1)) tour_max_cost = BP.extend_blocksize_initmax;

            if (((BP.extend_blocksize==true) && (tour_max_cost>=1) && (nofirst>1) && (supported==0)) || 
                    ((BP.extend_blocksize==true) && (BP.extend_blocksize_initmax>=0) && (nofirst>=1) && (supported==0)))
            {
                int kkorg = kk;
                for (i=jj-qastart+1;i<=qaend;i++) {
                    //cout << "c[" << i << "]=" << c[i] << endl;
                }
                while (kk-qastart+1<qlcindex) {
                    ucost = 2 * to_RR(ENUMCost(qc,jj-qastart+1,kk-jj+1,BP.pruning_prob,BP.init_radius,BP.init_mode));
                    if (ucost > tour_max_cost * BP.extend_blocksizemult  * heuristic_mult) {
                        break;
                    }
                    kk++;
                }
                kk = max(kk-1,kkorg);
                if (BP.verboselevel>=3) cout << "Reset blocksize: " << kkorg-jj+1 << "->" << kk-jj+1  << endl;
            }     // setting jj and kk is finished



            //modify probability and mod_init_radius
            double gmodprobs;
            double mod_init_radius;
            gmodprobs = BP.pruning_prob;
            mod_init_radius = BP.init_radius;
            char modprobflag=0;
            RR targetcost;
            targetcost = tour_max_cost * heuristic_mult;
            if (BP.modifyprobability==true) {
                if (kk-jj+1 < min(200,beta[0])) {
                    modprobflag=1;
                    targetcost = tour_max_cost;
                }
            }
            if (BP.blocktimelimit>0) {
                modprobflag=1;
                double persec = lattice_enum::enum_speed_bench(BP.numthreads) * 1000000.0;
                targetcost = BP.blocktimelimit * persec;
            }
            if ((modprobflag==1) && (supported==0) && (targetcost!=0)) {
                if (BP.verboselevel>=3) cout << "modprob: (p=" << BP.pruning_prob << ",alpha=" << BP.init_radius << ")";
                gmodprobs = BP.pruning_prob;
                double lfactor = 10.0;

                //changing log(prob)
                int whilecount=0;
                while (1) {

                    //if blocksize is less than 10, modify probability is not working
                    if (kk-jj+1 <= 10) break;

                    if (gmodprobs > 1.0) gmodprobs = 1.0;
                    double talpha;
                    if (BP.blocktimelimit>0) {
                        talpha = BP.init_radius;
                    } else {
                        talpha = max(1.0,exp(1.0 / (kk-jj+1) * log(2.0/gmodprobs)));
                    }
                    ucost = 2 * to_RR(ENUMCost(qc,jj-qastart+1,kk-jj+1,gmodprobs,talpha,BP.init_mode));
                    //if (kk-jj+1>100) cout << "prob=" << gmodprobs << "(" << talpha << ") ucost=" << ucost << "/" << targetcost << endl;

                    if (ucost < targetcost) {
                        if (gmodprobs > 0.2) {
                            gmodprobs = 0.2;
                            mod_init_radius = exp(1.0 / (kk-jj+1) * log(2.0/gmodprobs));
                            break;
                        }
                        if (0.5 * targetcost < ucost) {
                            mod_init_radius = talpha;
                            break;
                        }
                        //too small
                        gmodprobs *= lfactor; 
                    } else {
                        if ((ucost > 5.0 * targetcost) && (gmodprobs<0.1)) {
                            gmodprobs /= to_double(exp(0.5*log(ucost/targetcost))); 
                        } else {
                            gmodprobs /= lfactor; 
                            lfactor = lfactor/1.5;
                            gmodprobs *= lfactor;
                            if (lfactor < 1.2) lfactor = 1.2;
                        }
                    }
                    if (whilecount++>20) break;
                }
                if ((gmodprobs < 0.5*BP.pruning_prob) && (BP.blocktimelimit<=0)) {
                    gmodprobs = 0.5*BP.pruning_prob;
                    mod_init_radius = max(1.0,exp(1.0 / (kk-jj+1) * log(2.0/gmodprobs)));
                }  


                if (BP.verboselevel>=3) cout << " -> (" << gmodprobs << "," << mod_init_radius << ")" << endl;
            }   //end of modify probability



            //Need to process this block?
            char process_flag=1; 
            if (supported==1) process_flag=0;
            if (BP.ghskip==true) {
                double gh = lattice_tools::LatticeGH(qc+jj-1-qastart+1,kk-jj+1,INPUT_SQUARED);
                gh *= BP.ghskipcoeff;
                if (BP.verboselevel>=3) {
                    if (nolargeextract==0) {
                            cout << "|b*_i|=" << sqrt(c(jj)) << " a*GH=" << gh << " ";
                    } else {
                            cout << "|b*_i|=" << sqrt(qc[jj]) << " a*GH=" << gh << " ";
                    }
                }
                if (c(jj) < gh*gh) {
                    process_flag = 0;
                    if (BP.verboselevel>=3) cout << "skip this block" << endl;
                } else {
                    if (BP.verboselevel>=3) cout << endl;
                }
            } // end of "Need to process this block"

            //Decide whether preprocess is necessary or not?
            if ((BP.preprocess==true) && (process_flag==1)){
                char preprocess_flag = 0;
                if (BP.preprocess_at_least>0) {
                    double persec = lattice_enum::enum_speed_bench(BP.numthreads) * 1000000.0;
                    if (persec<=0) persec = 10000000;

                    ucost = to_RR(2 * ENUMCost(qc,jj-qastart+1,kk-jj+1,gmodprobs,mod_init_radius,BP.init_mode));
                    //Note: ENUMCost returns expected (upper bound of) # points in cylinder intersecion
                    //      The actual # processed nodes is approximately doubled 
                    if (BP.preprocess_at_least <= ucost/lattice_enum::enum_speed_bench(1)) {
                        preprocess_flag = 1;
                        if (BP.verboselevel>=3) cout << "Upper # nodes = " << ucost << " (do preprocess, init_sec=" << ucost / persec << ")" << endl;
                    } else {
                        if (BP.verboselevel>=3) cout << "Upper # nodes = " << ucost << " (skip preprocess)" << endl;
                    }
                }
                if (preprocess_flag == 1) {
                    //Do preprocess
                    double clim = InitRadius(qc,jj-qastart+1,kk-jj+1,mod_init_radius,BP.init_mode);
                    double persec = lattice_enum::enum_speed_bench(BP.numthreads) * 1000000;
                    if (BP.preprocess_strategy==1) {
                        if (nolargeextract==0) {
                            BKZpreprocessOPT(localB,&localU ,qmu,qc,qb,qB1,qlcindex,jj-qastart+1,kk-jj+1,BP.holdvecs,gmodprobs,sqrt(clim),persec,BP.preprocess_max_time,qasize,1,BP.verboselevel,0);
                        } else {
                            BKZpreprocessOPT(B,U ,qmu,qc,qb,qB1,qlcindex,jj,kk-jj+1,BP.holdvecs,gmodprobs,sqrt(clim),persec,BP.preprocess_max_time,qasize,1,BP.verboselevel,0);
                        }
                    } 
                    #ifdef __develop
                    //under developing
                    if (BP.preprocess_strategy==2) {
                        if (nolargeextract==0) {
                            BKZpreprocess2(localB,&localU ,qmu,qc,qb,qB1,qlcindex,jj-qastart+1,kk-jj+1,BP.holdvecs,gmodprobs,sqrt(clim),persec,BP.preprocess_max_time,qasize,1,BP.verboselevel,0);
                        } else {
                            BKZpreprocess2(B,U ,qmu,qc,qb,qB1,qlcindex,jj,kk-jj+1,BP.holdvecs,gmodprobs,sqrt(clim),persec,BP.preprocess_max_time,qasize,1,BP.verboselevel,0);
                        }
                    }
                    #endif
                }
                if (BP.preprocessonly==true) break;
                
            } else {
                //When no preprocess
                bdelta = Estimate_delta(qc+jj-qastart+1,kk-jj+1);
            }

             if (process_flag==1) {
                    parallel=1; //Default is single thread
                    if (BP.multithread==true) {
                        ucost = to_RR(2*ENUMCost(qc,jj-qastart+1,kk-jj+1,gmodprobs,mod_init_radius,BP.init_mode));
                        //cout << "decideMT: " << ucost << " " << BP.MTlimit << endl;
                        //cout << jj << " " << kk << " " << mod_init_radius << " " << gmodprobs << endl;
                        if (ucost >BP.MTlimit)  {
                            parallel = BP.numthreads;
                        } // If the expected cost was too small, enumerate in single thread 
                    }
                    bool optimizepf=false;
                    if (BP.optimizepf==true) {
                        optimizepf = BP.optimizepf;
                        ucost = to_RR(2*ENUMCost(qc,jj-qastart+1,kk-jj+1,gmodprobs,mod_init_radius,BP.init_mode));
                        //cout << "opf: " << ucost << " " << BP.optimizepf_at_least << endl;
                        if (ucost < BP.optimizepf_at_least * 1000000.0) optimizepf = false;
                    }
                    //main part of BKZ (ENUM+update)
                    if (nolargeextract==0) {
                        if (BKZCore(localB,&localU,qmu,qc,qb,qB1,qlcindex,jj-qastart+1,kk-jj+1,BP.holdvecs,gmodprobs,mod_init_radius,BP.init_mode,enumlim[0],qasize,parallel,enum_mode_find_shortest,BP.verboselevel,to_double(bdelta),optimizepf)!=0) {
                            update=1;
                        }         
                    } else {
                        if (BKZCore(B,U,qmu,qc,qb,qB1,qlcindex,jj-qastart+1,kk-jj+1,BP.holdvecs,gmodprobs,mod_init_radius,BP.init_mode,enumlim[0],qasize,parallel,enum_mode_find_shortest,BP.verboselevel,to_double(bdelta),optimizepf)!=0) {
                            update=1;
                        }         
    #ifdef _inc_vectorpool_cpp
                        //For recvec
                        if (BP.recvec==true) {
                            for (i=0;i<m;i++) {
                                Vectorpool::Addvector(BP.recvecdata,B[i]);
                            }
                        }
    #endif
                    }
                    if ((BP.extend_blocksize==true) & (nofirst==1)) {
                        if (kk-jj+1==BP.beta[0]) {
                            if (tour_max_cost< current_processed * 1000000) {
                                tour_max_cost = current_processed * 1000000;
                                if (BP.verboselevel>=3) cout << "new_tour_max=" << tour_max_cost << endl;
                                BP.return_tour_max_cost = to_quad_float(tour_max_cost);
                            }
                        }
                    }
                    total_enum_nodes += current_processed * 1000000;
             } //end of processflag==1

            //end of tour?
            if ((jj == iend-1) || (jj==BP.abortbyindex)) { ;
                jj = istart-1;
                nofirst++;

                if (nofirst>BP.tourlim) break;
                if (update==0) break;   //no update during the current tour
                if (BP.abortbysec>0) {
                    if (BP.abortbysec <  ((double)clock() - startbkztime) / CLOCKS_PER_SEC) break;
                }
                update=0;
    /*
                if ((BP.abortbyfec>0) && (jj<=BP.startindex)) {
                    if (nolargeextract==0) {
                        cout << "logFEC=" << log(FullENUMCost(c,istart,iend-istart+1)) << " target_logFEC=" << BP.abortbyfec << endl; ;
                        if (log(FullENUMCost(c,istart,iend-istart+1))< BP.abortbyfec) break;
                    } else {
                        cout << "logFEC=" << log(FullENUMCost(qc+istart-1,iend-istart+1)) << " target_logFEC=" << BP.abortbyfec << endl; ;
                        if (log(FullENUMCost(qc+istart-1,iend-istart+1))< BP.abortbyfec) break;
                    }
                }
     */ 
            }   //end of "end of tour section"

        } //end of tour

        //final update
        //cout << "approx? : " << isapproxcomputed << endl; 
        if ((isapproxcomputed==1) && (nolargeextract==0)) {
            for (i=0;i<qasize;i++) {
                for (j=0;j<qasize;j++) {
                    if (j==0) {
                        localB2[i] = localU[i][j] * B[j+qastart-1];
                    } else {
                        localB2[i] += localU[i][j] * B[j+qastart-1];
                    }
                }
            }
            //cout << B(1) << endl;
            for (i=0;i<qasize;i++) {
                B[i+qastart-1] = localB2[i];
            }
        }

        if (BP.record_min_fec==true) {
            mat_ZZ LL;
            LL = BP.min_fec_basis;
            //cout << "bpcheck:" << BP.min_fec_basis[m-1] << endl;
            int k;
            k  = 0;
            for (i=0;i<BP.min_fec_basis.NumRows();i++) {
                if (LengthOf(BP.min_fec_basis[i])!=0) {
                    LL[k++] =  BP.min_fec_basis[i];
                }
            }
            LL.SetDims(k,BP.min_fec_basis.NumCols());
            BP.min_fec_basis = LL;
    /*
            cout << "bpcheck2:" << BP.min_fec_basis[m-1] << endl;
            cout << "bpcheck_min=" << LatticeTools::LatticeApprox(LL) << endl;;
            cout <<  BP.min_fec_basis.NumRows() << endl;
            cout <<  BP.min_fec_basis.NumCols() << endl;
     */ 
        }

    #ifdef _inc_vectorpool_cpp
        if (BP.recvec==true) {
            Vectorpool::distinct(BP.recvecdata);
        }
    #endif
        BB = B;
        BB.SetDims(m, n);

        verbose = lllverboseback; 

        return m;

    }
}


#endif
