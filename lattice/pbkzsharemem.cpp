#ifndef _inc_pbkzsharemem
#define _inc_pbkzsharemem

namespace progressive_bkz {
    //shared memory
    int meminit=0;
    double** bkz_share_double;
    int** bkz_share_foundvec;
    long** bkz_share_long;
    quad_float** bkz_share_quad_float2;
    quad_float*** bkz_share_quad_float3;
    RR** bkz_share_RR2;
    RR*** bkz_share_RR3;
    
    clock_t start,end;  //for CPUtime
    double cputime;
    
    //ENUMspeed in Mnodes/sec
    //Not changed after first measuring
    double enext=2;
    double total_enum_nodes=0;
    
    double current_processed;    //in Mnodes
    
    //time
    double first_init_time=-1;
    double last_enum_time=-1;
    
    //used for caching bkz-simulator
    double** logfec;
    quad_float* temp_c;
    int current_exb=0;
    
    int bkz_verbose=0;  //outputlevel
    
    std::string simcachedir;
    
    void sharememalloc() {

        if (meminit!=0) return;
        
        int i;
        int numthread = omp_get_num_procs();

        simcachedir =  makefullpath(ReadConf("bkz.conf","simcache"));
        mkdir(simcachedir.c_str(),0777);
        

    
        bkz_share_double = new double*[15+numthread];
        bkz_share_long = new long*[10+numthread];
        bkz_share_foundvec = new int*[maxfoundvec];
        bkz_share_quad_float2 = new quad_float*[30];
        bkz_share_quad_float3 = new quad_float**[10];
	bkz_share_RR2 = new RR*[30];
        bkz_share_RR3 = new RR**[10];
        for (i=0;i<15+numthread;i++) bkz_share_double[i] = new double[latticemaxdim];
        for (i=0;i<10+numthread;i++) bkz_share_long[i] = new long[latticemaxdim];
        for (i=0;i<maxfoundvec;i++) bkz_share_foundvec[i] = new int[latticemaxdim];
        
        for (i=0;i<10;i++) {
            bkz_share_quad_float3[i] = new quad_float*[latticemaxdim];
	    bkz_share_RR3[i] = new RR*[latticemaxdim];
            //bkz_share_RR3[i].SetLength(latticemaxdim);
            for (int j=0;j<latticemaxdim;j++) {
                bkz_share_quad_float3[i][j] = new quad_float[latticemaxdim];
		// bkz_share_RR3[i][j].SetLength(latticemaxdim);
		bkz_share_RR3[i][j] = new RR[latticemaxdim];
            }

        }
        for (i=0;i<30;i++) {
            bkz_share_quad_float2[i] = new quad_float[latticemaxdim];
            bkz_share_RR2[i] = new RR[latticemaxdim];
        }
        cputime = 0;
        
        if (first_init_time==-1) first_init_time=clock();
        if (last_enum_time==-1) last_enum_time=clock();

        //simulator
        logfec = new double*[latticemaxdim+1];
        for (int i=0;i<latticemaxdim+1;i++) {
            logfec[i] = new double[i+1];
            for (int j=0;j<=i;j++) logfec[i][j] = -1;
        }
        temp_c = new quad_float[latticemaxdim+1];

        meminit = 1;
    }
    
}




#endif