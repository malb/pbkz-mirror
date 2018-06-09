#ifndef _inc_bkzconstants_hpp
#define _inc_bkzconstants_hpp

//Constant tables 

//Computing the volume of n-ball using boost library
template <typename T,typename T2> T vb(T2 dim,T radius) {
    T ret;
    ret = pow(boost::math::constants::pi<T>(),0.5*dim) / boost::math::tgamma<T>((T)dim*0.5 + 1.0) * pow(radius,dim);
    return ret;
}

template <typename T> T mypower(T a,double b) {
    return exp(log(a)*b);
}


template <typename T> class autoarray1d {
    std::vector<T> data;

    public:
    T& operator [] (int i) {
        if (i<0) {
            cout << "autoarray1d: range error" << endl;
            exit(0);
        }
        
        if (i>=data.size()) {
            data.resize(i+1);
        }
        return data[i];
    }
    int size() {return data.size();};
    int resize(int i) {data.resize(i);};

    void load(std::string& fname) {
        loadstdvector(data,fname);
    }

    void save(std::string& fname) {
        savestdvector(data,fname);
    }
    
};

template <typename T> class autoarray2d {
    autoarray1d<autoarray1d<T> > data;

    public:
    autoarray1d<T>& operator [] (int i) {
        if (i<0) {
            cout << "autoarray2d: range error" << endl;
            exit(0);
        }
        
        if (i>=data.size()) {
            data.resize(i+1);
        }
        return data[i];
    }
    int size() {return data.size();};
    int resize(int i) {data.resize(i);};

    void load(std::string& fname) {
        if (FileExists(fname)==false) {
            data.resize(0);
            return;
        }
        int fd = lock(fname);
        ifstream vstream;
        vstream.open(fname.c_str(),ios_base::in);
        if (vstream.fail()==true) {
            unlock(fd);
            return; 
        }
        int i,j;
        vstream >> i;
        data.resize(i);
        for (i=0;i<data.size();i++) {
            vstream >> j;
            data[i].resize(j);
            for (j=0;j<data[i].size();j++) {
                vstream >> data[i][j];
                //cout << i << " " << j << " " << data[i][j] << endl;
            }
        }
        vstream.close();
        unlock(fd);
    }

    void save(std::string& fname) {
        int fd;
        do {
                fd = lock(fname);
                usleep(100000);
        } while (fd<0);
        ofstream vstream;
        vstream.open(fname.c_str(),ios_base::trunc); 
        vstream << data.size() << endl;
        for (int i=0;i<(int)data.size();i++) {
            vstream << data[i].size() << "\t";
            for (int j=0;j<(int)data[i].size();j++) {
                vstream << data[i][j] << "\t";
            }
            vstream << endl;
        }
        vstream.close();
        unlock(fd);
    }
    
};

template <typename T> class autoarray3d {
    autoarray1d<autoarray1d<autoarray1d<T> > > data;

    public:
    autoarray1d<autoarray1d<T> >& operator [] (int i) {
        if (i<0) {
            cout << "autoarray3d: range error" << endl;
            exit(0);
        }
        
        if (i>=data.size()) {
            data.resize(i+1);
        }
        return data[i];
    }
    int size() {return data.size();};
    int resize(int i) {data.resize(i);};

    void load(std::string& fname,int checkdisp=0) {
        if (FileExists(fname)==false) {
            data.resize(0);
            return;
        }
        int fd = lock(fname);
        ifstream vstream;
        vstream.open(fname.c_str(),ios_base::in);
        if (vstream.fail()==true) {
            unlock(fd);
            return; 
        }
        int i,j,k;
        std::string tp;
        vstream >> tp;
        i = atoi(tp.c_str());
        if (checkdisp!=0) cout << "tp=" << tp << " i=" << i << endl;
        data.resize(i);
        for (i=0;i<data.size();i++) {
            vstream >> tp;
            j = atoi(tp.c_str());
            if (checkdisp!=0) cout << "tp=" << tp << " j=" << j << endl;
            data[i].resize(j);
            for (j=0;j<data[i].size();j++) {
                //vstream >> k;
                vstream >> tp;
                k = atoi(tp.c_str());
                if (checkdisp!=0) cout << "tp=" << tp << " k=" << k << endl;
                data[i][j].resize(k);
                for (k=0;k<data[i][j].size();k++) {
                    vstream >> tp;
                    data[i][j][k] = boost::lexical_cast<T>(tp);
                }
            }
        }
        vstream.close();
        unlock(fd);
    }

    void save(std::string& fname) {
        int fd;
        do {
                fd = lock(fname);
                usleep(100000);
        } while (fd<0);
        ofstream vstream;
        vstream.open(fname.c_str(),ios_base::trunc); 
        vstream << data.size() << endl;
        for (int i=0;i<(int)data.size();i++) {
            vstream << data[i].size() << "\t";
            for (int j=0;j<(int)data[i].size();j++) {
                vstream << data[i][j].size() << "\t";
                for (int k=0;k<(int)data[i][j].size();k++) {
                    vstream << data[i][j][k] << "\t";
                }
                vstream << endl;
            }
            vstream << endl;
        }
        vstream.close();
        unlock(fd);
    }
    
};

template <typename T> class autoarray4d {
    autoarray1d<autoarray1d<autoarray1d<autoarray1d<T> > > > data;

    public:
    autoarray1d<autoarray1d<autoarray1d<T> > >& operator [] (int i) {
        if (i<0) {
            cout << "autoarray4d: range error" << endl;
            exit(0);
        }
        
        if (i>=data.size()) {
            data.resize(i+1);
        }
        return data[i];
    }
    int size() {return data.size();};
    int resize(int i) {data.resize(i);};

    void load(std::string& fname) {
        if (FileExists(fname)==false) {
            data.resize(0);
            return;
        }
        int fd = lock(fname);
        ifstream vstream;
        vstream.open(fname.c_str(),ios_base::in);
        if (vstream.fail()==true) {
            unlock(fd);
            return; 
        }
        int i,j,k,m;
        vstream >> i;
        data.resize(i);
        for (i=0;i<data.size();i++) {
            vstream >> j;
            data[i].resize(j);
            for (j=0;j<data[i].size();j++) {
                vstream >> k;
                data[i][j].resize(k);
                for (k=0;k<data[i][j].size();k++) {
                    vstream >> m;
                    data[i][j][k].resize(m);
                    for (m=0;m<data[i][j][k].size();m++) {
                        vstream >> data[i][j][k][m];
                    }
                }
            }
        }
        vstream.close();
        unlock(fd);
    }

    void save(std::string& fname) {
        int fd;
        do {
                fd = lock(fname);
                usleep(100000);
        } while (fd<0);
        ofstream vstream;
        vstream.open(fname.c_str(),ios_base::trunc); 
        vstream << data.size() << endl;
        for (int i=0;i<(int)data.size();i++) {
            vstream << data[i].size() << "\t";
            for (int j=0;j<(int)data[i].size();j++) {
                vstream << data[i][j].size() << "\t";
                for (int k=0;k<(int)data[i][j].size();k++) {
                    vstream << data[i][j][k].size() << "\t";
                    for (int m=0;m<(int)data[i][j][k].size();m++) {
                        vstream << data[i][j][k][m] << "\t";
                    }
                }
                vstream << endl;
            }
            vstream << endl;
        }
        vstream.close();
        unlock(fd);
    }
    
};


template <typename T,typename T2> void vec_copy(std::vector<T>& a,autoarray1d<T2>& b) {
    a.resize(b.size());
    for (int i=0;i<b.size();i++) a[i] = b[i];
}

namespace bkzconstants{
    
    autoarray1d<bkzfloat> vol_unit_ball_table;
    autoarray1d<bkzfloat> ghconstant_table;
    autoarray1d<cpp_int> factorial_table;
    autoarray2d<bkzfloat> simlogfec_table;
    autoarray2d<cpp_int> binary_coeff_table;
 
    autoarray2d<bkzfloat> mainparam_table;
    autoarray3d<double> ealpha;
    autoarray3d<double> ualpha;
    autoarray3d<bkzfloat> eprob;
    
    autoarray3d<double> pfuppertable;
    autoarray3d<bkzfloat> pbkzsimlengthtable;

    autoarray4d<bkzfloat> simtime_enum;
    autoarray4d<double> simtime_lll;
    autoarray4d<double> simtime_loop;

    autoarray1d<bkzfloat> target_cost_lb_table;
    autoarray2d<bkzfloat> target_cost_lb_mult_table;
    autoarray1d<bkzfloat> target_cost_ub_table;

    bool initialized = false;
    std::string cachedir;    
    
    extern void loadcachetable();
    
    void initialize() {
        if (initialized == false) { 
            initialized = true;
            #ifdef _allow_cachefiles        
            cachedir =  makefullpath(ReadConf("bkz.conf","constantscache"));
            //cout << "bkzconst.cache: " << cachedir << endl;
            mkdirRecursive(cachedir.c_str(), 0777);
            loadcachetable();
            #endif
        }
    }

    
    bkzfloat vol_unit_ball(int dim) {

        initialize();
        if (dim<=0) return 0;
        if (dim==1) return 1;
        
        if (dim>100000) {
            cout << "vol_unit_ball: dimension error " << dim << endl;
            exit(0);
        }
        
        int cs = vol_unit_ball_table.size();
        if (dim >= cs) {
            #pragma omp critical        
            {
                vol_unit_ball_table.resize(dim+1);
                for (int i=cs;i<=dim;i++) vol_unit_ball_table[i] = 0;
            }
        }
        
        if (vol_unit_ball_table[dim]==0) {
            #pragma omp critical        
            {
                vol_unit_ball_table[dim] = vb<bkzfloat,bkzfloat>(dim,1.0);
            }
        }
        return vol_unit_ball_table[dim];
    }

    bkzfloat surface_unit_ball(int dim) {
        //Surface of n-ball
        return vol_unit_ball(dim) * dim;
    }

    bkzfloat ghconstant(int dim) {

        initialize();
        if (dim<=0) return 0;
        if (dim==1) return 1;
        
        if (dim>100000) {
            cout << "vol_unit_ball: dimension error " << dim << endl;
            exit(0);
        }
        
        int cs = ghconstant_table.size();
        if (dim >= cs) {
            #pragma omp critical        
            {
                ghconstant_table.resize(dim+1);
                for (int i=cs;i<=dim;i++) ghconstant_table[i] = 0;
            }
        }
        
        if (ghconstant_table[dim]==0) {
            #pragma omp critical        
            {
                ghconstant_table[dim] = pow(vb<bkzfloat,bkzfloat>(dim,1.0),-1.0/dim);
            }
        }
        return ghconstant_table[dim];
    }

cpp_int binary_coeff(int n,int k) {

        if (k<=0) return 1;
        if (k==1) return n;
        if (k>n) return 1;
        if (n==1) return 1;
        
        if (k > n/2) return binary_coeff(n,n-k);
        
        if ((n>100000) || (k>100000)) {
            cout << "binary_coeff error: n=" << n << " k=" << k << endl;
            exit(0);
        }
        
        int cs = binary_coeff_table.size();
        if (n >= cs) {
            #pragma omp critical        
            {
                binary_coeff_table.resize(n+1);
            }
        }
        
        cs = binary_coeff_table[n].size();
        if (k >= cs) {
            #pragma omp critical        
            {
                binary_coeff_table[n].resize(n/2+2);
            }
        }

        if (binary_coeff_table[n][k]==0) {
            #pragma omp critical        
            {
                binary_coeff_table[n][0] = 1;
                binary_coeff_table[n][1] = n;
                for (int kk=2;kk<=n/2;kk++) {
                    binary_coeff_table[n][kk] = binary_coeff_table[n][kk-1] * (n-kk+1)/kk;
                    //cout << "n=" << n << " k=" << kk << " " << binary_coeff_table[n][kk-1]  << " " << binary_coeff_table[n][kk] << endl;
                }
            }
        }
        
        return binary_coeff_table[n][k];
    }
    
    cpp_int factorial(int k) {

        initialize();
        if (k<=0) return 0;
        if (k==1) return 1;
        
        if (k>100000) {
            cout << "fractional: k error " << k << endl;
            exit(0);
        }
        
        if (factorial_table.size()<3) {
            factorial_table.resize(3);
            factorial_table[1] = 1;
            factorial_table[2] = 2;
        }
        int cs = factorial_table.size();
        if (k >= cs) {
            #pragma omp critical        
            {
                factorial_table.resize(k+1);
                for (int i=cs;i<=k;i++) factorial_table[i] = factorial_table[i-1] * i ;
            }
        }
        return factorial_table[k];
    }

    double last_load_time = -1;
    void loadcachetable() {
        std::string filename;
        initialize();
        #pragma omp critical        
        {
            if (last_load_time<0) {
                last_load_time = gettimeofday_sec(); 

                filename = cachedir + "/bkz.uvt.table";
                if (FileExists(filename)==true) {
                    vol_unit_ball_table.load(filename);
                }

                filename = cachedir + "/bkz.binary_coeff.table";
                if (FileExists(filename)==true) {
                    binary_coeff_table.load(filename);
                }

                filename = cachedir + "/bkz.factorial.table";
                if (FileExists(filename)==true) {
                    factorial_table.load(filename);
                }

                filename = cachedir + "/bkz.simlogfec.table";
                if (FileExists(filename)==true) {
                    simlogfec_table.load(filename);
                }

                filename = cachedir + "/bkz.mainparam";
                if (FileExists(filename)==true) {
                    mainparam_table.load(filename);
                }

                filename = cachedir + "/bkz.sim.lengths";
                if (FileExists(filename)==true) {
                    pbkzsimlengthtable.load(filename);
                }


                filename = cachedir + "/bkz.sim.ealpha";
                if (FileExists(filename)==true) {
                    ealpha.load(filename);
                }

                filename = cachedir + "/bkz.sim.ualpha";
                if (FileExists(filename)==true) {
                    ualpha.load(filename);
                }

                filename = cachedir + "/bkz.sim.eprob";
                if (FileExists(filename)==true) {
                    eprob.load(filename);
                }

                filename = cachedir + "/pf.upper.table";
                if (FileExists(filename)==true) {
                    pfuppertable.load(filename);
                }

                filename = cachedir + "/simtime.enum";
                if (FileExists(filename)==true) {
                    simtime_enum.load(filename);
                }

                filename = cachedir + "/simtime.lll";
                if (FileExists(filename)==true) {
                    simtime_lll.load(filename);
                }

                filename = cachedir + "/simtime.loop";
                if (FileExists(filename)==true) {
                    simtime_loop.load(filename);
                }

                filename = cachedir + "/simtime.cost.lbtable";
                if (FileExists(filename)==true) {
                    target_cost_lb_table.load(filename);
                }

                filename = cachedir + "/simtime.cost.ubtable";
                if (FileExists(filename)==true) {
                    target_cost_ub_table.load(filename);
                }

                filename = cachedir + "/simtime.cost.lbmulttable";
                if (FileExists(filename)==true) {
                    target_cost_lb_mult_table.load(filename);
                }

                filename = cachedir + "/ghconstant.table";
                if (FileExists(filename)==true) {
                    ghconstant_table.load(filename);
                }
            }
        }
    }

    double last_save_time = -1;
    
    void savecachetable(int force=0) {

        initialize();
        std::string filename;
        if ((gettimeofday_sec()-last_save_time<60) && (force==0)) return;
        last_save_time = gettimeofday_sec(); 

        #pragma omp critical
        {

            filename = cachedir + "/bkz.uvt.table";
            vol_unit_ball_table.save(filename);

            filename = cachedir + "/bkz.binary_coeff.table";
            binary_coeff_table.save(filename);

            filename = cachedir + "/bkz.factorial.table";
            factorial_table.save(filename);

            filename = cachedir + "/bkz.mainparam";
            mainparam_table.save(filename);

            filename = cachedir + "/bkz.simlogfec.table";
            simlogfec_table.save(filename);

            filename = cachedir + "/bkz.sim.ealpha";
            ealpha.save(filename);
            filename = cachedir + "/bkz.sim.ualpha";
            ualpha.save(filename);
            filename = cachedir + "/bkz.sim.eprob";
            eprob.save(filename);

            filename = cachedir + "/pf.upper.table";
            pfuppertable.save(filename);

            filename = cachedir + "/simtime.enum";
            simtime_enum.save(filename);

            filename = cachedir + "/simtime.lll";
            simtime_lll.save(filename);

            filename = cachedir + "/simtime.loop";
            simtime_loop.save(filename);

            filename = cachedir + "/simtime.cost.lbtable";
            target_cost_lb_table.save(filename);

            filename = cachedir + "/simtime.cost.ubtable";
            target_cost_ub_table.save(filename);

            filename = cachedir + "/simtime.cost.lbmulttable";
            target_cost_lb_mult_table.save(filename);

            filename = cachedir + "/ghconstant.table";
            ghconstant_table.save(filename);

            filename = cachedir + "/bkz.sim.lengths";
            pbkzsimlengthtable.save(filename);

        }
        
    }
}



#endif
