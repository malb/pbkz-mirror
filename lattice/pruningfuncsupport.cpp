#ifndef _inc_pruning_func_support_cpp
#define _inc_pruning_func_support_cpp

#include "computeevensimplex.cpp"
#include "samplingtools.cpp"

namespace pruning_func {

    double Rigid_lower_prob(double* R,int n) {
        EvenSimplex::Initialize();

        int i;
        double *R2;
        R2 = EvenSimplex::es_share_double;

        if (n%2==0) {
            for (i=1;i<=n/2;i++) R2[i]=R[2*i-1];

            //cout << "qf=" << EvenSimplex::EvenSimplexVolumeqf(n/2,R2,opt_surface_prob);
            //cout << " RR=" << EvenSimplex::EvenSimplexVolume(n/2,R2,opt_surface_prob) << endl;

            return to_double(EvenSimplex::EvenSimplexVolume(n/2,R2,opt_surface_prob));
        } else {
            for (i=1;i<=(n-1)/2;i++) R2[i]=R[2*i-1];
            return to_double(EvenSimplex::EvenSimplexVolume((n-1)/2,R2,opt_surface_prob)); 
        }    
    }

    double Rigid_upper_prob(double* R,int n) {
        EvenSimplex::Initialize();
        int i;
        double *R2;
        R2 = EvenSimplex::es_share_double;
        if (n%2==0) {
            for (i=1;i<=n/2;i++) R2[i]=R[2*i];
            return to_double(EvenSimplex::EvenSimplexVolume(n/2,R2,opt_surface_prob));
        } else {
            for (i=1;i<=(n-1)/2;i++) R2[i]=R[2*i];
            R2[(n+1)/2] = 1;    
            return to_double(EvenSimplex::EvenSimplexVolume((n+1)/2,R2,opt_surface_prob));  //upper bound
        }
    }

    RR Rigid_lower_cost(double* R,quad_float* c,int n,double radius_c) {
        EvenSimplex::Initialize();
        RR sum;
        RR t,t2;
        int i,k;
        double* R2;
        R2 = EvenSimplex::es_share_double;
        t2 = 1;
        for (k=2;k<=n;k+=2) {
          if (R[k-1]!=0) {
            for (i=1;i<=k/2;i++) R2[i] = R[2*i-1]/R[k-1];
            t = EvenSimplex::EvenSimplexVolume(k/2,R2,opt_volume_prob);
            t *= VolumeBall(k,sqrt(R[k-1]));
            t2 *= radius_c * radius_c / to_RR(sqrt(c[n-k+1]*c[n-k+2]));
            //cout << "Factor(" << k << ")=" << t*t2 << endl;
            sum += t * t2;
          }
        }

        if (sum<0) {
            cout << "lower" << endl;
            for (i=1;i<=n;i++) {
                cout << R[i] << endl;
            }
            t2 = 1;
            for (k=2;k<=n;k+=2) {
                if (R[k-1]!=0) {
                    for (i=1;i<=k/2;i++) R2[i] = R[2*i-1]/R[k-1];
                    t = EvenSimplex::EvenSimplexVolume(k/2,R2,opt_volume_prob);
                    cout << t << endl;
                    t *= VolumeBall(k,sqrt(R[k-1]));
                    cout << t << endl;
                    t2 *= radius_c * radius_c / to_RR(sqrt(c[n-k+1]*c[n-k+2]));
                    cout << t2 << endl;
                    cout << "Factor(" << k << ")=" << t*t2 << endl;
                    sum += t * t2;
                }
            }
            exit(0);

        }
        return sum;
    }

    RR Rigid_lower_cost_double(double* R,double* c,int n,double radius_c) {
        EvenSimplex::Initialize();
        RR sum;
        RR t,t2;
        int i,k;
        double* R2;
        R2 = EvenSimplex::es_share_double;
        t2 = 1;
        for (k=2;k<=n;k+=2) {
          if (R[k-1]!=0) {
            for (i=1;i<=k/2;i++) R2[i] = R[2*i-1]/R[k-1];
            R2[i] = 1.0;
            t = EvenSimplex::EvenSimplexVolume(k/2,R2,opt_volume_prob);
            t *= VolumeBall(k,sqrt(R[k-1]));
            t2 *= radius_c * radius_c / sqrt(c[n-k+1]) / sqrt(c[n-k+2]);
            //cout << "Factor(" << k << ")=" << t*t2 << endl;
            sum += t * t2;
          }
        }
        return sum;
    }


    RR Rigid_upper_cost(double* R,quad_float* c,int n,double radius_c) {
        //pruning_func[1,...,n]
        //c[1,...,n]
        EvenSimplex::Initialize();
        RR sum;
        RR t,t2;
        int i,k;
        double* R2;
        R2 = EvenSimplex::es_share_double;
        t2 = 1;
        sum = to_RR(0);
        for (k=2;k<=n;k+=2) {
          if (R[k]!=0) {
           for (i=1;i<=k/2;i++) {
               R2[i] = R[2*i]/R[k];
           }
            t = EvenSimplex::EvenSimplexVolume(k/2,R2,opt_volume_prob);
            //cout << k << " " << t << " ";
            t *= VolumeBall(k,sqrt(R[k]));
            t2 *= radius_c * radius_c / to_RR(sqrt(c[n-k+1]*c[n-k+2]));
            //cout << t << " " << t2 << " ";
            sum += t * t2;

            if ((k>n/2) && (t*t2 < sum*0.0001)) break;  //output the approximated value
            //cout << sum << endl;
          }
        }
        

        if (sum<0) {
            //bug??
            cout << "upper" << endl;
            for (i=1;i<=n;i++) {
                cout << R[i] << endl;
            }
            t2 = 1;
            sum = to_RR(0);
            for (k=2;k<=n;k+=2) {
                if (R[k]!=0) {
                    for (i=1;i<=k/2;i++) R2[i] = R[2*i]/R[k];
                    t = EvenSimplex::EvenSimplexVolume(k/2,R2,opt_volume_prob);
                    t *= VolumeBall(k,sqrt(R[k]));
                    t2 *= radius_c * radius_c / to_RR(sqrt(c[n-k+1]*c[n-k+2]));
                    cout << "k=" << k << " t=" << t << " t2=" << t2 << endl;
                    sum += t * t2;
                    cout << sum << endl;
                }
            }
            exit(0);

        }

        return sum; //for SVP, enumeration cost is halved by symmetry
    }

    quad_float Rigid_upper_cost_pcvf(quad_float* V,quad_float* c,int n,double radius_c,quad_float* factors=0) {
        //pruning_func[1,...,n]
        //c[1,...,n]
        //pcvf=pre-computed volume factor
        //radius is not squared
        quad_float sum;
        quad_float t,t2;
        int k;
        t2 = 1;
        for (k=2;k<n;k+=2) {
            t2 *= radius_c * radius_c / sqrt(to_double(c[n-k+1])*to_double(c[n-k+2]));
            //cout << "Factor-pcvf(" << k << ")=" << V[k/2]*t2 << endl;
            if (factors==0) {
                sum += V[k/2] * t2;
            } else {
                factors[k/2] = V[k/2] * t2;
                sum += factors[k/2];
            }
        }
        return sum;
    }

    RR Rigid_upper_cost_double(double* R,double* c,int n,double radius_c) {
        EvenSimplex::Initialize();
        RR sum;
        RR t,t2;
        int i,k;
        double* R2;
        R2 = EvenSimplex::es_share_double;
        t2 = 1;
        for (k=2;k<=n;k+=2) {
          if (R[k]!=0) {
           for (i=1;i<=k/2;i++) R2[i] = R[2*i]/R[k];
            R2[i] = 1.0;
            t = EvenSimplex::EvenSimplexVolume(k/2,R2,opt_volume_prob);
            t *= VolumeBall(k,sqrt(R[k]));
            t2 *= radius_c * radius_c / sqrt(c[n-k+1]*c[n-k+2]);
            sum += t * t2;
          }
        }
        return sum;
    }

    double Approx_prob(double* pf,int n) {
        //pf[1...n]
        if (n%2==0) {
            double ret = Rigid_upper_prob(pf,n) * sampling_tools::Ratio_prob(pf,n); 
            //cout << "ap_ret=" << ret << endl; 
            return ret;
        } else {
            pf[n+1]=1.0;
            pf[n+2]=1.0;
            return Rigid_upper_prob(pf,n) * sampling_tools::Ratio_prob_odd(pf,n);
        }
    }

}


#endif
