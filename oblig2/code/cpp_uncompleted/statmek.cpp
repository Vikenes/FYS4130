#include<iostream>
#include<vector>
#include<cstdlib>
#include<math.h>
#include<complex>
#include<fstream>


#define PI 3.14159265358979323
const int q=3;      // number of spin states
const int L=16;      // Linear system size 
const double T=1;   // Temperature in units of J 
const int J=4;

const int N=L;      // Tot. number of spins 
const double pconnect=1-exp(-J/T);  // Connection prob. 

const int NCLUSTERS=1;  // # of clusters in one step
const int NESTEPS=10000;// # of equil. MC steps
const int NMSTEPS=10000;// # of measurement MC steps
const int NBINS=10;     // # of measurement bins

std::vector<int> S(N);  // Spin array 
std::vector<int> M(q);  // Number of spins in different states 
std::vector<std::complex<double> > W(q); // order param weights.

std::complex<double> m_0(0.,0.);
std::vector<std::complex<double>> m_0r(L);
std::vector<std::complex<double>> m_r(L);


enum dirs{RIGHT,LEFT};
int indx(int x){return x;}
int xpos(int i){return i%L;}

int Nbr(int i,int dir)
{
    int x=xpos(i);
    switch(dir)
    {
        case RIGHT: return indx((x+1)%L);
        case LEFT:  return indx((x-1+L)%L);
    }
    return 0;
}


void FlipandBuildFrom(int s)
{
    int oldstate(S[s]),newstate((S[s]+1)%q);

    S[s]=newstate;
    M[oldstate]--; M[newstate]++;

    for(int dir=0; dir<2; dir++)
    {
        int j=Nbr(s,dir);
        if(S[j]==oldstate)
            if (rand()/(RAND_MAX+1.) < pconnect){FlipandBuildFrom(j);}
    }
}

int main()
{   
    double m_avg=0;
    for(int s=0; s<q; s++){
        W[s]=std::complex<double>(cos(2*PI*s/q),sin(2*PI*s/q));
    }
    for(int i=0; i<N; i++) S[i]=0;
    for(int s=1; s<q; s++) M[s]=0;
    M[0]=N;
    srand((unsigned) time(0));

    for(int t=0; t<NESTEPS; t++)
        for(int c=0; c<NCLUSTERS; c++)
        {
            FlipandBuildFrom(rand()%N);
        }
    
    for(int n=0; n<NBINS; n++)
    {
        std::complex<double> m(0.,0.);
        double m1=0, m2=0, m4=0;

        for(int t=0; t<NMSTEPS; t++)
        {
            for(int c=0; c<NCLUSTERS; c++) FlipandBuildFrom(rand()%N);

            m_0 += std::complex<double> (cos(2*PI*S[0]/q),-sin(2*PI*S[0]/q));
            for(int r=0; r<L; r++){
                m_r[r] += std::complex<double>(cos(2*PI*S[r]/q),sin(2*PI*S[r]/q));
                m_0r[r]+= std::complex<double>(cos(2*PI*(S[r]-S[0])/q),sin(2*PI*(S[r]-S[0])/q));
            }


            std::complex<double> tm(0., 0.);
            for(int s=0; s<q; s++){
                tm+=W[s]*double(M[s]);
                }
            tm/=N;
            double tm1=abs(tm);
            double tm2=tm1*tm1;
            m+=tm; m1+=tm1; m2+=tm2; m4+=tm2*tm2;
        }

        m/=NMSTEPS; m1/=NMSTEPS; m2/=NMSTEPS; m4/=NMSTEPS;
        m_avg += m1;

    }
    std::cout << "\n" << "m: " << m_avg/NBINS << std::endl;

    std::cout << "huh" << std::endl; 
    
    std::ofstream corr_values ("corr_params_TJ_025.txt");
    corr_values << "m0r_re, m0r_im, m0mr_re, m0mr_im" << std::endl;

    m_0/=(NMSTEPS*NBINS);
    for(int r=0; r<L; r++){
        m_r[r] /= (NMSTEPS*NBINS);
        m_0r[r]/= (NMSTEPS*NBINS);
        corr_values << real(m_0r[r]) << ", " << imag(m_0r[r]) << ", ";
        corr_values << real(m_0*m_r[r]) << ", " << imag(m_0*m_r[r]) << std::endl;
    }
    corr_values.close();
    // std::cout << m_r[2] << std::endl; //<< " " << m1 << " " << m2 << " " << m4 << std::endl;
}