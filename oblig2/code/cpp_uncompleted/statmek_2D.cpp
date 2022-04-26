#include<iostream>
#include<vector>
#include<cstdlib>
#include<math.h>
#include<complex>


#define PI 3.14159265358979323
const int q=3;      // number of spin states
const int L=16;      // Linear system size 
const double T=2;   // Temperature in units of J 
const int J=4;

const int N=L*L;      // Tot. number of spins 
const double pconnect=1-exp(-J/T);  // Connection prob. 

const int NCLUSTERS=1;  // # of clusters in one step
const int NESTEPS=10000;// # of equil. MC steps
const int NMSTEPS=10000;// # of measurement MC steps
const int NBINS=10;     // # of measurement bins

std::vector<std::vector<int> > 
    S(L, std::vector<int>(L));  // Spin array 

std::vector<int> M(q);  // Number of spins in different states 
std::vector<std::complex<double> > W(q); // order param weights.

enum dirs{RIGHT,LEFT,UP,DOWN};
std::vector<int> indices(2);
int indx(int x){return x;}
int indy(int y){return y;}
int xpos(int i){return i%L;}
int ypos(int j){return j%L;}

int Nbr(int i, int j, int dir)
{
    int x=xpos(i);
    int y=ypos(j);
    switch(dir)
    {
        case RIGHT: indices[0] = indx((x+1)%L); 
                    indices[1] = indy(y%L);
                    break;
        case LEFT:  
                    indices[0] = indx((x-1+L)%L); 
                    indices[1] = indy(y%L);
                    break;
        case UP:    
                    indices[0] = indx(x%L); 
                    indices[1] = indy((y+1)%L);
                    break;
        case DOWN:  
                    indices[0] = indx(x%L); 
                    indices[1] = indy((y-1+L)%L);
                    break;
    }
    return 0;
}


void FlipandBuildFrom(int s, int v)
{
    int oldstate(S[s][v]),newstate((S[s][v]+1)%q);

    S[s][v]=newstate;
    M[oldstate]--; M[newstate]++;

    for(int dir=0; dir<4; dir++)
    {
        Nbr(s,v,dir);
        int xn=indices[0];
        int yn=indices[1];
        if(S[xn][yn]==oldstate)
            if (rand()/(RAND_MAX+1.) < pconnect){FlipandBuildFrom(xn,yn);}
    }
}

int main()
{   
    double m_avg=0;
    for(int s=0; s<q; s++)
        W[s]=std::complex<double>(cos(2*PI*s/q),sin(2*PI*s/q));
    for(int i=0; i<L; i++){
        for(int j=0; j<L; j++){
        S[i][j]=0;
        }
    } 
    
    for(int s=1; s<q; s++) M[s]=0;
    M[0]=N;
    srand((unsigned) time(0));

    for(int t=0; t<NESTEPS; t++)
        for(int c=0; c<NCLUSTERS; c++)
        {
            FlipandBuildFrom(rand()%L, rand()%L);
        }
    
    for(int n=0; n<NBINS; n++)
    {
        std::complex<double> m(0.,0.);
        double m1=0, m2=0, m4=0;

        for(int t=0; t<NMSTEPS; t++)
        {
            for(int c=0; c<NCLUSTERS; c++) FlipandBuildFrom(rand()%L, rand()%L);
            std::complex<double> tm(0., 0.);
            for(int s=0; s<q; s++){tm+=W[s]*double(M[s]);}
            tm/=N;
            double tm1=abs(tm);
            double tm2=tm1*tm1;
            m+=tm; m1+=tm1; m2+=tm2; m4+=tm2*tm2;
        }

        m/=NMSTEPS; m1/=NMSTEPS; m2/=NMSTEPS; m4/=NMSTEPS;
        m_avg += m1;
        // std::cout << m << " " << m1 << " " << m2 << " " << m4 << std::endl;
    }
    std::cout << "\n" << "m: " << m_avg/NBINS << std::endl;
}