#include <stdio.h>
#include <stdlib.h>
float* PTRANS( float* a, float* b, float* c, float* d, float* e, int n,  float* y)
{
    float* mu = malloc((n)*sizeof(float));
    float* alpha = malloc((n-1)*sizeof(float));
    float* beta = malloc((n-2)*sizeof(float));
    float* zeta = malloc((n)*sizeof(float));
    float* gamma = malloc((n)*sizeof(float));
    //3
    mu[0] = d[0]; alpha[0]=(a[0]/mu[0]);
    beta[0]=(b[0]/mu[0]); zeta[0]=(y[0]/mu[0]);
    //4
    gamma[1] = c[1]; mu[1] = d[1] -alpha[0]*gamma[1];
    alpha[1] = (a[1] -beta[0]*gamma[1])/mu[1];
    beta[1] = b[1]/mu[1]; zeta[1] = (y[1]-zeta[0]*gamma[1])/mu[1];
    //5
    for(int i = 2; i<n-2;i++)
    {
        gamma[i] = c[i] - alpha[i-2]*e[i];
        mu[i] = d[i] - b[i-2]*e[i] -alpha[i-1]*gamma[i];
        alpha[i] = (a[i]-beta[i-1]*gamma[i])/mu[i];
        beta[i] = b[i]/mu[i];
        zeta[i] = (y[i] - zeta[i-2]*e[i] - zeta[i-1]*gamma[i]);
    }
    y[n-2] = c[n-2] - alpha[n-4]*e[n-2];
    mu[n-2] = d[n-2] - beta[n-4]*e[n-2] - alpha[n-3]*gamma[n-2];
    alpha[n-2] = (a[n-2] -beta[n-3]*gamma[n-2])/mu[n-2];
    gamma[n-1] = c[n-1] - alpha[n-3]*e[n-1];
    mu[n-1] = d[n-1] -beta[n-3]*e[n-1] - alpha[n-2]*gamma[n-1];
    zeta[n-2] = (y[n-2] - zeta[n-3]*e[n-2] - zeta[n-3]*gamma[n-2])/mu[n-2];
    zeta[n-1] = (y[n-1]-zeta[n-2]*e[n-1] -zeta[n-2]*gamma[n-1])/mu[n-1];
    
    //6
    float* x = malloc(n*sizeof(float));
    x[n-1] = zeta[n-1]; x[n-2] = zeta[n-2] -alpha[n-2]*x[n-1];
    int i = n-3;
    while(i>-1)
    {
        x[i] = zeta[i] - alpha[i]*x[i+1] - beta[i]*x[i+2];
        i--;
    }
    return x;

}
int main() 
{
    float d[]= {15,12,19,21,11,2,4};
    float a[]= {-8,-4,-9,6,8,4};
    float b[]= {-6,-4,4,7,3};
    float c[]= {-2,-4,-9,10,-2,2};
    float e[]= {-6,-1,9,10,-2};
    float y[]= {300,0,0,0,1,2,6};
    float *x = PTRANS(a,b,c,d,e,7,y);
    for(int i = 0; i<7; i++)
    {
        printf("%f\n",x[i]);
    }
    return 0;
}