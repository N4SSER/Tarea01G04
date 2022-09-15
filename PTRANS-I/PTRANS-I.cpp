#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include<ctime>
using namespace std;
float detPentMatrix(vector<vector<float>> a, int n)
{

	float p1 = a[0][0];
	float p2 = a[1][1]*p1 - a[1][0]*a[0][1];
	float r2 = a[0][1];
	float s2 = a[1][0];
	float p3 = a[2][2] * p2 - a[2][1] * a[1][2] * p1 - a[2][0] * a[0][2] * a[1][1] + a[2][0] * a[1][2] * r2 + a[2][1] * a[0][2] * s2;
	float r3 = a[1][2] * p1 - a[0][2] * s2;
	float s3 = a[2][1] * p1 - a[2][0] * r2;
	float p4 = a[3][3] * p3 - a[3][2] * a[2][3] * p2 - a[3][1] * a[1][3] * (a[2][2] * p1 - a[2][0] * a[0][2]) + a[3][1] * a[2][3] * r3 + a[3][2] * a[1][3] * s3;

	vector<float> p;
	p.push_back(p1);
	p.push_back(p2);
	p.push_back(p3);
	p.push_back(p4);

	vector<float> r;
	r.push_back(1);
	r.push_back(r2);
	r.push_back(r3);

	vector<float> s;
	s.push_back(1);
	s.push_back(s2);
	s.push_back(s3);



	for (int i = 4; i < n; i++)
	{
		float newR = a[i - 2][i - 1] * (float)p[i - 3] - a[i - 3][i - 1] * (float)s[i - 2];
		r.push_back(newR);
		float newS = a[i - 1][i - 2] * (float)p[i - 3] - a[i - 1][i - 3] * (float)r[i - 2];
		s.push_back(newS);
		float newP = a[i][i] * (float)p[i - 1] - a[i][i - 1] * a[i - 1][i] * (float)p[i - 2] - a[i][i - 2] * a[i - 2][i] * (a[i - 1][i - 1] * (float)p[i - 3] - a[i - 1][i - 3] * a[i - 3][i - 1] * (float)p[i - 4]) + a[i][i - 2] * a[i - 1][i] * (float)r[i - 1] + a[i][i - 1] * a[i - 2][i] * (float)s[i - 1];
		p.push_back(newP);

	}

	return p.back();
}

float* PTRANS_I( float* a, float* b, float* c, float* d, float* e, int n,  float* y)
{
        float* mu =  (float*)malloc((n)*sizeof(float));
        float* alpha = (float*)malloc((n-1)*sizeof(float));
        float* beta = (float*)malloc((n-2)*sizeof(float));
        float* zeta = (float*)malloc((n)*sizeof(float));
        float* gamma = (float*)malloc((n)*sizeof(float));
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
        zeta[i] = (y[i] - zeta[i-2]*e[i] - zeta[i-1]*gamma[i])/mu[i];
    }
    y[n-2] = c[n-2] - alpha[n-4]*e[n-2];
    mu[n-2] = d[n-2] - beta[n-4]*e[n-2] - alpha[n-3]*gamma[n-2];
    alpha[n-2] = (a[n-2] -beta[n-3]*gamma[n-2])/mu[n-2];
    gamma[n-1] = c[n-1] - alpha[n-3]*e[n-1];
    mu[n-1] = d[n-1] -beta[n-3]*e[n-1] - alpha[n-2]*gamma[n-1];
    zeta[n-2] = (y[n-2] - zeta[n-3]*e[n-2] - zeta[n-3]*gamma[n-2])/mu[n-2];
    zeta[n-1] = (y[n-1]-zeta[n-2]*e[n-1] -zeta[n-2]*gamma[n-1])/mu[n-1];
    
    //6
    float* x = (float*)malloc(n*sizeof(float));
    x[n-1] = zeta[n-1]; x[n-2] = zeta[n-2] -alpha[n-2]*x[n-1];
    int i = n-3;
    while(i>-1)
    {
        x[i] = zeta[i] - alpha[i]*x[i+1] - beta[i]*x[i+2];
        i--;
    }
    return x;
}

float* PTRANS(vector<vector<float>> m, int n, float* y)
{
  
    if(detPentMatrix(m,n) == 0)
    {
        cout<<"No es invertible";
        return nullptr;
        
    }
        
    float* a =  (float*)malloc((n)*sizeof(float));
    float* b = (float*)malloc((n)*sizeof(float));
    float* c = (float*)malloc((n)*sizeof(float));
    float* d = (float*)malloc((n)*sizeof(float));
    float* e = (float*)malloc((n)*sizeof(float));
  
      for(int i =0; i<n; i++)
    {
        a[i] = 0;
        b[i] = 0;
        c[i] = 0;
        d[i] = 0;
        e[i] = 0;
    }
    for(int i = 0; i<n;i++)
    {
        d[i] = m[i][i];
        if(i-1>=0)
        {
            a[i] = m[i-1][i];
            c[i+1] = m[i][i-1];
            
        }
         if(i-2>=0)
        {
            b[i] = m[i-2][i];
            e[i+2] = m[i][i-2];
        }

       
        
    }
    return PTRANS_I(a,b,c,d,e,4,y);
}

vector<vector<float>> mat(int n)
{
     vector<vector<float>> m(n, vector<float>(n) );
     for(int i =0; i<n;i++)
    {
        for(int j =0; j<n;j++)
            m[i][j] = 0;
    }
    for(int i =0; i<n;i++)
    {
         for(int j =0; j<n;j++)
         {
            m[i][i] = 6;
            if(i-1>=0)
            {
                m[i-1][i] = -4;
                m[i][i-1] = -4;;
            }
         if(i-2>=0)
            {
                m[i-2][i] = 1;
                m[i][i-2] = 1;;
            }
         }
    }
    m[0][0] =9;
    m[n-2][n-2] = 5;
    m[n-1][n-1] = -2;
    m[n-2][n-1] = -2;
    return m;
}
float* result(int n)
{
    float* y = (float*)malloc(n*sizeof(float));
    for(int i =0; i<n;i++)
        y[i] = 1;
    return y;
}
int main(int argc, char *argv[]) 
{
    std::clock_t start;
    double duration;

    start = std::clock();

    int n=atoi(argv[1]);
    vector<vector<float>> matrix = mat(n);
  

     float *x = PTRANS(matrix,n,result(n));



    duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;

    std::cout<<"n =  "<<n<<": "<< duration <<'s';
    return 0;
}