
#include<iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
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

int main() 
{
    
   int m[4][4] = {
	    {15, -4, 1, 0},
		{ -4, 6, -4, 1},
		{ 1,-4, 6, -4},
		{ 0, 1,-4, 6}};
  vector<vector<float>> matrix(4, vector<float>(4) );
  for(int i =0; i<4; i++)
  {
	for(int j =0; j<4; j++)
	{
		matrix[i][j] = m[i][j];
	}
  }
    cout<<detPentMatrix(matrix,4)<<endl;

    return 0;
}