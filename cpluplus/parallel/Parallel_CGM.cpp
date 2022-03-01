#include<iostream>
#include <iomanip>
#include<cmath>
#include "Parallel_CGM.h"
using namespace std;

// c = a'*b
float_t Parallel_INNER_PRODUCT(vector<float_t> &a,vector<float_t>&b){
    float_t c = 0;
    int_t N = a.size();
#pragma omp parallel for reduction(+:c)
    for(int_t i=0;i<N;i++){
        c+=a[i]*b[i];
    }
    return c;
}

// r = b-A*x
vector<float_t>  Parallel_Residual(vector<vector<float_t>>&A,vector<float_t> &x,vector<float_t>& b){
    float_t temp = 0;
    int_t N = A.size();
    vector<float_t> r(N,0);
#pragma omp parallel for private(temp)
    for(int_t i=0;i<N;i++){
        temp = 0;
#pragma omp parallel for reduction(+:temp)
        for(int_t j=0;j<N;j++){
            temp += A[i][j]*x[j];
        }
        r[i] = -temp + b[i];
    }
    return r;
}

// c=A*b
vector<float_t>  Parallel_MATRIX_VECTOR_PRODUCT(vector<vector<float_t>>&A, vector<float_t> &x){
	int_t N = A.size();
	int_t M = A[0].size();
    vector<float_t> c(N,0);
#pragma omp parallel for
    for(int_t i=0;i<N;i++){
        float_t c_temp=0;
#pragma omp parallel for reduction(+:c_temp)
        for(int_t j=0;j<M;j++){
        	c_temp += A[i][j]*x[j];
        }
        c[i] = c_temp;
    }
    return c;
}

// vector update
vector<float_t> Parallel_saxpy(float_t a, vector<float_t> &x, vector<float_t> &y){
	int_t N = x.size();
	vector<float_t> r(N,0);
#pragma omp parallel for
	for(int_t j=0;j<N;j++){
		r[j] = a * x[j] + y[j];
	}
	return r;
}
/*
vector<float_t>  MATRIX_VECTOR_PRODUCT(vector<vector<float_t>>*A, vector<float_t> &x){
	int_t N = A->size();
	int_t M = A->at(0).size();
    vector<float_t> c(N,0);
    for(int_t i=0;i<N;i++){
        for(int_t j=0;j<M;j++){
        	c[i] += A->at(i).at(j)*x[j];
        }
    }
    return c;
}

// C=A*B
//template <class T>
vector<vector<float_t>> MATRIX_MATRIX_MULTIPLY(vector<vector<float_t>> &a, vector<vector<float_t>> &b)
{
    const int_t n = a.size();
    const int_t m = a[0].size();
    const int_t p = b[0].size();

    vector<vector<float_t>> c(n, vector<float_t>(p, 0));
    for (auto j = 0; j < p; ++j)
    {
        for (auto k = 0; k < m; ++k)
        {
            for (auto i = 0; i < n; ++i)
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}
*/


vector<float_t> Parallel_CG(vector<vector<float_t>>& A,vector<float_t>& b){
	int_t N = A.size();
	vector<float_t> r(N,-1);
	vector<float_t> p(N,0);
	vector<float_t> q(N,0);
	vector<float_t> x(N,0);

    float_t norm_b = sqrt(Parallel_INNER_PRODUCT(b,b));
    int count;
    float_t alpha, beta, rho1, rho, error;
    alpha = beta = rho1 = rho = error = 1;
    r = Parallel_Residual(A,x,b);
    for(count = 0;count<N;count++){
    
        rho = Parallel_INNER_PRODUCT(r,r);
        beta = rho/rho1;
        error = sqrt(rho)/norm_b;
		//cout << error << endl;
		if(error < TOL){
			break;
		}

        p= Parallel_saxpy(beta, p, r);
		q = Parallel_MATRIX_VECTOR_PRODUCT(A,p);

		alpha = rho/ Parallel_INNER_PRODUCT(p,q);
		x= Parallel_saxpy(alpha,p,x);
		r= Parallel_saxpy(-alpha,q,r);

        rho1 = rho;
    }
    //cout << "Serial CG reaches convergence at " << count << " iterations." << endl;
    return x;
}
